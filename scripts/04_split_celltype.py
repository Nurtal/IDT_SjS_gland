"""
04_split_celltype.py — Phase 1.3 : extrait les 10 modules cell-type de la
SjD Map à partir de `node_to_celltype.tsv`.

Pour chaque cell-type ∈ {SGEC, TH1, TH17, TFH, TREG, BCELL, PLASMA, M1, M2, PDC},
écrit dans `01_disease_map/celltype_modules/<celltype>/` :

    <celltype>_nodes.tsv       — manifeste des nœuds (core + ports + phenotype)
    <celltype>_reactions.tsv   — réactions strictement contenues
    <celltype>_metrics.json    — métriques topologiques

Et au niveau parent :

    celltype_modules_summary.json  — agrégé tous cell-types
    celltype_modules_report.md     — rapport humain + Gate 1.3

Usage :
    python scripts/04_split_celltype.py [--cache-dir 01_disease_map/cache]
                                          [--n2c 01_disease_map/node_to_celltype.tsv]
                                          [--output-dir 01_disease_map/celltype_modules]
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.celltype_module import (  # noqa: E402
    CELL_TYPES,
    CellTypeModule,
    GateResult,
    collect_extra_and_phenotype,
    evaluate_gate,
    extract_module,
    index_by_celltype,
    load_node_to_celltype,
    write_module_nodes_tsv,
    write_module_reactions_tsv,
)
from lib.map_audit import build_graph  # noqa: E402
from lib.minerva_api import (  # noqa: E402
    load_or_fetch_elements,
    load_or_fetch_reactions,
    make_session,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Rapport markdown
# ---------------------------------------------------------------------------


def render_report(
    summary: dict[str, Any],
    gates: dict[str, GateResult],
) -> str:
    L: list[str] = []
    L.append("# Cell-type Modules Report — Phase 1.3")
    L.append("")
    L.append("> Sortie de `scripts/04_split_celltype.py`. Extraction des sous-cartes "
             "à partir de `node_to_celltype.tsv` (Phase 1.2bis raffiné R6c).")
    L.append("")

    # Synthèse globale
    L.append("## 1. Synthèse globale")
    L.append("")
    L.append(f"- Cell-types extraits : **{len(CELL_TYPES)}**")
    L.append(f"- Réactions sources    : **{summary['n_reactions_total']}**")
    L.append(f"- Nœuds sources        : **{summary['n_elements_total']}**")
    L.append(f"- Modules valides Gate 1.3 : "
             f"**{summary['n_modules_pass']}/{len(CELL_TYPES)}**")
    L.append("")

    # Tableau par module
    L.append("## 2. Métriques par module")
    L.append("")
    L.append("| Cell-type | n_core | n_ports | n_pheno | n_total | n_reactions | "
             "n_edges | n_loops | density |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for ct in CELL_TYPES:
        m = summary["per_celltype"][ct]
        L.append(
            f"| {ct} | {m['n_core_nodes']} | {m['n_extra_ports']} | "
            f"{m['n_phenotype_outputs']} | {m['n_total_nodes']} | "
            f"{m['n_reactions']} | {m['n_edges']} | {m['n_feedback_loops']} | "
            f"{m['density']} |"
        )
    L.append("")

    # Gate 1.3
    L.append("## 3. Gate 1.3")
    L.append("")
    L.append("Critères : ≥80 nœuds, ≥50 réactions, ≥1 boucle de feedback, "
             "≥1 phénotype connecté.")
    L.append("")
    L.append("| Cell-type | ≥80 nœuds | ≥50 réactions | ≥1 loop | ≥1 phénotype | Statut |")
    L.append("|---|---|---|---|---|---|")
    for ct in CELL_TYPES:
        g = gates[ct]
        marks = [
            "✅" if g.pass_min_nodes else "❌",
            "✅" if g.pass_min_reactions else "❌",
            "✅" if g.pass_feedback else "❌",
            "✅" if g.pass_phenotype else "❌",
        ]
        status = "**PASS**" if g.pass_all else "**FAIL**"
        L.append(f"| {ct} | {marks[0]} | {marks[1]} | {marks[2]} | {marks[3]} | {status} |")
    L.append("")

    # Recommandations
    L.append("## 4. Recommandations")
    L.append("")
    failed = [ct for ct in CELL_TYPES if not gates[ct].pass_all]
    if not failed:
        L.append("- ✅ Tous les modules passent Gate 1.3 — passage Phase 1.4 autorisé.")
    else:
        L.append("- Modules à enrichir avant Phase 1.4 :")
        for ct in failed:
            g = gates[ct]
            issues = []
            if not g.pass_min_nodes:
                issues.append(f"<80 nœuds ({g.n_nodes})")
            if not g.pass_min_reactions:
                issues.append(f"<50 réactions ({g.n_reactions})")
            if not g.pass_feedback:
                issues.append("aucune boucle feedback détectée")
            if not g.pass_phenotype:
                issues.append("aucun phénotype connecté")
            L.append(f"    - **{ct}** : {', '.join(issues)} → curation manuelle "
                     f"(littérature + CellDesigner)")
    L.append("")
    L.append("- Revue expert SjD signée requise sur `<celltype>_nodes.tsv` "
             "avant Phase 1.4.")
    L.append("")

    L.append("---")
    L.append("")
    L.append("_Généré par `scripts/04_split_celltype.py`_")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-dir", type=Path,
                   default=ROOT / "01_disease_map" / "cache")
    p.add_argument("--n2c", type=Path,
                   default=ROOT / "01_disease_map" / "node_to_celltype.tsv")
    p.add_argument("--output-dir", type=Path,
                   default=ROOT / "01_disease_map" / "celltype_modules")
    p.add_argument("--force-refresh", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stdout,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # --- Données MINERVA ---
    logger.info("Chargement données MINERVA depuis %s", args.cache_dir)
    session = make_session() if args.force_refresh else None
    elements = load_or_fetch_elements(args.cache_dir, session, args.force_refresh)
    reactions = load_or_fetch_reactions(args.cache_dir, session, args.force_refresh)
    elements_by_id = {int(e["id"]): e for e in elements if e.get("id") is not None}
    logger.info("Données : %d éléments, %d réactions", len(elements), len(reactions))

    # --- Graphe global ---
    logger.info("Construction du graphe global…")
    full_graph = build_graph(elements, reactions)

    # --- Assignations cell-type ---
    logger.info("Chargement %s", args.n2c)
    rows = load_node_to_celltype(args.n2c)
    logger.info("Assignations chargées : %d lignes", len(rows))
    by_ct = index_by_celltype(rows)
    extra_set, phen_set = collect_extra_and_phenotype(rows)
    logger.info("EXTRA: %d, PHENOTYPE: %d", len(extra_set), len(phen_set))

    assignments_by_node: dict[int, list] = defaultdict(list)
    for r in rows:
        assignments_by_node[r.node_id].append(r)

    # --- Extraction par cell-type ---
    per_ct_metrics: dict[str, dict[str, Any]] = {}
    gates: dict[str, GateResult] = {}
    for ct in CELL_TYPES:
        core = by_ct[ct]
        if not core:
            logger.warning("Cell-type %s : aucun nœud assigné — skip", ct)
            continue

        logger.info("Module %s : %d nœuds core", ct, len(core))
        module = extract_module(
            celltype=ct,
            core_nodes=core,
            extra_set=extra_set,
            phenotype_set=phen_set,
            reactions=reactions,
            full_graph=full_graph,
        )
        per_ct_metrics[ct] = module.metrics
        gates[ct] = evaluate_gate(module)

        # --- Sorties par module ---
        ct_dir = args.output_dir / ct
        ct_dir.mkdir(parents=True, exist_ok=True)
        write_module_nodes_tsv(
            module, elements_by_id, assignments_by_node,
            ct_dir / f"{ct}_nodes.tsv",
        )
        write_module_reactions_tsv(
            module, elements_by_id,
            ct_dir / f"{ct}_reactions.tsv",
        )
        with (ct_dir / f"{ct}_metrics.json").open("w", encoding="utf-8") as fh:
            json.dump(module.metrics, fh, indent=2, ensure_ascii=False)

        logger.info(
            "  → %s : nodes=%d reactions=%d edges=%d loops=%d gate=%s",
            ct, module.metrics["n_total_nodes"], module.metrics["n_reactions"],
            module.metrics["n_edges"], module.metrics["n_feedback_loops"],
            "PASS" if gates[ct].pass_all else "FAIL",
        )

    # --- Synthèse globale ---
    summary = {
        "n_elements_total": len(elements),
        "n_reactions_total": len(reactions),
        "per_celltype": per_ct_metrics,
        "gates": {ct: {
            "n_nodes": g.n_nodes,
            "n_reactions": g.n_reactions,
            "n_feedback_loops": g.n_feedback_loops,
            "n_phenotype_outputs": g.n_phenotype_outputs,
            "pass_min_nodes": g.pass_min_nodes,
            "pass_min_reactions": g.pass_min_reactions,
            "pass_feedback": g.pass_feedback,
            "pass_phenotype": g.pass_phenotype,
            "pass_all": g.pass_all,
        } for ct, g in gates.items()},
        "n_modules_pass": sum(1 for g in gates.values() if g.pass_all),
    }

    json_path = args.output_dir / "celltype_modules_summary.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    md_path = args.output_dir / "celltype_modules_report.md"
    md_path.write_text(render_report(summary, gates), encoding="utf-8")

    # --- Récap stdout ---
    print()
    print("=" * 70)
    print("  PHASE 1.3 — Cell-type modules")
    print("=" * 70)
    for ct in CELL_TYPES:
        if ct not in gates:
            continue
        g = gates[ct]
        m = per_ct_metrics[ct]
        flag = "PASS" if g.pass_all else "FAIL"
        print(f"  {ct:<8} nodes={m['n_total_nodes']:>4} "
              f"rxn={m['n_reactions']:>4} loops={m['n_feedback_loops']:>4} "
              f"pheno={m['n_phenotype_outputs']:>2}  → {flag}")
    print("-" * 70)
    print(f"  Modules PASS : {summary['n_modules_pass']}/{len(CELL_TYPES)}")
    print(f"  Rapport      : {md_path}")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
