"""
03_dissociate.py — Phase 1.2 : assigne chaque nœud de la SjD Map à un ou
plusieurs cell-types selon les règles R1–R7 (cf. dissociation_rules.md).

Sorties produites dans `01_disease_map/` :
    node_to_celltype.tsv         — format long, multi-assignment
    extracellular_nodes.tsv      — sous-ensemble R1
    unassigned_nodes.tsv         — sous-ensemble R7
    dissociation_summary.json    — statistiques agrégées
    dissociation_report.md       — rapport humain + Gate 1.2

Usage :
    python scripts/03_dissociate.py [--cache-dir 01_disease_map/cache]
                                     [--output-dir 01_disease_map]
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.dissociator import (  # noqa: E402
    CELL_TYPES,
    EXTRA,
    PHENOTYPE,
    UNASSIGNED,
    Assignment,
    collect_extracellular,
    dissociate,
    summarize,
)
from lib.map_audit import build_graph  # noqa: E402
from lib.minerva_api import (  # noqa: E402
    load_or_fetch_elements,
    load_or_fetch_reactions,
    make_session,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sorties TSV
# ---------------------------------------------------------------------------


def write_node_to_celltype(
    assignments: dict[int, Assignment],
    path: Path,
) -> None:
    fields = [
        "node_id", "node_name", "node_type", "compartment_id",
        "celltype", "confidence", "rule", "evidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for nid in sorted(assignments.keys()):
            for row in assignments[nid].to_rows():
                w.writerow(row)


def write_extracellular(
    rows: list[dict[str, Any]],
    path: Path,
) -> None:
    fields = [
        "node_id", "node_name", "node_type", "compartment_id",
        "canonical_form", "g_r_p_variants",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in sorted(rows, key=lambda r: r["canonical_form"]):
            w.writerow(row)


def write_unassigned(
    assignments: dict[int, Assignment],
    path: Path,
) -> None:
    fields = ["node_id", "node_name", "node_type", "compartment_id"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for nid, a in sorted(assignments.items()):
            if a.rule != "R7":
                continue
            w.writerow({
                "node_id": nid,
                "node_name": a.name,
                "node_type": a.type,
                "compartment_id": a.compartment_id,
            })


# ---------------------------------------------------------------------------
# Rapport markdown
# ---------------------------------------------------------------------------


def render_report(summary: dict[str, Any]) -> str:
    L: list[str] = []
    L.append("# Dissociation Report — SjD Map → cell-types")
    L.append("")
    L.append("> **Phase 1.2 ROADMAP** — sortie de `scripts/03_dissociate.py`.")
    L.append("> Règles définies dans `dissociation_rules.md`.")
    L.append("")

    # Synthèse
    L.append("## 1. Couverture globale")
    L.append("")
    L.append(f"- Nœuds totaux : **{summary['n_total_nodes']}**")
    L.append(f"- Extracellulaires (R1) : **{summary['n_extracellular']}**")
    L.append(f"- Phénotypes (R2) : **{summary['n_phenotype']}**")
    L.append(f"- Inassignables (R7) : **{summary['n_unassigned']}**")
    L.append(f"- Assignables (hors EXTRA/PHENOTYPE) : **{summary['n_assignable']}**")
    L.append(f"- Assignés à ≥1 cell-type : **{summary['n_assigned_to_celltype']}** "
             f"({summary['coverage_pct']}%)")
    L.append("")

    # Par règle
    L.append("## 2. Distribution par règle")
    L.append("")
    L.append("| Règle | n_nœuds |")
    L.append("|---|---|")
    for rule, n in sorted(summary["by_rule"].items()):
        L.append(f"| {rule} | {n} |")
    L.append("")

    # Par confiance
    L.append("## 3. Distribution par confiance")
    L.append("")
    L.append("| Confidence | n_nœuds |")
    L.append("|---|---|")
    for conf, n in summary["by_confidence"].items():
        L.append(f"| {conf or '(R7)'} | {n} |")
    L.append("")

    # Par cell-type
    L.append("## 4. Effectifs par cell-type (clones potentiels)")
    L.append("")
    L.append("| Cell-type | n_nœuds | seuil 80 |")
    L.append("|---|---|---|")
    by_ct = summary["by_celltype"]
    threshold = summary["celltype_threshold"]
    for ct in (*CELL_TYPES, EXTRA, PHENOTYPE):
        n = by_ct.get(ct, 0)
        if ct in {EXTRA, PHENOTYPE}:
            mark = "—"
        else:
            mark = "✅" if n >= threshold else "❌"
        L.append(f"| {ct} | {n} | {mark} |")
    L.append("")

    # Gate 1.2
    L.append("## 5. Gate 1.2")
    L.append("")
    L.append(f"- Couverture ≥90% (assignables) : "
             f"**{'PASS' if summary['gate_coverage_pass'] else 'FAIL'}** "
             f"({summary['coverage_pct']}%)")
    if summary["underprovisioned_celltypes"]:
        L.append(f"- Cell-types <80 nœuds : "
                 f"**FAIL** → {', '.join(summary['underprovisioned_celltypes'])}")
    else:
        L.append("- Tous les cell-types ont ≥80 nœuds : **PASS**")
    L.append("- Revue expert : **EN ATTENTE** (manuel)")
    L.append("")

    # Recommandations
    L.append("## 6. Recommandations Phase 1.3")
    L.append("")
    if summary["underprovisioned_celltypes"]:
        L.append("- Cell-types sous-provisionnés à enrichir manuellement (curation littérature) :")
        for ct in summary["underprovisioned_celltypes"]:
            L.append(f"    - **{ct}** (n={by_ct.get(ct, 0)}) → ajouter marqueurs "
                     f"et nœuds depuis Reactome/UniProt curation manuelle.")
    if summary["n_unassigned"] > 0:
        L.append(f"- Réviser les **{summary['n_unassigned']} nœuds inassignables** "
                 f"(`unassigned_nodes.tsv`) — soit ajouter un marqueur, soit retirer.")
    if summary["coverage_pct"] < 90.0:
        L.append("- Couverture <90% → ré-itérer les règles R5/R6 ou enrichir R3/R4.")
    L.append("")

    L.append("---")
    L.append("")
    L.append("_Généré par `scripts/03_dissociate.py`_")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache-dir", type=Path,
                   default=ROOT / "01_disease_map" / "cache")
    p.add_argument("--output-dir", type=Path,
                   default=ROOT / "01_disease_map")
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

    logger.info("Chargement données MINERVA depuis %s", args.cache_dir)
    session = make_session() if args.force_refresh else None
    elements = load_or_fetch_elements(args.cache_dir, session, args.force_refresh)
    reactions = load_or_fetch_reactions(args.cache_dir, session, args.force_refresh)
    logger.info("Données : %d éléments, %d réactions", len(elements), len(reactions))

    logger.info("Construction du graphe…")
    graph = build_graph(elements, reactions)

    logger.info("Application des règles de dissociation R1–R7…")
    assignments = dissociate(elements, graph)
    summary = summarize(assignments)

    # --- TSV ---
    n2c_path = args.output_dir / "node_to_celltype.tsv"
    extra_path = args.output_dir / "extracellular_nodes.tsv"
    unass_path = args.output_dir / "unassigned_nodes.tsv"
    write_node_to_celltype(assignments, n2c_path)
    write_extracellular(collect_extracellular(assignments, elements), extra_path)
    write_unassigned(assignments, unass_path)

    # --- JSON ---
    json_path = args.output_dir / "dissociation_summary.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    # --- Rapport ---
    md_path = args.output_dir / "dissociation_report.md"
    md_path.write_text(render_report(summary), encoding="utf-8")

    # Récap stdout
    print()
    print("=" * 60)
    print("  DISSOCIATION — SjD Map")
    print("=" * 60)
    print(f"  Nœuds totaux    : {summary['n_total_nodes']}")
    print(f"  EXTRA / PHENO   : {summary['n_extracellular']} / {summary['n_phenotype']}")
    print(f"  Assignables     : {summary['n_assignable']}")
    print(f"  Couverture      : {summary['coverage_pct']}%  "
          f"(gate {'PASS' if summary['gate_coverage_pass'] else 'FAIL'})")
    print(f"  Inassignables   : {summary['n_unassigned']}")
    print(f"  Sous-provision  : {summary['underprovisioned_celltypes']}")
    print(f"  Rapport         : {md_path}")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
