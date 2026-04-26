"""
05_intercellular_edges.py — Phase 1.4 : matérialise les edges intercellulaires
entre les modules cell-type.

Pour chaque entrée de la table curée `lib.intercellular.INTERCELLULAR_EDGES` :

    1. Vérifie que le ligand existe dans la SjD Map (compartiments EXTRA
       21555/21629) — sinon SKIP avec raison.
    2. Vérifie que le récepteur existe dans le module **cible** (rôle=core
       dans `<CT>_nodes.tsv`) — sinon SKIP avec raison.
    3. Émet une ligne par couple (source_celltype, target_celltype, ligand,
       receptor) avec annotations CellPhoneDB/OmniPath/SjS-specific.

Sorties (dans `01_disease_map/`) :
    intercellular_edges.tsv          — 1 ligne par edge instancié
    intercellular_edges_skipped.tsv  — entrées curées non instanciables
    intercellular_edges_summary.json — agrégats + Gate 1.4
    intercellular_edges_report.md    — rapport + Gate 1.4

Usage :
    python scripts/05_intercellular_edges.py
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.intercellular import (  # noqa: E402
    INTERCELLULAR_EDGES,
    IntercellularEdge,
    mandatory_axes_covered,
)

logger = logging.getLogger(__name__)


CELL_TYPES = ("SGEC", "TH1", "TH17", "TFH", "TREG",
              "BCELL", "PLASMA", "M1", "M2", "PDC")


# ---------------------------------------------------------------------------
# Chargement des inventaires
# ---------------------------------------------------------------------------


def load_extracellular_names(path: Path) -> dict[str, list[int]]:
    """Lit `extracellular_nodes.tsv` → {name → [node_id]}."""
    out: dict[str, list[int]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            name = (r.get("node_name") or "").strip()
            if name:
                try:
                    out[name].append(int(r["node_id"]))
                except (ValueError, KeyError):
                    pass
    return dict(out)


def load_module_core_names(modules_dir: Path) -> dict[str, dict[str, list[int]]]:
    """Lit `<CT>/<CT>_nodes.tsv` → {celltype → {name → [node_id]}}, role=core."""
    out: dict[str, dict[str, list[int]]] = {}
    for ct in CELL_TYPES:
        path = modules_dir / ct / f"{ct}_nodes.tsv"
        if not path.exists():
            logger.warning("Module manquant : %s", path)
            out[ct] = {}
            continue
        idx: dict[str, list[int]] = defaultdict(list)
        with path.open(encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader:
                if r.get("role") != "core":
                    continue
                name = (r.get("node_name") or "").strip()
                if name:
                    try:
                        idx[name].append(int(r["node_id"]))
                    except (ValueError, KeyError):
                        pass
        out[ct] = dict(idx)
    return out


# ---------------------------------------------------------------------------
# Matching name (tolérant aux variantes courantes)
# ---------------------------------------------------------------------------


def _name_variants(name: str) -> list[str]:
    """Variantes acceptables pour un nom MINERVA (suffixes/préfixes courants)."""
    base = name.strip()
    variants = {base, base.upper(), base.lower()}
    # Receptor MINERVA names parfois listés comme 'IL21R' vs 'IL-21R' vs 'IL21RA'
    if base.endswith("R"):
        variants.add(base + "A")
        variants.add(base[:-1] + "RA")
    if base.endswith("RA"):
        variants.add(base[:-1])
    return list(variants)


def _resolve_in_index(name: str, index: dict[str, list[int]]) -> list[int]:
    for v in _name_variants(name):
        if v in index:
            return index[v]
    return []


# ---------------------------------------------------------------------------
# Instanciation
# ---------------------------------------------------------------------------


def instantiate_edges(
    extra_index: dict[str, list[int]],
    core_indices: dict[str, dict[str, list[int]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Pour chaque edge curé, émet une ligne par couple (source, target).

    Returns:
        (instantiated, skipped)
    """
    instantiated: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for edge in INTERCELLULAR_EDGES:
        ligand_ids = _resolve_in_index(edge.ligand, extra_index)
        if not ligand_ids:
            skipped.append({
                "ligand": edge.ligand, "receptor": edge.receptor,
                "reason": "ligand_not_in_extracellular",
                "source_celltypes": ",".join(sorted(edge.source_celltypes)),
                "target_celltypes": ",".join(sorted(edge.target_celltypes)),
                "evidence": edge.evidence,
            })
            continue

        for src in sorted(edge.source_celltypes):
            for tgt in sorted(edge.target_celltypes):
                receptor_ids = _resolve_in_index(
                    edge.receptor, core_indices.get(tgt, {})
                )
                if not receptor_ids:
                    skipped.append({
                        "ligand": edge.ligand, "receptor": edge.receptor,
                        "reason": f"receptor_not_in_core({tgt})",
                        "source_celltypes": src,
                        "target_celltypes": tgt,
                        "evidence": edge.evidence,
                    })
                    continue
                instantiated.append({
                    "ligand": edge.ligand,
                    "receptor": edge.receptor,
                    "source_celltype": src,
                    "target_celltype": tgt,
                    "mechanism": edge.mechanism,
                    "ligand_node_ids": ",".join(str(i) for i in ligand_ids),
                    "receptor_node_ids": ",".join(str(i) for i in receptor_ids),
                    "in_cellphonedb": edge.in_cellphonedb,
                    "in_omnipath": edge.in_omnipath,
                    "sjs_specific": edge.sjs_specific,
                    "evidence": edge.evidence,
                })

    return instantiated, skipped


# ---------------------------------------------------------------------------
# Sorties
# ---------------------------------------------------------------------------


def write_edges_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "ligand", "receptor", "source_celltype", "target_celltype",
        "mechanism", "ligand_node_ids", "receptor_node_ids",
        "in_cellphonedb", "in_omnipath", "sjs_specific", "evidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (
            x["source_celltype"], x["target_celltype"], x["ligand"], x["receptor"],
        )):
            r2 = dict(r)
            for k in ("in_cellphonedb", "in_omnipath", "sjs_specific"):
                r2[k] = "1" if r[k] else "0"
            w.writerow(r2)


def write_skipped_tsv(rows: list[dict[str, Any]], path: Path) -> None:
    fields = [
        "ligand", "receptor", "source_celltypes", "target_celltypes",
        "reason", "evidence",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["reason"], x["ligand"])):
            row = {f: r.get(f, "") for f in fields}
            w.writerow(row)


def render_report(
    summary: dict[str, Any],
    coverage: dict[str, bool],
) -> str:
    L: list[str] = []
    L.append("# Intercellular Edges Report — Phase 1.4")
    L.append("")
    L.append("> Sortie de `scripts/05_intercellular_edges.py` — table curée "
             "`scripts/lib/intercellular.py` croisée avec les modules Phase 1.3.")
    L.append("")

    L.append("## 1. Synthèse globale")
    L.append("")
    L.append(f"- Entrées curées : **{summary['n_curated']}**")
    L.append(f"- Edges instanciés : **{summary['n_instantiated']}**")
    L.append(f"- Entrées skippées : **{summary['n_skipped']}**")
    L.append("")

    L.append("## 2. Couverture des axes obligatoires SjD (Gate 1.4)")
    L.append("")
    L.append("| Axe | Statut |")
    L.append("|---|---|")
    for axis, ok in coverage.items():
        L.append(f"| {axis} | {'✅' if ok else '❌'} |")
    L.append("")

    L.append("## 3. Distribution par mécanisme")
    L.append("")
    L.append("| Mécanisme | n_edges |")
    L.append("|---|---|")
    for mech, n in summary["by_mechanism"].items():
        L.append(f"| {mech} | {n} |")
    L.append("")

    L.append("## 4. Distribution par paire (source → target)")
    L.append("")
    L.append("| Source | Target | n_edges |")
    L.append("|---|---|---|")
    for pair_key, n in sorted(
        summary["by_pair"].items(), key=lambda x: -x[1],
    )[:30]:
        s, t = pair_key.split("->", 1)
        L.append(f"| {s} | {t} | {n} |")
    L.append("")
    if len(summary["by_pair"]) > 30:
        L.append(f"_… et {len(summary['by_pair']) - 30} autres paires._")
        L.append("")

    L.append("## 5. Gate 1.4")
    L.append("")
    gate_n = summary["n_instantiated"] >= 30
    gate_axes = all(coverage.values())
    sjs_count = summary["n_sjs_specific"]
    L.append(f"- ≥30 edges intercellulaires : "
             f"**{'PASS' if gate_n else 'FAIL'}** ({summary['n_instantiated']})")
    L.append(f"- 4 axes SjD obligatoires couverts : "
             f"**{'PASS' if gate_axes else 'FAIL'}**")
    L.append(f"- Edges sjs_specific : {sjs_count}")
    L.append("- Revue expert SjD : **EN ATTENTE** (manuel)")
    L.append("")

    L.append("## 6. Sources / provenance")
    L.append("")
    L.append("| Provenance | n_edges |")
    L.append("|---|---|")
    L.append(f"| in_cellphonedb=1 | {summary['n_in_cellphonedb']} |")
    L.append(f"| in_omnipath=1 | {summary['n_in_omnipath']} |")
    L.append(f"| sjs_specific=1 | {summary['n_sjs_specific']} |")
    L.append("")

    L.append("---")
    L.append("")
    L.append("_Généré par `scripts/05_intercellular_edges.py`_")
    return "\n".join(L)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--extra-tsv", type=Path,
                   default=ROOT / "01_disease_map" / "extracellular_nodes.tsv")
    p.add_argument("--modules-dir", type=Path,
                   default=ROOT / "01_disease_map" / "celltype_modules")
    p.add_argument("--output-dir", type=Path,
                   default=ROOT / "01_disease_map")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stdout,
    )

    logger.info("Chargement extracellulaire : %s", args.extra_tsv)
    extra_index = load_extracellular_names(args.extra_tsv)
    logger.info("EXTRA names indexés : %d", len(extra_index))

    logger.info("Chargement modules : %s", args.modules_dir)
    core_indices = load_module_core_names(args.modules_dir)
    for ct, idx in core_indices.items():
        logger.info("  %s : %d noms core", ct, len(idx))

    logger.info("Instanciation des edges curés (%d entrées)…",
                len(INTERCELLULAR_EDGES))
    instantiated, skipped = instantiate_edges(extra_index, core_indices)

    # --- Stats ---
    by_mechanism = Counter(e["mechanism"] for e in instantiated)
    by_pair = Counter(
        (e["source_celltype"], e["target_celltype"]) for e in instantiated
    )
    n_cpdb = sum(1 for e in instantiated if e["in_cellphonedb"])
    n_op = sum(1 for e in instantiated if e["in_omnipath"])
    n_sjs = sum(1 for e in instantiated if e["sjs_specific"])
    coverage = mandatory_axes_covered(instantiated)

    summary = {
        "n_curated": len(INTERCELLULAR_EDGES),
        "n_instantiated": len(instantiated),
        "n_skipped": len(skipped),
        "n_in_cellphonedb": n_cpdb,
        "n_in_omnipath": n_op,
        "n_sjs_specific": n_sjs,
        "by_mechanism": dict(by_mechanism.most_common()),
        "by_pair": {f"{s}->{t}": n for (s, t), n in by_pair.most_common()},
        "mandatory_coverage": coverage,
        "gate_min_edges_pass": len(instantiated) >= 30,
        "gate_axes_pass": all(coverage.values()),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    edges_path = args.output_dir / "intercellular_edges.tsv"
    skip_path = args.output_dir / "intercellular_edges_skipped.tsv"
    json_path = args.output_dir / "intercellular_edges_summary.json"
    md_path = args.output_dir / "intercellular_edges_report.md"

    write_edges_tsv(instantiated, edges_path)
    write_skipped_tsv(skipped, skip_path)
    with json_path.open("w", encoding="utf-8") as fh:
        # by_pair tuple keys → string
        json.dump({**summary, "by_pair": summary["by_pair"]},
                  fh, indent=2, ensure_ascii=False)
    md_path.write_text(render_report(summary, coverage), encoding="utf-8")

    # --- Récap stdout ---
    print()
    print("=" * 70)
    print("  PHASE 1.4 — Intercellular edges")
    print("=" * 70)
    print(f"  Entrées curées      : {summary['n_curated']}")
    print(f"  Edges instanciés    : {summary['n_instantiated']}")
    print(f"  Skippés             : {summary['n_skipped']}")
    print(f"  CellPhoneDB / OP    : {n_cpdb} / {n_op}")
    print(f"  SjS-spécifiques     : {n_sjs}")
    print("-" * 70)
    print("  Axes obligatoires :")
    for ax, ok in coverage.items():
        print(f"    {'✅' if ok else '❌'}  {ax}")
    print("-" * 70)
    gate_pass = summary["gate_min_edges_pass"] and summary["gate_axes_pass"]
    print(f"  Gate 1.4 : {'PASS' if gate_pass else 'FAIL'}")
    print(f"  Rapport  : {md_path}")
    print("=" * 70)

    return 0


if __name__ == "__main__":
    sys.exit(main())
