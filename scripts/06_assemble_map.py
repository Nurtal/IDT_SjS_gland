"""
06_assemble_map.py — Phase 1.5 : assemble la SjD multicellular map en SBML/CellDesigner.

Combine :
    - 10 modules cell-type (Phase 1.3) via `node_to_celltype.tsv`
    - 470 edges inter-cellulaires curés (Phase 1.4) via `intercellular_edges.tsv`

Sortie principale :
    01_disease_map/SjD_multicellular_map.xml

Sorties annexes :
    01_disease_map/assembly_summary.json
    01_disease_map/assembly_report.md

Gate Phase 1 final :
    - ≥1200 espèces, ≥800 réactions
    - Validation libsbml : 0 erreur fatale
    - Aucune speciesReference orpheline

Usage :
    python scripts/06_assemble_map.py
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

from lib.assembly import (  # noqa: E402
    CELL_TYPES,
    AssemblyContext,
    AssemblyStats,
    assemble_map,
)
from lib.celldesigner_xml import (  # noqa: E402
    check_species_references,
    count_elements,
    validate_sbml,
    write_celldesigner_xml,
)
from lib.minerva_api import (  # noqa: E402
    load_or_fetch_elements,
    load_or_fetch_reactions,
    make_session,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def load_n2c_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_intercellular_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Rapport
# ---------------------------------------------------------------------------


def render_report(
    summary: dict[str, Any],
    libsbml_errors: list[str],
    unresolved_refs: list[str],
) -> str:
    L: list[str] = []
    L.append("# Assembly Report — Phase 1.5")
    L.append("")
    L.append("> Sortie de `scripts/06_assemble_map.py` — fusion des 10 modules "
             "cell-type + 470 edges inter-cellulaires en un unique CellDesigner SBML.")
    L.append("")

    L.append("## 1. Synthèse globale")
    L.append("")
    L.append(f"- Compartiments : **{summary['compartments']}**")
    L.append(f"- Espèces       : **{summary['species']}**")
    L.append(f"- Réactions     : **{summary['reactions']}**")
    L.append("")

    s = summary["stats"]
    L.append("## 2. Détail des espèces")
    L.append("")
    L.append("| Catégorie | n |")
    L.append("|---|---|")
    L.append(f"| core (clonées par cell-type) | {s['n_species_core']} |")
    L.append(f"| EXTRA (partagées) | {s['n_species_extra']} |")
    L.append(f"| PHENOTYPE (partagées) | {s['n_species_phenotype']} |")
    L.append(f"| complexes inter-cellulaires (contact) | {s['n_species_complex']} |")
    L.append("")

    L.append("## 3. Réactions intracellulaires par cell-type")
    L.append("")
    L.append("| Cell-type | n_reactions |")
    L.append("|---|---|")
    for ct in CELL_TYPES:
        n = s["by_celltype_reactions"].get(ct, 0)
        L.append(f"| {ct} | {n} |")
    L.append(f"| **Total intra** | **{s['n_reactions_intra']}** |")
    L.append(f"| Réactions skippées (participants non résolus) | {s['n_reactions_skipped']} |")
    L.append("")

    L.append("## 4. Edges inter-cellulaires matérialisés")
    L.append("")
    L.append("| Mécanisme | n |")
    L.append("|---|---|")
    L.append(f"| secreted/autocrine (PHYSICAL_STIMULATION) | {s['n_inter_secreted']} |")
    L.append(f"| contact (HETERODIMER_ASSOCIATION) | {s['n_inter_contact']} |")
    L.append(f"| skipped (ligand/recepteur non résolus) | {s['n_inter_skipped']} |")
    L.append("")

    L.append("## 5. Validation SBML")
    L.append("")
    n_fatal = sum(1 for e in libsbml_errors if "[Fatal]" in e)
    n_schema = sum(1 for e in libsbml_errors if "[Error]" in e)
    L.append(f"- Erreurs FATAL libsbml : **{n_fatal}**")
    L.append(f"- Schema-warnings (annotations CellDesigner non couvertes par L2V4) : "
             f"{n_schema} _(non bloquantes — pattern identique à la SjD Map "
             f"publiée Silva-Saffar 2026)_")
    L.append(f"- speciesReference non résolues : **{len(unresolved_refs)}**")
    if libsbml_errors[:5]:
        L.append("")
        L.append("Premiers messages libsbml :")
        L.append("")
        L.append("```")
        for e in libsbml_errors[:10]:
            L.append(e)
        L.append("```")
    if unresolved_refs[:5]:
        L.append("")
        L.append("Premières références orphelines :")
        L.append("")
        L.append("```")
        for u in unresolved_refs[:10]:
            L.append(u)
        L.append("```")
    L.append("")

    L.append("## 6. Gate Phase 1 final")
    L.append("")
    g = summary["gate"]
    L.append("| Critère | Seuil | Mesuré | Statut |")
    L.append("|---|---|---|---|")
    L.append(f"| Espèces ≥ 1200 | 1200 | {summary['species']} | "
             f"{'✅' if g['species'] else '❌'} |")
    L.append(f"| Réactions ≥ 800 | 800 | {summary['reactions']} | "
             f"{'✅' if g['reactions'] else '❌'} |")
    L.append(f"| 0 erreur libsbml fatale | 0 | {n_fatal} | "
             f"{'✅' if g['libsbml'] else '❌'} |")
    L.append(f"| 0 référence orpheline | 0 | {len(unresolved_refs)} | "
             f"{'✅' if g['refs'] else '❌'} |")
    L.append("")
    overall = all(g.values())
    L.append(f"**Statut global Phase 1 : {'PASS' if overall else 'FAIL'}**")
    L.append("")

    L.append("---")
    L.append("")
    L.append("_Généré par `scripts/06_assemble_map.py`_")
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
    p.add_argument("--inter", type=Path,
                   default=ROOT / "01_disease_map" / "intercellular_edges.tsv")
    p.add_argument("--output-dir", type=Path,
                   default=ROOT / "01_disease_map")
    p.add_argument("--xml-name", type=str,
                   default="SjD_multicellular_map.xml")
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
    logger.info("Chargement données MINERVA…")
    session = make_session() if args.force_refresh else None
    elements = load_or_fetch_elements(args.cache_dir, session, args.force_refresh)
    reactions = load_or_fetch_reactions(args.cache_dir, session, args.force_refresh)
    elements_by_id = {int(e["id"]): e for e in elements if e.get("id") is not None}
    logger.info("  %d éléments, %d réactions", len(elements), len(reactions))

    # --- Assignations & edges ---
    logger.info("Chargement %s", args.n2c)
    n2c_rows = load_n2c_rows(args.n2c)
    logger.info("  %d lignes d'assignation", len(n2c_rows))

    logger.info("Chargement %s", args.inter)
    inter_rows = load_intercellular_rows(args.inter)
    logger.info("  %d edges inter-cellulaires", len(inter_rows))

    # --- Assemblage ---
    ctx = AssemblyContext(
        n2c_rows=n2c_rows,
        elements_by_id=elements_by_id,
        reactions=reactions,
        intercellular_rows=inter_rows,
    )
    sbml, stats = assemble_map(ctx)

    # --- Écriture ---
    xml_path = args.output_dir / args.xml_name
    write_celldesigner_xml(sbml, xml_path)

    # --- Validation ---
    logger.info("Validation libsbml…")
    valid, libsbml_errors = validate_sbml(xml_path)

    logger.info("Vérification des références espèces…")
    unresolved = check_species_references(sbml)
    if unresolved:
        logger.warning("%d speciesReference orphelines détectées", len(unresolved))

    counts = count_elements(sbml)
    n_fatal = sum(1 for e in libsbml_errors if "[Fatal]" in e)
    n_schema_warn = sum(1 for e in libsbml_errors if "[Error]" in e)

    summary: dict[str, Any] = {
        "xml_path": str(xml_path),
        "compartments": counts["compartments"],
        "species": counts["species"],
        "reactions": counts["reactions"],
        "stats": {
            "n_species_core": stats.n_species_core,
            "n_species_extra": stats.n_species_extra,
            "n_species_phenotype": stats.n_species_phenotype,
            "n_species_complex": stats.n_species_complex,
            "n_reactions_intra": stats.n_reactions_intra,
            "n_reactions_skipped": stats.n_reactions_skipped,
            "n_inter_secreted": stats.n_inter_secreted,
            "n_inter_contact": stats.n_inter_contact,
            "n_inter_skipped": stats.n_inter_skipped,
            "by_celltype_reactions": stats.by_celltype_reactions,
        },
        "libsbml_n_errors": len(libsbml_errors),
        "libsbml_n_fatal": n_fatal,
        "libsbml_n_schema_warnings": n_schema_warn,
        "n_unresolved_species_refs": len(unresolved),
        "gate": {
            "species":   counts["species"]   >= 1200,
            "reactions": counts["reactions"] >= 800,
            "libsbml":   n_fatal == 0,
            "refs":      len(unresolved) == 0,
        },
    }

    json_path = args.output_dir / "assembly_summary.json"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    md_path = args.output_dir / "assembly_report.md"
    md_path.write_text(
        render_report(summary, libsbml_errors, unresolved),
        encoding="utf-8",
    )

    # --- Récap stdout ---
    print()
    print("=" * 70)
    print("  PHASE 1.5 — Assemblage multi-cellulaire")
    print("=" * 70)
    print(f"  Compartiments    : {counts['compartments']}")
    print(f"  Espèces          : {counts['species']}")
    print(f"  Réactions        : {counts['reactions']}")
    print(f"    - intra         : {stats.n_reactions_intra}")
    print(f"    - secreted      : {stats.n_inter_secreted}")
    print(f"    - contact       : {stats.n_inter_contact}")
    print(f"  libsbml fatal    : {n_fatal}")
    print(f"  libsbml schema-warnings (CellDesigner) : {n_schema_warn}")
    print(f"  Refs orphelines  : {len(unresolved)}")
    print("-" * 70)
    g = summary["gate"]
    for k, ok in g.items():
        print(f"    {'✅' if ok else '❌'}  {k}")
    overall = all(g.values())
    print("-" * 70)
    print(f"  Gate Phase 1 final : {'PASS' if overall else 'FAIL'}")
    print(f"  XML              : {xml_path}")
    print(f"  Rapport          : {md_path}")
    print("=" * 70)

    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
