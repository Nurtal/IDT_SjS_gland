"""
02_audit_map.py — Audit topologique et sémantique de la SjD Map.

Phase 1.1 de la ROADMAP : caractérise la carte avant dissociation cell-type.

Sorties :
    01_disease_map/audit_summary.json   — synthèse machine-lisible
    01_disease_map/audit_report.md      — rapport humain
    01_disease_map/audit_celltype_hits.tsv — marqueurs cell-type détectés
    01_disease_map/audit_intercellular_pivots.tsv — cytokines clés détectées

Usage :
    python scripts/02_audit_map.py [--cache-dir 01_disease_map/cache]
                                    [--output-dir 01_disease_map]

Le script lit les caches JSON produits par 01_download_sjd_map.py.
Si les caches sont absents, il appelle l'API MINERVA (cf. 01_download_sjd_map.py).
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

from lib.map_audit import (
    CELLTYPE_MARKERS,
    SJD_CYTOKINES_OF_INTEREST,
    full_audit,
)
from lib.minerva_api import (
    load_or_fetch_elements,
    load_or_fetch_reactions,
    make_session,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Sérialisation JSON robuste (Counter, set...)
# ---------------------------------------------------------------------------


class _AuditJSONEncoder(json.JSONEncoder):
    def default(self, obj: Any) -> Any:
        if isinstance(obj, set):
            return sorted(obj)
        try:
            return super().default(obj)
        except TypeError:
            return str(obj)


# ---------------------------------------------------------------------------
# Génération du rapport markdown
# ---------------------------------------------------------------------------


def _markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    out = ["| " + " | ".join(headers) + " |"]
    out.append("|" + "|".join(["---"] * len(headers)) + "|")
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(out)


def render_markdown_report(audit: dict[str, Any]) -> str:
    topo = audit["topology"]
    el_ann = audit["element_annotations"]
    rxn_ann = audit["reaction_annotations"]
    comps = audit["compartments"]
    hits = audit["celltype_marker_hits"]
    pivots = audit["intercellular_pivots"]

    lines: list[str] = []
    lines.append("# Audit topologique & sémantique — SjD Map")
    lines.append("")
    lines.append("> **Phase 1.1 ROADMAP** — caractérisation de la SjD Map "
                 "avant dissociation cell-type (Phase 1.2).")
    lines.append("")

    # === Topologie ===
    lines.append("## 1. Topologie")
    lines.append("")
    lines.append(_markdown_table(
        ["Métrique", "Valeur"],
        [
            ["Nombre de nœuds", topo["n_nodes"]],
            ["Nombre d'edges", topo["n_edges"]],
            ["Composantes faibles", topo["n_weakly_connected_components"]],
            ["Plus grosse composante faible", topo["largest_weak_component_size"]],
            ["Composantes fortes", topo["n_strongly_connected_components"]],
            ["Plus grosse composante forte", topo["largest_strong_component_size"]],
            ["Nœuds isolés (degree 0)", topo["n_isolated_nodes"]],
            ["Nœuds sources (in=0)", topo["n_source_nodes"]],
            ["Nœuds puits (out=0)", topo["n_sink_nodes"]],
            [f"Hubs (degree ≥ {topo['hub_threshold']})", topo["n_hubs"]],
        ],
    ))
    lines.append("")
    lines.append("### Top 20 hubs")
    lines.append("")
    if topo["top_20_hubs"]:
        rows = [
            [h["id"], h["name"], h["type"], h["in_degree"], h["out_degree"], h["total_degree"]]
            for h in topo["top_20_hubs"]
        ]
        lines.append(_markdown_table(
            ["id", "name", "type", "in", "out", "total"],
            rows,
        ))
    else:
        lines.append("_aucun hub_")
    lines.append("")

    # === Annotations ===
    lines.append("## 2. Couverture des annotations")
    lines.append("")
    lines.append(f"**Éléments** : {el_ann['total_elements']}")
    lines.append("")
    lines.append(_markdown_table(
        ["Annotation", "Nb éléments"],
        [[k, v] for k, v in el_ann["by_annotation"].items()],
    ))
    lines.append("")
    lines.append(f"- Avec ≥1 annotation externe : **{el_ann['with_any_annotation']}** "
                 f"({el_ann['coverage_pct_any']}%)")
    lines.append(f"- Sans aucune annotation : **{el_ann['without_annotation']}**")
    lines.append(f"- Avec PubMed ID : **{el_ann['with_pmid']}**")
    lines.append("")
    lines.append(f"**Réactions** : {rxn_ann['total_reactions']}")
    lines.append(f"- Avec PubMed ID : **{rxn_ann['with_pmid']}** "
                 f"({rxn_ann['pmid_coverage_pct']}%)")
    lines.append("")

    # === Gate 1.1 ===
    gate_any = el_ann["coverage_pct_any"] >= 80.0
    gate_pmid = rxn_ann["pmid_coverage_pct"] >= 50.0
    lines.append("### Gate 1.1")
    lines.append("")
    lines.append(f"- ≥80% éléments avec annotation : **{'PASS' if gate_any else 'FAIL'}** "
                 f"({el_ann['coverage_pct_any']}%)")
    lines.append(f"- ≥50% réactions avec PMID : **{'PASS' if gate_pmid else 'FAIL'}** "
                 f"({rxn_ann['pmid_coverage_pct']}%)")
    lines.append("")

    # === Compartiments ===
    lines.append("## 3. Compartiments")
    lines.append("")
    rows = [
        [str(cid), info["n_species"], ", ".join(list(info["by_type"].keys())[:4])]
        for cid, info in sorted(comps.items(), key=lambda x: -x[1]["n_species"])
    ]
    lines.append(_markdown_table(
        ["compartmentId", "n_species", "top types"],
        rows,
    ))
    lines.append("")

    # === Marqueurs cell-type ===
    lines.append("## 4. Marqueurs cell-type détectés")
    lines.append("")
    lines.append("Pour chaque cell-type cible, nombre de marqueurs canoniques "
                 "(cf. CONVENTIONS.md) trouvés par correspondance de nom :")
    lines.append("")
    rows = []
    for ct, found in hits.items():
        total_markers = len(CELLTYPE_MARKERS[ct])
        rows.append([
            ct,
            len(found),
            total_markers,
            f"{round(100*len(found)/total_markers, 1)}%",
            ", ".join(h["name"] for h in found[:6]),
        ])
    lines.append(_markdown_table(
        ["Cell-type", "trouvés", "attendus", "%", "exemples"],
        rows,
    ))
    lines.append("")

    # === Pivots intercellulaires ===
    lines.append("## 5. Pivots intercellulaires (cytokines/chemokines clés)")
    lines.append("")
    lines.append(f"**{len(pivots)} / {len(SJD_CYTOKINES_OF_INTEREST)}** "
                 "cytokines d'intérêt SjD trouvées dans la carte.")
    lines.append("")
    if pivots:
        rows = [[p["id"], p["name"], p["type"], p["matched_cytokine"]] for p in pivots]
        lines.append(_markdown_table(
            ["id", "name", "type", "matched"],
            rows,
        ))
    lines.append("")

    found_cytokines = {p["matched_cytokine"] for p in pivots}
    missing = SJD_CYTOKINES_OF_INTEREST - found_cytokines
    if missing:
        lines.append(f"**Cytokines clés manquantes** ({len(missing)}) : "
                     + ", ".join(sorted(missing)))
        lines.append("")

    # === Pathway hints ===
    pathway_hints = audit["pathway_hints_from_notes"]
    lines.append("## 6. Indices de pathway extraits des notes")
    lines.append("")
    lines.append(f"**{len(pathway_hints)}** éléments avec ≥1 identifiant "
                 "Reactome / KEGG / GO dans leurs notes.")
    lines.append("")

    # === Recommandations Phase 1.2 ===
    lines.append("## 7. Recommandations pour Phase 1.2 (dissociation)")
    lines.append("")

    # Compartiment
    n_compartments = len([c for c in comps if c is not None])
    if n_compartments >= 5:
        lines.append(f"- ✅ **{n_compartments} compartiments** présents → "
                     "approche compartiment-driven viable comme première passe.")
    else:
        lines.append(f"- ⚠️ Seulement **{n_compartments} compartiments** → "
                     "approche compartiment-driven insuffisante seule.")

    # Markers
    well_covered = [ct for ct, h in hits.items() if len(h) >= 3]
    poorly_covered = [ct for ct, h in hits.items() if len(h) < 2]
    if well_covered:
        lines.append(f"- ✅ Marqueurs bien représentés : {', '.join(well_covered)}")
    if poorly_covered:
        lines.append(f"- ⚠️ Marqueurs peu représentés : {', '.join(poorly_covered)} "
                     "→ envisager extension littérature ou retrait du cell-type.")

    # PDC
    pdc_hits = len(hits.get("PDC", []))
    if pdc_hits < 2:
        lines.append("- ⚠️ pDC sous-représentées → décision Phase 1.3 : "
                     "exclure ou compléter manuellement (TLR7/9, IRF7, IL3RA).")

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"_Généré par `scripts/02_audit_map.py`_")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Sérialisation TSV
# ---------------------------------------------------------------------------


def write_celltype_hits_tsv(
    hits: dict[str, list[dict[str, Any]]], path: Path
) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["celltype", "matched_marker", "node_id", "node_name", "node_type"])
        for ct, items in hits.items():
            for h in items:
                w.writerow([
                    ct, h["matched_marker"], h["id"], h["name"], h["type"],
                ])


def write_pivots_tsv(pivots: list[dict[str, Any]], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["matched_cytokine", "node_id", "node_name", "node_type"])
        for p in pivots:
            w.writerow([
                p["matched_cytokine"], p["id"], p["name"], p["type"],
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=ROOT / "01_disease_map" / "cache",
        help="Répertoire de cache JSON MINERVA",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "01_disease_map",
        help="Répertoire de sortie pour le rapport",
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Ignore le cache et ré-interroge MINERVA",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Active les logs DEBUG",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
        stream=sys.stdout,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Chargement des données MINERVA (cache : %s)", args.cache_dir)
    session = make_session() if args.force_refresh else None
    elements = load_or_fetch_elements(args.cache_dir, session, args.force_refresh)
    reactions = load_or_fetch_reactions(args.cache_dir, session, args.force_refresh)
    logger.info("Données : %d éléments, %d réactions", len(elements), len(reactions))

    logger.info("Lancement de l'audit...")
    audit = full_audit(elements, reactions)

    # Sérialisation
    json_path = args.output_dir / "audit_summary.json"
    md_path = args.output_dir / "audit_report.md"
    hits_path = args.output_dir / "audit_celltype_hits.tsv"
    pivots_path = args.output_dir / "audit_intercellular_pivots.tsv"

    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, ensure_ascii=False, cls=_AuditJSONEncoder)
    logger.info("JSON écrit : %s", json_path)

    md_path.write_text(render_markdown_report(audit), encoding="utf-8")
    logger.info("Rapport markdown écrit : %s", md_path)

    write_celltype_hits_tsv(audit["celltype_marker_hits"], hits_path)
    logger.info("Marqueurs cell-type écrits : %s", hits_path)

    write_pivots_tsv(audit["intercellular_pivots"], pivots_path)
    logger.info("Pivots intercellulaires écrits : %s", pivots_path)

    # Petit récap stdout
    topo = audit["topology"]
    el_ann = audit["element_annotations"]
    print()
    print("=" * 60)
    print("  AUDIT — SjD Map")
    print("=" * 60)
    print(f"  Nœuds : {topo['n_nodes']}    Edges : {topo['n_edges']}")
    print(f"  Composantes faibles : {topo['n_weakly_connected_components']}")
    print(f"  Annotation coverage : {el_ann['coverage_pct_any']}%")
    print(f"  Pivots cytokines    : {len(audit['intercellular_pivots'])}/"
          f"{len(SJD_CYTOKINES_OF_INTEREST)}")
    print(f"  Rapport : {md_path}")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
