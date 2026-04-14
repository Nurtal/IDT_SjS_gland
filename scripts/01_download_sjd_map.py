"""
01_download_sjd_map.py — Reconstruit SjD_Map_original.xml depuis l'API MINERVA.

Étapes :
    1. Interroge l'API MINERVA (avec cache local)
    2. Construit tous les compartiments
    3. Construit toutes les espèces (species)
    4. Construit toutes les réactions
    5. Écrit 01_disease_map/SjD_Map_original.xml
    6. Valide avec libsbml
    7. Affiche un résumé

Usage :
    python scripts/01_download_sjd_map.py [--force-refresh] [--output PATH]

Options :
    --force-refresh   Ignore le cache JSON et ré-interroge l'API MINERVA
    --output PATH     Chemin de sortie (défaut : 01_disease_map/SjD_Map_original.xml)
    --verbose         Active les logs DEBUG
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Résolution des imports relatifs depuis scripts/
ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from lib.minerva_api import (
    make_session,
    load_or_fetch_elements,
    load_or_fetch_reactions,
    summarize_elements,
    summarize_reactions,
    get_phenotype_nodes,
    COMPARTMENT_NAMES,
)
from lib.celldesigner_xml import (
    make_sbml_root,
    add_compartment,
    make_species,
    add_species,
    make_reaction,
    add_reaction,
    write_celldesigner_xml,
    validate_sbml,
    count_elements,
    check_species_references,
    MINERVA_TYPE_TO_CD,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Nœuds phénotypiques attendus (d'après la publication SjD Map)
# ---------------------------------------------------------------------------

EXPECTED_PHENOTYPES = {
    "Inflammation",
    "Regulated_Necrosis",
    "Phagocytosis",
    "Chemotaxis/Infiltration",
    "Fibrosis",
    "B_Cell_Activation/Survival",
    "Lymphoid_organ_development",
    "MHC_Class_2_Activation",
    "Angiogenesis",
    "Matrix_degradation",
    "MHC_Class_1_Activation",
    "Apoptosis",
    "T_Cell_Activation/Differentiation",
    "Cell_Proliferation/Survival",
}

# ---------------------------------------------------------------------------
# Construction de la carte
# ---------------------------------------------------------------------------

def _collect_compartment_ids(elements: list[dict]) -> set[int | None]:
    """Collecte tous les compartmentId uniques présents dans les éléments."""
    return {e.get("compartmentId") for e in elements}


def build_compartments(
    sbml_root,
    elements: list[dict],
) -> dict[int | None, str]:
    """
    Ajoute tous les compartiments au document SBML.

    Returns:
        {minerva_compartment_id → sbml_compartment_id}
    """
    comp_ids = _collect_compartment_ids(elements)
    mapping: dict[int | None, str] = {}

    # Ordre déterministe : d'abord None/Default, puis par ID numérique croissant
    ordered = sorted(comp_ids, key=lambda x: (x is not None, x or 0))

    for cid in ordered:
        comp_el = add_compartment(sbml_root, cid, module_prefix="")
        mapping[cid] = comp_el.get("id")
        logger.debug("Compartiment : %s → %s", cid, comp_el.get("id"))

    logger.info("%d compartiments créés", len(mapping))
    return mapping


def build_species(
    sbml_root,
    elements: list[dict],
    comp_map: dict[int | None, str],
) -> dict[int, str]:
    """
    Ajoute toutes les espèces au document SBML.

    Returns:
        {minerva_element_id (int) → sbml_species_id (str)}
    """
    id_map: dict[int, str] = {}
    skipped = 0

    for elem in elements:
        eid = elem.get("id")
        if eid is None:
            skipped += 1
            continue

        comp_id = elem.get("compartmentId")
        sbml_comp = comp_map.get(comp_id)
        if sbml_comp is None:
            # Compartiment inconnu → assignation au compartiment Default
            sbml_comp = comp_map.get(None, "Default")
            logger.debug(
                "Élément %s (%s) : compartimentId %s inconnu, assigné à Default",
                eid, elem.get("name"), comp_id,
            )

        sp_el = make_species(elem, sbml_comp, module_prefix="")
        add_species(sbml_root, sp_el)
        id_map[int(eid)] = sp_el.get("id")

    logger.info(
        "%d espèces créées, %d ignorées (id manquant)",
        len(id_map), skipped,
    )
    return id_map


def build_reactions(
    sbml_root,
    reactions: list[dict],
    id_map: dict[int, str],
) -> tuple[int, int]:
    """
    Ajoute toutes les réactions au document SBML.

    Returns:
        (nb_created, nb_skipped)
    """
    created = 0
    skipped = 0

    for rxn in reactions:
        rxn_el = make_reaction(rxn, id_map, module_prefix="")
        if rxn_el is None:
            skipped += 1
            continue
        add_reaction(sbml_root, rxn_el)
        created += 1

    logger.info(
        "%d réactions créées, %d ignorées (espèces non résolues)",
        created, skipped,
    )
    return created, skipped


# ---------------------------------------------------------------------------
# Validation et rapport
# ---------------------------------------------------------------------------

def _check_phenotype_nodes(elements: list[dict]) -> tuple[set[str], set[str]]:
    """
    Vérifie la présence des 14 nœuds phénotypiques attendus.

    Returns:
        (found_names, missing_names)
    """
    phenotypes = get_phenotype_nodes(elements)
    found_names = {p.get("name", "").strip() for p in phenotypes}

    # Correspondance souple (insensible à la casse, espace/underscore)
    def normalize(s: str) -> str:
        return s.lower().replace(" ", "_").replace("/", "_").replace("-", "_")

    found_norm = {normalize(n) for n in found_names}
    missing = {
        exp for exp in EXPECTED_PHENOTYPES
        if normalize(exp) not in found_norm
    }
    return found_names, missing


def print_report(
    elements: list[dict],
    reactions: list[dict],
    sbml_root,
    output_path: Path,
    valid: bool,
    sbml_errors: list[str],
) -> None:
    """Affiche un rapport de construction structuré."""
    print("\n" + "=" * 60)
    print("  RAPPORT — SjD_Map_original.xml")
    print("=" * 60)

    # Résumé API
    el_summary = summarize_elements(elements)
    rx_summary = summarize_reactions(reactions)

    print(f"\n[API MINERVA]")
    print(f"  Éléments récupérés   : {el_summary['total']}")
    print(f"  Réactions récupérées : {rx_summary['total']}")
    print(f"  Types d'éléments     : {el_summary['by_type']}")
    print(f"  Types de réactions   : {rx_summary['by_type']}")

    # Résumé XML
    counts = count_elements(sbml_root)
    print(f"\n[XML produit]")
    print(f"  Compartiments : {counts['compartments']}")
    print(f"  Espèces       : {counts['species']}")
    print(f"  Réactions     : {counts['reactions']}")
    print(f"  Taille fichier: {output_path.stat().st_size / 1024:.1f} Ko")

    # Nœuds phénotypiques
    found, missing = _check_phenotype_nodes(elements)
    print(f"\n[Nœuds Phenotype]")
    print(f"  Trouvés ({len(found)}) : {sorted(found)}")
    if missing:
        print(f"  MANQUANTS ({len(missing)}) : {sorted(missing)}")
    else:
        print("  Tous les phenotypes attendus sont présents.")

    # Validation SBML
    print(f"\n[Validation SBML]")
    if valid:
        print("  OK — aucune erreur fatale libsbml")
    else:
        print(f"  ECHEC — {len([e for e in sbml_errors if 'Error' in e or 'Fatal' in e])} erreur(s)")
        for err in sbml_errors[:20]:
            print(f"    {err}")
        if len(sbml_errors) > 20:
            print(f"    ... ({len(sbml_errors) - 20} erreurs supplémentaires)")

    # Vérification des références d'espèces
    unresolved = check_species_references(sbml_root)
    print(f"\n[Références espèces]")
    if unresolved:
        print(f"  {len(unresolved)} référence(s) non résolue(s) :")
        for u in unresolved[:10]:
            print(f"    {u}")
    else:
        print("  Toutes les speciesReference sont résolues.")

    # Gate 1 — critères de validation obligatoires
    print(f"\n[Gate 1 — Critères de passage]")
    gates = {
        "Species >= 800"    : counts["species"] >= 800,
        "Reactions >= 590"  : counts["reactions"] >= 590,
        "SBML valide"       : valid,
        "Refs résolues"     : len(unresolved) == 0,
        "14 phénotypes"     : len(missing) == 0,
        "Fichier > 1 Mo"    : output_path.stat().st_size > 1_000_000,
    }
    all_pass = True
    for label, ok in gates.items():
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {label}")
        if not ok:
            all_pass = False

    print()
    if all_pass:
        print("=> Gate 1 : VALIDÉ — prêt pour l'étape d'audit (02_parse_and_audit.py)")
    else:
        print("=> Gate 1 : ECHEC — corriger les points ci-dessus avant de continuer")
    print("=" * 60 + "\n")


# ---------------------------------------------------------------------------
# Point d'entrée
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Télécharge le SjD Map depuis MINERVA et reconstruit le XML CellDesigner."
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Ignore le cache JSON et ré-interroge l'API MINERVA",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "01_disease_map" / "SjD_Map_original.xml",
        help="Chemin du fichier XML de sortie",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=ROOT / "01_disease_map" / "cache",
        help="Répertoire de cache JSON",
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

    # 1. Récupération des données depuis MINERVA (ou cache)
    logger.info("=== Étape 1/5 : Récupération des données MINERVA ===")
    session = make_session()
    elements = load_or_fetch_elements(args.cache_dir, session, args.force_refresh)
    reactions = load_or_fetch_reactions(args.cache_dir, session, args.force_refresh)

    logger.info(
        "%d éléments et %d réactions chargés",
        len(elements), len(reactions),
    )

    # Dimensions du canvas (utiliser les bounds max des éléments)
    max_x = max((e.get("bounds", {}) or {}).get("x", 0) + (e.get("bounds", {}) or {}).get("w", 60)
                for e in elements if e.get("bounds"))
    max_y = max((e.get("bounds", {}) or {}).get("y", 0) + (e.get("bounds", {}) or {}).get("h", 25)
                for e in elements if e.get("bounds"))
    width = max(max_x * 1.1, 10815.0)
    height = max(max_y * 1.1, 6000.0)
    logger.info("Canvas : %.0f × %.0f px", width, height)

    # 2. Création de la racine SBML
    logger.info("=== Étape 2/5 : Création de la structure SBML ===")
    sbml_root = make_sbml_root(
        model_id="SjD_Map",
        model_name="Sjogren's Disease Map",
        width=width,
        height=height,
    )

    # 3. Compartiments
    logger.info("=== Étape 3/5 : Construction des compartiments ===")
    comp_map = build_compartments(sbml_root, elements)

    # 4. Espèces
    logger.info("=== Étape 4/5 : Construction des espèces ===")
    id_map = build_species(sbml_root, elements, comp_map)

    # 5. Réactions
    logger.info("=== Étape 5/5 : Construction des réactions ===")
    n_rxn, n_skip = build_reactions(sbml_root, reactions, id_map)

    if n_skip > 0:
        logger.warning(
            "%d réaction(s) ignorée(s) — espèces source/cible non résolues. "
            "Vérifier les éléments avec id manquant dans le cache.",
            n_skip,
        )

    # Écriture du fichier
    write_celldesigner_xml(sbml_root, args.output)

    # Validation SBML
    valid, sbml_errors = validate_sbml(args.output)

    # Rapport final
    print_report(elements, reactions, sbml_root, args.output, valid, sbml_errors)

    return 0 if valid else 1


if __name__ == "__main__":
    sys.exit(main())
