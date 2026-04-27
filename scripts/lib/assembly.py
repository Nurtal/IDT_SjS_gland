"""
assembly.py — Assemblage final du multi-cellular CellDesigner SBML map (Phase 1.5).

Fusionne les 10 modules cell-type extraits en Phase 1.3 + les 470 edges
inter-cellulaires curés en Phase 1.4 en un unique fichier SBML/CellDesigner
compatible avec CaSQ et BMA pour la conversion booléenne (Phase 2.1).

Stratégie :
    1. Compartiments :
        - 1 compartiment intracellulaire par (cell-type × MINERVA-compartment)
          ex. SGEC_Cytoplasm, SGEC_Nucleus, SGEC_EndoplasmicReticulum, SGEC_Default,
              TH1_Cytoplasm, … (40 au total).
        - 3 compartiments partagés : Extracellular (21555), Secreted (21629),
          Phenotypes (21540).

    2. Espèces :
        - Pour chaque (node_id, celltype) ∈ node_to_celltype.tsv, celltype ∈ 10 réels :
          clone l'espèce avec préfixe `<CT>_s<id>` dans le compartiment intracellulaire.
        - Pour chaque node_id ∈ EXTRA : 1 instance unique partagée.
        - Pour chaque node_id ∈ PHENOTYPE : 1 instance unique partagée.

    3. Réactions intracellulaires :
        - Pour chaque réaction MINERVA, pour chaque cell-type C : émettre une
          copie préfixée si tous les participants sont dans (core_C ∪ EXTRA ∪
          PHENOTYPE) ET au moins 1 ∈ core_C.

    4. Réactions inter-cellulaires (Zerrouk 2024) :
        - mechanism="secreted" / "autocrine" :
            EXTRA[ligand] → PHYSICAL_STIMULATION → <tgt_ct>_<receptor>
        - mechanism="contact" :
            EXTRA[ligand] + <tgt_ct>_<receptor> → HETERODIMER_ASSOCIATION
            → complex synthétique en Extracellular.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from lxml import etree

from lib.celldesigner_xml import (  # noqa: E402
    MINERVA_TYPE_TO_CD,
    add_compartment,
    add_reaction,
    add_species,
    add_species_alias,
    make_extracellular_species,
    make_heterodimer_reaction,
    make_physical_stimulation_reaction,
    make_reaction,
    make_sbml_root,
    make_species,
    _compartment_id,
    _safe_id,
)

logger = logging.getLogger(__name__)


CELL_TYPES: tuple[str, ...] = (
    "SGEC", "TH1", "TH17", "TFH", "TREG",
    "BCELL", "PLASMA", "M1", "M2", "PDC",
)

# Compartiments MINERVA
COMP_CYTOPLASM = 20513
COMP_NUCLEUS = 21231
COMP_ER = 20730
COMP_EXTRACELLULAR = 21555
COMP_SECRETED = 21629
COMP_PHENOTYPES = 21540
INTRACELLULAR_COMPS: tuple[int | None, ...] = (
    COMP_CYTOPLASM, COMP_NUCLEUS, COMP_ER, None,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _resolve_id(entry: dict[str, Any]) -> int | None:
    eid = entry.get("aliasId")
    if eid is None:
        elem = entry.get("element") or {}
        eid = elem.get("id")
    return int(eid) if eid is not None else None


def _reaction_participants(rxn: dict[str, Any]) -> set[int]:
    out: set[int] = set()
    for key in ("reactants", "products", "modifiers"):
        for entry in rxn.get(key) or []:
            rid = _resolve_id(entry)
            if rid is not None:
                out.add(rid)
    return out


def _intracell_compartment(elem: dict[str, Any]) -> int | None:
    """Choisit le compartiment intracellulaire cible pour un élément core."""
    cid = elem.get("compartmentId")
    if cid in (COMP_CYTOPLASM, COMP_NUCLEUS, COMP_ER):
        return cid
    return None  # → "Default"


# ---------------------------------------------------------------------------
# Contexte d'assemblage
# ---------------------------------------------------------------------------


@dataclass
class AssemblyStats:
    n_compartments: int = 0
    n_species_core: int = 0
    n_species_extra: int = 0
    n_species_phenotype: int = 0
    n_species_complex: int = 0
    n_reactions_intra: int = 0
    n_reactions_skipped: int = 0
    n_inter_secreted: int = 0
    n_inter_contact: int = 0
    n_inter_skipped: int = 0
    by_celltype_reactions: dict[str, int] = field(default_factory=dict)


@dataclass
class AssemblyContext:
    n2c_rows: list[dict[str, Any]]
    elements_by_id: dict[int, dict[str, Any]]
    reactions: list[dict[str, Any]]
    intercellular_rows: list[dict[str, Any]]

    core_by_ct: dict[str, set[int]] = field(default_factory=dict)
    extra_ids: set[int] = field(default_factory=set)
    phenotype_ids: set[int] = field(default_factory=set)
    id_map: dict[tuple[str, int], str] = field(default_factory=dict)
    extra_id_map: dict[int, str] = field(default_factory=dict)
    pheno_id_map: dict[int, str] = field(default_factory=dict)

    def index_assignments(self) -> None:
        self.core_by_ct = {ct: set() for ct in CELL_TYPES}
        for row in self.n2c_rows:
            ct = (row.get("celltype") or "").strip()
            try:
                nid = int(row["node_id"])
            except (KeyError, ValueError, TypeError):
                continue
            if ct in self.core_by_ct:
                self.core_by_ct[ct].add(nid)
            elif ct == "EXTRA":
                self.extra_ids.add(nid)
            elif ct == "PHENOTYPE":
                self.phenotype_ids.add(nid)


# ---------------------------------------------------------------------------
# Assemblage
# ---------------------------------------------------------------------------


def _build_compartments(
    sbml: etree._Element,
    stats: AssemblyStats,
) -> None:
    """Compartiments per-cell-type + 3 partagés."""
    for ct in CELL_TYPES:
        for cid in INTRACELLULAR_COMPS:
            add_compartment(sbml, cid, module_prefix=ct)
            stats.n_compartments += 1
    for cid in (COMP_EXTRACELLULAR, COMP_SECRETED, COMP_PHENOTYPES):
        add_compartment(sbml, cid, module_prefix="")
        stats.n_compartments += 1


def _build_species(
    sbml: etree._Element,
    ctx: AssemblyContext,
    stats: AssemblyStats,
) -> None:
    extra_comp_id = _compartment_id(COMP_EXTRACELLULAR, "")
    secreted_comp_id = _compartment_id(COMP_SECRETED, "")
    pheno_comp_id = _compartment_id(COMP_PHENOTYPES, "")

    def _emit(elem: dict[str, Any], comp_id: str, prefix: str) -> str:
        species = make_species(elem, comp_id, module_prefix=prefix)
        add_species(sbml, species)
        sid = species.get("id")
        cd_class = MINERVA_TYPE_TO_CD.get(elem.get("type", ""), "PROTEIN")
        add_species_alias(
            sbml, sid, comp_id,
            bounds=elem.get("bounds") or {},
            species_class=cd_class,
        )
        return sid

    # 1) Espèces core par cell-type
    for ct in CELL_TYPES:
        for nid in sorted(ctx.core_by_ct.get(ct, set())):
            elem = ctx.elements_by_id.get(nid)
            if not elem:
                continue
            cid = _intracell_compartment(elem)
            comp_sbml_id = _compartment_id(cid, ct)
            ctx.id_map[(ct, nid)] = _emit(elem, comp_sbml_id, ct)
            stats.n_species_core += 1

    # 2) EXTRA partagés
    for nid in sorted(ctx.extra_ids):
        elem = ctx.elements_by_id.get(nid)
        if not elem:
            continue
        comp_id = (
            secreted_comp_id
            if elem.get("compartmentId") == COMP_SECRETED
            else extra_comp_id
        )
        ctx.extra_id_map[nid] = _emit(elem, comp_id, "")
        stats.n_species_extra += 1

    # 3) PHENOTYPE partagés
    for nid in sorted(ctx.phenotype_ids):
        elem = ctx.elements_by_id.get(nid)
        if not elem:
            continue
        ctx.pheno_id_map[nid] = _emit(elem, pheno_comp_id, "")
        stats.n_species_phenotype += 1


def _build_intracellular_reactions(
    sbml: etree._Element,
    ctx: AssemblyContext,
    stats: AssemblyStats,
) -> None:
    """Réplique chaque réaction MINERVA pour chaque cell-type compatible."""
    stats.by_celltype_reactions = {ct: 0 for ct in CELL_TYPES}
    for rxn in ctx.reactions:
        participants = _reaction_participants(rxn)
        if not participants:
            stats.n_reactions_skipped += 1
            continue
        emitted_any = False
        for ct in CELL_TYPES:
            local_id_map: dict[int, str] = {}
            ok = True
            has_core = False
            for nid in participants:
                if (ct, nid) in ctx.id_map:
                    local_id_map[nid] = ctx.id_map[(ct, nid)]
                    has_core = True
                elif nid in ctx.extra_id_map:
                    local_id_map[nid] = ctx.extra_id_map[nid]
                elif nid in ctx.pheno_id_map:
                    local_id_map[nid] = ctx.pheno_id_map[nid]
                else:
                    ok = False
                    break
            if not ok or not has_core:
                continue
            reaction = make_reaction(rxn, local_id_map, module_prefix=ct)
            if reaction is None:
                continue
            add_reaction(sbml, reaction)
            stats.n_reactions_intra += 1
            stats.by_celltype_reactions[ct] += 1
            emitted_any = True
        if not emitted_any:
            stats.n_reactions_skipped += 1


def _parse_int_list(raw: str) -> list[int]:
    out: list[int] = []
    for tok in (raw or "").split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.append(int(tok))
        except ValueError:
            continue
    return out


def _build_intercellular_edges(
    sbml: etree._Element,
    ctx: AssemblyContext,
    stats: AssemblyStats,
) -> None:
    extra_comp_id = _compartment_id(COMP_EXTRACELLULAR, "")
    for i, edge in enumerate(ctx.intercellular_rows):
        ligand_ids = _parse_int_list(edge.get("ligand_node_ids", ""))
        receptor_ids = _parse_int_list(edge.get("receptor_node_ids", ""))
        if not ligand_ids or not receptor_ids:
            stats.n_inter_skipped += 1
            continue
        ligand_nid = ligand_ids[0]
        receptor_nid = receptor_ids[0]
        src_ct = edge.get("source_celltype", "")
        tgt_ct = edge.get("target_celltype", "")
        mechanism = edge.get("mechanism", "secreted")
        ligand = edge.get("ligand", "L")
        receptor = edge.get("receptor", "R")

        src_sbml = ctx.extra_id_map.get(ligand_nid)
        tgt_sbml = ctx.id_map.get((tgt_ct, receptor_nid))
        if not src_sbml or not tgt_sbml:
            stats.n_inter_skipped += 1
            continue

        rid = _safe_id(f"inter_{src_ct}_{tgt_ct}_{ligand}_{receptor}_{i}")
        rname = f"{ligand}({src_ct})->{receptor}({tgt_ct})"

        if mechanism == "contact":
            complex_id = _safe_id(
                f"complex_{src_ct}_{tgt_ct}_{ligand}_{receptor}_{i}"
            )
            cx_species = make_extracellular_species(
                f"{ligand}_{receptor}_complex",
                extra_comp_id,
                species_id=complex_id,
            )
            add_species(sbml, cx_species)
            add_species_alias(
                sbml, complex_id, extra_comp_id,
                bounds={}, species_class="COMPLEX",
            )
            stats.n_species_complex += 1
            rxn = make_heterodimer_reaction(
                rid, src_sbml, tgt_sbml, complex_id,
                reaction_name=rname,
            )
            add_reaction(sbml, rxn)
            stats.n_inter_contact += 1
        else:  # secreted | autocrine
            rxn = make_physical_stimulation_reaction(
                rid, src_sbml, tgt_sbml, reaction_name=rname,
            )
            add_reaction(sbml, rxn)
            stats.n_inter_secreted += 1


def assemble_map(ctx: AssemblyContext) -> tuple[etree._Element, AssemblyStats]:
    """Construit l'arbre SBML/CellDesigner complet et renvoie (root, stats)."""
    stats = AssemblyStats()
    ctx.index_assignments()

    sbml = make_sbml_root(
        "SjD_multicellular",
        "SjD multicellular Boolean map",
    )

    logger.info("Compartiments…")
    _build_compartments(sbml, stats)

    logger.info("Espèces (core + EXTRA + PHENOTYPE)…")
    _build_species(sbml, ctx, stats)
    logger.info(
        "  → core=%d  extra=%d  pheno=%d",
        stats.n_species_core, stats.n_species_extra, stats.n_species_phenotype,
    )

    logger.info("Réactions intracellulaires…")
    _build_intracellular_reactions(sbml, ctx, stats)
    logger.info(
        "  → intracell=%d  skipped=%d",
        stats.n_reactions_intra, stats.n_reactions_skipped,
    )

    logger.info("Edges inter-cellulaires…")
    _build_intercellular_edges(sbml, ctx, stats)
    logger.info(
        "  → secreted/autocrine=%d  contact=%d  skipped=%d",
        stats.n_inter_secreted, stats.n_inter_contact, stats.n_inter_skipped,
    )

    return sbml, stats
