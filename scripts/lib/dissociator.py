"""
dissociator.py — Assignation des nœuds de la SjD Map à des cell-types.

Implémente les règles R1–R7 spécifiées dans
`01_disease_map/dissociation_rules.md` (Phase 1.2 ROADMAP).

Approche hybride à priorité décroissante :
    R1  Compartiment extracellulaire / sécrétoire → EXTRA      [HIGH]
    R2  Type Phenotype                            → PHENOTYPE  [HIGH]
    R3  Marqueur cell-type exclusif               → mono / lignée [HIGH]
    R4  Pathway-driven (Reactome / KEGG / GO)     → multi      [MEDIUM]
    R5  Propagation par voisinage (random walk)   → multi      [LOW]
    R6  Default fallback (signaling générique)    → 6 cell-types [LOW]
    R7  Inassignable                              → UNASSIGNED

L'algorithme procède en 2 passes :
    Passe 1 : R1, R2, R3, R4 (règles autonomes par nœud)
    Passe 2 : R5 (voisinage), R6 (fallback), R7 (résiduel)

Multi-assignment autorisé : un même nœud peut être cloné dans plusieurs
cell-types (ex. STAT1 → STAT1_TH1, STAT1_SGEC, ...).

Sortie : dict {node_id: Assignment} où Assignment contient celltypes,
confidence, rule, evidence pour audit.
"""

from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any

import networkx as nx

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Cell-types et marqueurs (cf. dissociation_rules.md §3 R3)
# ---------------------------------------------------------------------------

CELL_TYPES: tuple[str, ...] = (
    "SGEC", "TH1", "TH17", "TFH", "TREG",
    "BCELL", "PLASMA", "M1", "M2", "PDC",
)

# Pseudo-cell-types
EXTRA = "EXTRA"
PHENOTYPE = "PHENOTYPE"
UNASSIGNED = "UNASSIGNED"

# R1 : compartiments extracellulaires / sécrétoires
EXTRACELLULAR_COMPARTMENTS: set[int] = {21555, 21629}

# R3 : marqueurs exclusifs par lignée → set de cell-types
EXCLUSIVE_MARKERS: dict[str, set[str]] = {
    # SGEC-only
    "AQP5": {"SGEC"}, "AQP3": {"SGEC"}, "MUC5B": {"SGEC"}, "MUC7": {"SGEC"},
    "KRT7": {"SGEC"}, "KRT8": {"SGEC"}, "KRT18": {"SGEC"}, "KRT19": {"SGEC"},
    "EPCAM": {"SGEC"}, "CFTR": {"SGEC"}, "CDH1": {"SGEC"},
    "CLDN1": {"SGEC"}, "CLDN3": {"SGEC"}, "CLDN4": {"SGEC"}, "OCLN": {"SGEC"},
    "ZO1": {"SGEC"}, "TJP1": {"SGEC"},
    "AMY1A": {"SGEC"}, "PRB1": {"SGEC"}, "STATH": {"SGEC"},
    # T-lineage commun (TCR signaling)
    "CD3D": {"TH1", "TH17", "TFH", "TREG"},
    "CD3E": {"TH1", "TH17", "TFH", "TREG"},
    "CD3G": {"TH1", "TH17", "TFH", "TREG"},
    "CD4": {"TH1", "TH17", "TFH", "TREG"},
    "CD2": {"TH1", "TH17", "TFH", "TREG"},
    "LCK": {"TH1", "TH17", "TFH", "TREG"},
    "ZAP70": {"TH1", "TH17", "TFH", "TREG"},
    "LAT": {"TH1", "TH17", "TFH", "TREG"},
    "LCP2": {"TH1", "TH17", "TFH", "TREG"},
    "CD28": {"TH1", "TH17", "TFH", "TREG"},
    # TH1-only (les cytokines IFNG sont surtout R1 si extracellulaires)
    "TBX21": {"TH1"}, "STAT4": {"TH1"},
    # TH17-only
    "RORC": {"TH17"}, "RORA": {"TH17"}, "IL23R": {"TH17"}, "AHR": {"TH17"},
    # TFH-only
    "CXCR5": {"TFH"}, "BCL6": {"TFH"}, "ICOS": {"TFH"}, "PDCD1": {"TFH"},
    # TREG-only
    "FOXP3": {"TREG"}, "IL2RA": {"TREG"}, "IKZF2": {"TREG"}, "CTLA4": {"TREG"},
    # B-lineage (BCELL + PLASMA)
    "CD19": {"BCELL", "PLASMA"},
    "MS4A1": {"BCELL"},
    "CD79A": {"BCELL", "PLASMA"},
    "CD79B": {"BCELL", "PLASMA"},
    "CR2": {"BCELL"},
    "BLNK": {"BCELL"},
    "BANK1": {"BCELL"},
    # BCELL-only
    "BTK": {"BCELL"}, "BLK": {"BCELL"},
    "TNFRSF13B": {"BCELL"},  # TACI
    "TNFRSF13C": {"BCELL"},  # BAFFR
    "FCRL3": {"BCELL"},
    # PLASMA-only
    "XBP1": {"PLASMA"}, "IRF4": {"PLASMA"}, "PRDM1": {"PLASMA"},
    "SDC1": {"PLASMA"}, "MZB1": {"PLASMA"}, "JCHAIN": {"PLASMA"},
    "TNFRSF17": {"PLASMA"},  # BCMA
    # Macrophage-lineage (M1 + M2)
    "CD68": {"M1", "M2"}, "FCGR1A": {"M1", "M2"}, "FCGR2A": {"M1", "M2"},
    "FCGR3A": {"M1", "M2"}, "CSF1R": {"M1", "M2"}, "MERTK": {"M1", "M2"},
    "ITGAM": {"M1", "M2"},
    # M1-only
    "NOS2": {"M1"}, "CCR7": {"M1"},
    # M2-only
    "CD163": {"M2"}, "MRC1": {"M2"}, "ARG1": {"M2"}, "MSR1": {"M2"},
    # PDC-only
    "CLEC4C": {"PDC"}, "IL3RA": {"PDC"}, "LILRA4": {"PDC"},
    "IRF7": {"PDC"}, "TLR7": {"PDC"}, "TLR9": {"PDC"},
}

# R4 : table pathway → cell-types
ALL_LYMPHOCYTES = {"TH1", "TH17", "TFH", "TREG", "BCELL", "PLASMA"}
ALL_T = {"TH1", "TH17", "TFH", "TREG"}
ALL_INNATE = {"M1", "M2", "PDC", "SGEC"}
ALL_CT = set(CELL_TYPES)

PATHWAY_TO_CELLTYPES: dict[str, set[str]] = {
    # Reactome
    "R-HSA-1280215": ALL_CT,                           # Cytokine signaling
    "R-HSA-877300": ALL_CT,                            # IFN signaling
    "R-HSA-877253": {"SGEC", "PDC", "BCELL", "M1", "M2"},  # IFN α/β
    "R-HSA-877312": {"SGEC", "M1", "M2"},              # IFN γ
    "R-HSA-983705": {"BCELL", "PLASMA"},               # BCR signaling
    "R-HSA-202403": ALL_T,                             # TCR signaling
    "R-HSA-168256": ALL_CT,                            # Immune system / Adaptive
    "R-HSA-168249": ALL_INNATE,                        # Innate immune system
    "R-HSA-168928": ALL_INNATE,                        # DDX58/IFIH1 RNA sensing
    "R-HSA-449147": ALL_INNATE,                        # TLR signaling
    "R-HSA-983168": ALL_CT - {"TREG"},                 # Antigen processing/MHC-I
    "R-HSA-2132295": {"M1", "M2", "BCELL", "PDC", "SGEC"},  # MHC-II
    "R-HSA-1474244": {"SGEC"},                         # ECM organisation
    "R-HSA-372790": ALL_CT,                            # Signaling by GPCR
    "R-HSA-451927": ALL_T,                             # IL-2 family
    "R-HSA-512988": ALL_CT,                            # IL-6 signaling
    "R-HSA-9020702": {"M1", "M2", "PDC"},              # Interleukin-1 family
    "R-HSA-446652": ALL_CT,                            # Interleukins
    "R-HSA-2871837": {"BCELL", "PLASMA"},              # FCERI signaling
    # KEGG (préfixe hsa)
    "hsa04060": ALL_CT,   # cytokine-cytokine receptor
    "hsa04062": ALL_CT,   # chemokine signaling
    "hsa04630": ALL_CT,   # JAK-STAT
    "hsa04064": ALL_CT,   # NF-kB
    "hsa04668": ALL_CT,   # TNF
    "hsa04620": ALL_INNATE,  # TLR signaling
    "hsa04660": ALL_T,    # TCR
    "hsa04662": {"BCELL", "PLASMA"},  # BCR
    "hsa04640": ALL_CT,   # Hematopoietic cell lineage
    "hsa04612": ALL_CT,   # Antigen processing/presentation
    "hsa04658": {"TH1"},  # Th1/Th2 differentiation
    "hsa04659": {"TH17"}, # Th17 differentiation
    "hsa04672": {"BCELL", "PLASMA"},  # Intestinal IgA production
    "hsa04621": ALL_INNATE,  # NOD-like receptor
    "hsa04622": ALL_INNATE,  # RIG-I-like receptor
    "hsa04610": {"M1", "M2"},  # Complement
}

# R6 : fallback default — 6 cell-types principaux
DEFAULT_FALLBACK_CELLTYPES: set[str] = {
    "SGEC", "TH1", "TH17", "BCELL", "M1", "M2",
}

# Confiance ordonnée
CONFIDENCE_RANK: dict[str, int] = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}

# Patterns pathway depuis les notes (memo : déjà dans map_audit)
PATHWAY_PATTERNS: list[re.Pattern] = [
    re.compile(r"R-HSA-\d+", re.I),
    re.compile(r"hsa\d{5}", re.I),
    re.compile(r"GO:\d{7}"),
]


# ---------------------------------------------------------------------------
# Structure d'assignation
# ---------------------------------------------------------------------------


@dataclass
class Assignment:
    """Résultat d'assignation pour un nœud."""
    node_id: int
    name: str
    type: str
    compartment_id: int | None
    celltypes: set[str] = field(default_factory=set)
    confidence: str = ""        # HIGH / MEDIUM / LOW
    rule: str = ""              # R1..R7
    evidence: str = ""          # texte court (marker, pathway id, voisin...)

    def to_rows(self) -> list[dict[str, Any]]:
        """Une ligne par cell-type assigné (format long)."""
        if not self.celltypes:
            return [{
                "node_id": self.node_id,
                "node_name": self.name,
                "node_type": self.type,
                "compartment_id": self.compartment_id,
                "celltype": UNASSIGNED,
                "confidence": "",
                "rule": "R7",
                "evidence": "",
            }]
        return [{
            "node_id": self.node_id,
            "node_name": self.name,
            "node_type": self.type,
            "compartment_id": self.compartment_id,
            "celltype": ct,
            "confidence": self.confidence,
            "rule": self.rule,
            "evidence": self.evidence,
        } for ct in sorted(self.celltypes)]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _normalise_name(name: str) -> str:
    """Normalisation pour matching marqueur : majuscules, alphanum only."""
    return re.sub(r"[^A-Z0-9]", "", (name or "").upper())


def _extract_pathway_ids(notes: str) -> list[str]:
    """Extrait les identifiants Reactome/KEGG/GO depuis le texte des notes."""
    if not notes:
        return []
    out: list[str] = []
    for pat in PATHWAY_PATTERNS:
        out.extend(pat.findall(notes))
    return out


# ---------------------------------------------------------------------------
# Règles R1–R4 (passe 1)
# ---------------------------------------------------------------------------


def rule_R1_extracellular(elem: dict[str, Any]) -> Assignment | None:
    cid = elem.get("compartmentId")
    if cid is None or int(cid) not in EXTRACELLULAR_COMPARTMENTS:
        return None
    return Assignment(
        node_id=int(elem["id"]),
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=int(cid),
        celltypes={EXTRA},
        confidence="HIGH",
        rule="R1",
        evidence=f"compartment={cid}",
    )


def rule_R2_phenotype(elem: dict[str, Any]) -> Assignment | None:
    if (elem.get("type") or "") != "Phenotype":
        return None
    return Assignment(
        node_id=int(elem["id"]),
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes={PHENOTYPE},
        confidence="HIGH",
        rule="R2",
        evidence="type=Phenotype",
    )


def rule_R3_marker(elem: dict[str, Any]) -> Assignment | None:
    raw = elem.get("name") or ""
    norm = _normalise_name(raw)
    upper = raw.upper().strip()
    # Match exact ou alphanum-only ; teste aussi le nom brut majuscules
    for candidate in (upper, norm):
        if candidate in EXCLUSIVE_MARKERS:
            cts = EXCLUSIVE_MARKERS[candidate]
            return Assignment(
                node_id=int(elem["id"]),
                name=raw,
                type=elem.get("type") or "",
                compartment_id=elem.get("compartmentId"),
                celltypes=set(cts),
                confidence="HIGH",
                rule="R3",
                evidence=f"marker={candidate}",
            )
    return None


def rule_R4_pathway(elem: dict[str, Any]) -> Assignment | None:
    notes = elem.get("notes") or ""
    pids = _extract_pathway_ids(notes)
    if not pids:
        return None
    cts: set[str] = set()
    matched_ids: list[str] = []
    for pid in pids:
        # Recherche d'abord match exact, sinon préfixe (R-HSA-1280215...)
        if pid in PATHWAY_TO_CELLTYPES:
            cts |= PATHWAY_TO_CELLTYPES[pid]
            matched_ids.append(pid)
            continue
        for known, k_cts in PATHWAY_TO_CELLTYPES.items():
            if pid.startswith(known) or known.startswith(pid):
                cts |= k_cts
                matched_ids.append(pid)
                break
    if not cts:
        return None
    return Assignment(
        node_id=int(elem["id"]),
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes=cts,
        confidence="MEDIUM",
        rule="R4",
        evidence=";".join(sorted(set(matched_ids))[:5]),
    )


# ---------------------------------------------------------------------------
# Règles R5–R7 (passe 2)
# ---------------------------------------------------------------------------


def rule_R5_neighbor(
    elem: dict[str, Any],
    graph: nx.DiGraph,
    assignments: dict[int, Assignment],
    min_concordant: int = 2,
) -> Assignment | None:
    nid = int(elem["id"])
    if not graph.has_node(nid):
        return None
    # Voisinage non orienté
    neighbors = set(graph.predecessors(nid)) | set(graph.successors(nid))
    if not neighbors:
        return None

    # Collecte des cell-types des voisins déjà assignés (R1 EXTRA et R2 PHENOTYPE
    # ne propagent pas leur assignation — ce sont des pseudo-cell-types)
    votes: Counter = Counter()
    for nb in neighbors:
        a = assignments.get(nb)
        if a is None or not a.celltypes:
            continue
        if a.celltypes & {EXTRA, PHENOTYPE}:
            continue
        weight = CONFIDENCE_RANK.get(a.confidence, 1)
        for ct in a.celltypes:
            if ct in ALL_CT:
                votes[ct] += weight

    if not votes:
        return None

    # Sélectionner les cell-types ayant ≥ min_concordant voisins concordants
    concordant = [ct for ct, c in votes.items() if c >= min_concordant]
    if not concordant:
        return None

    return Assignment(
        node_id=nid,
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes=set(concordant),
        confidence="LOW",
        rule="R5",
        evidence=f"neighbors_vote={dict(votes.most_common(5))}",
    )


INTRACELLULAR_COMPARTMENTS: set[int] = {20513, 21231, 20730}  # Cell / nucleus / ER


def rule_R6_fallback(
    elem: dict[str, Any],
    graph: nx.DiGraph,
) -> Assignment | None:
    etype = elem.get("type") or ""
    if etype not in {"Protein", "Complex", "Simple molecule", "SimpleMolecule", "RNA", "Gene", "Ion"}:
        return None
    nid = int(elem["id"])
    cid = elem.get("compartmentId")
    is_intracellular = cid is not None and int(cid) in INTRACELLULAR_COMPARTMENTS

    has_edges = False
    if graph.has_node(nid):
        has_edges = (graph.in_degree(nid) + graph.out_degree(nid)) > 0

    # Connecté OU intracellulaire (signaling présumé partagé) déclenche R6.
    # Reste R7 : isolés sans compartiment intracellulaire (satellites incertains).
    if not has_edges and not is_intracellular:
        return None

    evidence = "default_signaling_fallback"
    if not has_edges:
        evidence += f"_isolated_in_compartment_{cid}"

    return Assignment(
        node_id=nid,
        name=elem.get("name") or "",
        type=etype,
        compartment_id=cid,
        celltypes=set(DEFAULT_FALLBACK_CELLTYPES),
        confidence="LOW",
        rule="R6",
        evidence=evidence,
    )


def rule_R7_unassigned(elem: dict[str, Any]) -> Assignment:
    return Assignment(
        node_id=int(elem["id"]),
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes=set(),
        confidence="",
        rule="R7",
        evidence="no_rule_matched",
    )


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def dissociate(
    elements: list[dict[str, Any]],
    graph: nx.DiGraph,
) -> dict[int, Assignment]:
    """
    Applique les règles R1–R7 dans l'ordre et renvoie {node_id: Assignment}.

    Args:
        elements : éléments MINERVA (depuis cache JSON).
        graph    : graphe networkx (depuis lib.map_audit.build_graph).

    Returns:
        dict {node_id (int): Assignment}
    """
    assignments: dict[int, Assignment] = {}
    pending: list[dict[str, Any]] = []

    # --- Passe 1 : R1, R2, R3, R4 ---
    for el in elements:
        if el.get("id") is None:
            continue
        a = (
            rule_R1_extracellular(el)
            or rule_R2_phenotype(el)
            or rule_R3_marker(el)
            or rule_R4_pathway(el)
        )
        if a is not None:
            assignments[a.node_id] = a
        else:
            pending.append(el)

    logger.info(
        "Passe 1 — R1:%d, R2:%d, R3:%d, R4:%d  (pending=%d)",
        sum(1 for a in assignments.values() if a.rule == "R1"),
        sum(1 for a in assignments.values() if a.rule == "R2"),
        sum(1 for a in assignments.values() if a.rule == "R3"),
        sum(1 for a in assignments.values() if a.rule == "R4"),
        len(pending),
    )

    # --- Passe 2 : R5 itératif (jusqu'à stabilité), puis R6, puis R7 ---
    still_pending = pending
    iteration = 0
    while True:
        iteration += 1
        new_assignments: list[Assignment] = []
        next_pending: list[dict[str, Any]] = []
        for el in still_pending:
            a = rule_R5_neighbor(el, graph, assignments)
            if a is not None:
                new_assignments.append(a)
            else:
                next_pending.append(el)
        for a in new_assignments:
            assignments[a.node_id] = a
        logger.info(
            "Passe R5 itération %d : +%d assignés, reste %d",
            iteration, len(new_assignments), len(next_pending),
        )
        if not new_assignments or iteration >= 5:
            break
        still_pending = next_pending

    # R6 fallback
    after_r6: list[dict[str, Any]] = []
    for el in still_pending:
        a = rule_R6_fallback(el, graph)
        if a is not None:
            assignments[a.node_id] = a
        else:
            after_r6.append(el)
    logger.info(
        "R6 fallback : +%d assignés, reste %d",
        sum(1 for a in assignments.values() if a.rule == "R6"),
        len(after_r6),
    )

    # R7 inassignables
    for el in after_r6:
        a = rule_R7_unassigned(el)
        assignments[a.node_id] = a

    return assignments


# ---------------------------------------------------------------------------
# Synthèse
# ---------------------------------------------------------------------------


def summarize(assignments: dict[int, Assignment]) -> dict[str, Any]:
    """Statistiques pour rapport."""
    n_total = len(assignments)
    by_rule: Counter = Counter(a.rule for a in assignments.values())
    by_confidence: Counter = Counter(a.confidence for a in assignments.values())

    by_celltype: Counter = Counter()
    for a in assignments.values():
        for ct in a.celltypes:
            by_celltype[ct] += 1

    n_extra = by_celltype.get(EXTRA, 0)
    n_phen = by_celltype.get(PHENOTYPE, 0)
    n_unassigned = sum(1 for a in assignments.values() if not a.celltypes)
    n_assignable = n_total - n_extra - n_phen
    n_assigned_real = sum(
        1 for a in assignments.values()
        if a.celltypes and not a.celltypes & {EXTRA, PHENOTYPE}
    )

    coverage_pct = round(
        100 * n_assigned_real / n_assignable, 2
    ) if n_assignable else 0.0

    # Cell-types respectant le seuil de 80 nœuds
    threshold = 80
    underprovisioned = [
        ct for ct in CELL_TYPES if by_celltype.get(ct, 0) < threshold
    ]

    return {
        "n_total_nodes": n_total,
        "n_extracellular": n_extra,
        "n_phenotype": n_phen,
        "n_unassigned": n_unassigned,
        "n_assignable": n_assignable,
        "n_assigned_to_celltype": n_assigned_real,
        "coverage_pct": coverage_pct,
        "by_rule": dict(by_rule.most_common()),
        "by_confidence": dict(by_confidence.most_common()),
        "by_celltype": dict(by_celltype.most_common()),
        "celltype_threshold": threshold,
        "underprovisioned_celltypes": underprovisioned,
        "gate_coverage_pass": coverage_pct >= 90.0,
        "gate_celltype_size_pass": len(underprovisioned) == 0,
    }


def collect_extracellular(
    assignments: dict[int, Assignment],
    elements: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Pour chaque nœud R1-EXTRA, collecte ses variantes (G/R/P) par nom."""
    by_name: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for el in elements:
        nm = (el.get("name") or "").strip().upper()
        if nm:
            by_name[nm].append(el)

    out: list[dict[str, Any]] = []
    for nid, a in assignments.items():
        if a.rule != "R1":
            continue
        nm = (a.name or "").strip().upper()
        siblings = [
            f"{e['id']}({e.get('type','?')})"
            for e in by_name.get(nm, [])
            if int(e["id"]) != nid
        ]
        out.append({
            "node_id": nid,
            "node_name": a.name,
            "node_type": a.type,
            "compartment_id": a.compartment_id,
            "canonical_form": nm,
            "g_r_p_variants": ",".join(siblings),
        })
    return out
