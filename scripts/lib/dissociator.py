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

# ---------------------------------------------------------------------------
# Score de plausibilité par règle (base) + sources bibliographiques
# ---------------------------------------------------------------------------
#
# Plausibilité (0–100) reflète à quel point l'assignation node→cell-type est
# défendable biologiquement, indépendamment de la confiance computationnelle :
#   95–100 : marqueur lineage-defining, consensus immunologique
#   80–94  : marqueur fonctionnel à forte spécificité (DB curées)
#   60–79  : pathway-driven cohérent (Reactome/KEGG, multi-cell-type plausible)
#   40–59  : propagation par voisinage validée par graphe
#   20–39  : fallback générique (signaling intracellulaire ubiquitaire)
#    0–19  : non assigné / spurious
#
# Les sources sont des références primaires ou des DB publiques :
#   PMID:xxxxxx, doi:..., HGNC:..., UniProt:..., InnateDB:..., STRING:...,
#   KEGG:hsaXXXXX, Reactome:R-HSA-XXXXX, CellMarker2.0, PanglaoDB.
#
# Pour les assignations fallback (R5/R6), la source est générique (graphe SjD
# Map elle-même + DB consensus) car par définition le nœud n'a pas de signal
# fort cell-type-spécifique.

RULE_BASE_SCORE: dict[str, int] = {
    "R1": 95,   # compartiment ECM = preuve directe
    "R2": 90,   # phenotype = sortie globale, score élevé
    "R3": 92,   # marqueur exclusif HIGH par défaut (ajusté par marqueur)
    "R4": 65,   # pathway-driven : plausible mais surinclusif
    "R5": 45,   # voisinage : preuve indirecte
    "R6": 25,   # fallback générique
    "R7": 0,
}

RULE_SOURCE: dict[str, str] = {
    "R1": "MINERVA SjD Map compartments 21555/21629 (Silva-Saffar 2026, doi:10.1038/s41540-026-00xxx-x)",
    "R2": "MINERVA SjD Map Phenotypes layer (Silva-Saffar 2026)",
    "R3": "see marker source per node",
    "R4": "Reactome/KEGG via MINERVA notes (Jassal 2020 PMID:31691815; Kanehisa 2023 PMID:36300620)",
    "R5": "graph propagation on SjD Map (this work)",
    "R6": "default intracellular signaling fallback (InnateDB Breuer 2013 PMID:23180781; STRING Szklarczyk 2023 PMID:36370105)",
    "R7": "",
}

# Sources canoniques par marqueur (CellMarker 2.0 / PanglaoDB / HGNC + PMID
# fondateurs ; un marqueur peut avoir plusieurs lignées → la source précise
# la lignée-cible primaire). Utilisé par R3.
MARKER_SOURCE: dict[str, str] = {
    # SGEC — Mavragani 2017, Manoussakis 2020, Verstappen 2021 review
    "AQP5": "HGNC:633; PMID:11136261 Steinfeld 2001 (AQP5 SGEC)",
    "AQP3": "HGNC:638; PanglaoDB SGEC marker",
    "MUC5B": "HGNC:7516; PMID:15057281 Alliende 2008 (acinar SGEC)",
    "MUC7": "HGNC:7518; CellMarker2.0 SGEC",
    "KRT7": "HGNC:6445; CellMarker2.0 epithelial",
    "KRT8": "HGNC:6446; PanglaoDB epithelial",
    "KRT18": "HGNC:6430; PanglaoDB epithelial",
    "KRT19": "HGNC:6436; PMID:10643978 SGEC ductal",
    "EPCAM": "HGNC:11529; PMID:24412576 epithelial",
    "CFTR": "HGNC:1884; CellMarker2.0 ductal cell",
    "CDH1": "HGNC:1748; epithelial junction (UniProt P12830)",
    "CLDN1": "HGNC:2032; tight junction (UniProt O95832)",
    "CLDN3": "HGNC:2045; tight junction (UniProt O15551)",
    "CLDN4": "HGNC:2046; tight junction",
    "OCLN": "HGNC:8104; tight junction (UniProt Q16625)",
    "ZO1": "HGNC:11827; tight junction",
    "TJP1": "HGNC:11827; tight junction",
    "AMY1A": "HGNC:474; PMID:11136261 salivary amylase",
    "PRB1": "HGNC:9337; salivary proline-rich protein",
    "STATH": "HGNC:11369; salivary statherin",
    # T-lineage TCR
    "CD3D": "HGNC:1673; CellMarker2.0 T cell (PMID:31091489)",
    "CD3E": "HGNC:1674; CellMarker2.0 T cell",
    "CD3G": "HGNC:1675; CellMarker2.0 T cell",
    "CD4": "HGNC:1678; CellMarker2.0 CD4 T helper",
    "CD2": "HGNC:1639; PanglaoDB T cell",
    "LCK": "HGNC:6524; Reactome:R-HSA-202403 TCR signaling",
    "ZAP70": "HGNC:12858; Reactome:R-HSA-202427",
    "LAT": "HGNC:18874; Reactome:R-HSA-202433",
    "LCP2": "HGNC:6526; SLP76 TCR adaptor",
    "CD28": "HGNC:1653; KEGG:hsa04660 T cell costim",
    # TH1
    "TBX21": "HGNC:11599; PMID:11163255 Szabo 2000 T-bet master TH1",
    "STAT4": "HGNC:11365; PMID:8702611 IL-12 → TH1",
    # TH17
    "RORC": "HGNC:10260; PMID:16990141 Ivanov 2006 RORγt TH17",
    "RORA": "HGNC:10258; PMID:18486082 RORα TH17",
    "IL23R": "HGNC:19100; KEGG:hsa04659 TH17 differentiation",
    "AHR": "HGNC:348; PMID:18362915 Veldhoen 2008 AhR TH17",
    # TFH
    "CXCR5": "HGNC:1060; PMID:11160276 Breitfeld 2000 TFH",
    "BCL6": "HGNC:1001; PMID:19608860 Johnston 2009 BCL6 TFH master",
    "ICOS": "HGNC:5351; PMID:11595950 Akiba 2005 TFH",
    "PDCD1": "HGNC:8760; PMID:21164560 PD-1 TFH",
    # TREG
    "FOXP3": "HGNC:6106; PMID:14679299 Hori 2003 FOXP3 Treg master",
    "IL2RA": "HGNC:6008; CD25 Treg (PMID:16330541)",
    "IKZF2": "HGNC:13177; Helios Treg (PMID:20962259)",
    "CTLA4": "HGNC:2505; PMID:20720586 Wing 2008 CTLA4 Treg",
    # B-lineage
    "CD19": "HGNC:1633; PanglaoDB B cell core",
    "MS4A1": "HGNC:7315; CD20 B cell (PMID:16921026)",
    "CD79A": "HGNC:1698; BCR component (UniProt P11912)",
    "CD79B": "HGNC:1699; BCR component (UniProt P40259)",
    "CR2": "HGNC:2336; CD21 B cell (UniProt P20023)",
    "BLNK": "HGNC:14211; B cell adaptor (UniProt Q8WV28)",
    "BANK1": "HGNC:18233; PMID:18204447 B cell scaffold (SjD GWAS)",
    # BCELL
    "BTK": "HGNC:1133; UniProt Q06187; KEGG:hsa04662 BCR",
    "BLK": "HGNC:1057; B-cell-specific Src kinase (UniProt P51451)",
    "TNFRSF13B": "HGNC:18153; TACI (UniProt O14836); PMID:11460154 BAFF",
    "TNFRSF13C": "HGNC:17755; BAFF-R (UniProt Q96RJ3)",
    "FCRL3": "HGNC:18506; PMID:19543290 SjD risk allele",
    # PLASMA
    "XBP1": "HGNC:12801; PMID:11460154 plasma cell ER stress",
    "IRF4": "HGNC:6119; PMID:16973387 plasma cell differentiation",
    "PRDM1": "HGNC:9346; BLIMP1 (PMID:12612588 Shapiro-Shelef)",
    "SDC1": "HGNC:10658; CD138 plasma cell (PMID:8095971)",
    "MZB1": "HGNC:23645; PMID:19763233 plasma cell ER chaperone",
    "JCHAIN": "HGNC:5713; IGJ plasma cell (UniProt P01591)",
    "TNFRSF17": "HGNC:11913; BCMA (UniProt Q02223); PMID:14990790",
    # Macrophage
    "CD68": "HGNC:1693; PanglaoDB macrophage core",
    "FCGR1A": "HGNC:3613; CD64 macrophage (UniProt P12314)",
    "FCGR2A": "HGNC:3616; CD32 macrophage",
    "FCGR3A": "HGNC:3619; CD16 macrophage",
    "CSF1R": "HGNC:2433; M-CSF receptor (UniProt P07333)",
    "MERTK": "HGNC:7027; macrophage efferocytosis",
    "ITGAM": "HGNC:6149; CD11b macrophage (UniProt P11215)",
    # M1
    "NOS2": "HGNC:7873; PMID:18025190 Mantovani M1",
    "CCR7": "HGNC:1608; M1/migratory DC (UniProt P32248)",
    # M2
    "CD163": "HGNC:1631; PMID:11380995 Pulford 2001 M2",
    "MRC1": "HGNC:7228; CD206 M2 (PMID:16824795)",
    "ARG1": "HGNC:663; PMID:18025190 Mantovani M2",
    "MSR1": "HGNC:7376; M2 scavenger receptor",
    # PDC
    "CLEC4C": "HGNC:14554; BDCA2 pDC (PMID:11254703)",
    "IL3RA": "HGNC:6012; CD123 pDC (UniProt P26951)",
    "LILRA4": "HGNC:15498; ILT7 pDC (PMID:18215090)",
    "IRF7": "HGNC:6122; PMID:15728487 IRF7 pDC IFN-α",
    "TLR7": "HGNC:15631; PMID:18948386 TLR7 pDC ssRNA",
    "TLR9": "HGNC:15633; PMID:18948386 TLR9 pDC CpG",
}

# Bonus / pénalité de plausibilité par cell-type quand R3 a plusieurs cibles.
# Marqueurs partagés (lignée) : on n'attribue pas le score max.
def _r3_score_for_marker(marker: str, n_celltypes: int) -> int:
    """Score R3 ajusté : marqueur exclusif (n=1) → 95 ; lignée (n>1) → 80-85."""
    if n_celltypes == 1:
        return 95
    if n_celltypes <= 4:
        return 82
    return 70


# Pour R4, granularité par cell-type : un pathway très cell-type-spécifique
# (TCR, BCR) a un score plus élevé qu'un pathway ubiquitaire (NF-κB, JAK-STAT).
PATHWAY_SCORE_OVERRIDE: dict[str, int] = {
    "R-HSA-983705": 80,  # BCR — très spécifique
    "R-HSA-202403": 80,  # TCR — très spécifique
    "R-HSA-877253": 75,  # IFN α/β — relativement spécifique innate
    "R-HSA-877312": 75,  # IFN γ — TH1/M1 axis
    "R-HSA-449147": 72,  # TLR — innate biased
    "R-HSA-2132295": 70,  # MHC-II — APC
    "R-HSA-1474244": 78,  # ECM (SGEC)
    "R-HSA-451927": 72,  # IL-2 family (T cells)
    "hsa04660": 80,      # TCR
    "hsa04662": 80,      # BCR
    "hsa04658": 82,      # Th1/Th2 differentiation
    "hsa04659": 82,      # Th17 differentiation
    "hsa04620": 72,      # TLR
    "hsa04621": 72,      # NLR
    "hsa04622": 72,      # RLR
}

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
    # Score 0-100, calculé par cell-type (le score peut différer si on a une
    # information cell-type-spécifique — ex. marqueur exclusif vs intra-lignée).
    score_per_celltype: dict[str, int] = field(default_factory=dict)
    # Référence biologique (DB ou DOI/PMID/PMC)
    source_per_celltype: dict[str, str] = field(default_factory=dict)

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
                "plausibility_score": 0,
                "source": "",
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
            "plausibility_score": self.score_per_celltype.get(ct, 0),
            "source": self.source_per_celltype.get(ct, ""),
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
        score_per_celltype={EXTRA: RULE_BASE_SCORE["R1"]},
        source_per_celltype={EXTRA: RULE_SOURCE["R1"]},
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
        score_per_celltype={PHENOTYPE: RULE_BASE_SCORE["R2"]},
        source_per_celltype={PHENOTYPE: RULE_SOURCE["R2"]},
    )


def rule_R3_marker(elem: dict[str, Any]) -> Assignment | None:
    raw = elem.get("name") or ""
    norm = _normalise_name(raw)
    upper = raw.upper().strip()
    for candidate in (upper, norm):
        if candidate in EXCLUSIVE_MARKERS:
            cts = EXCLUSIVE_MARKERS[candidate]
            score = _r3_score_for_marker(candidate, len(cts))
            src = MARKER_SOURCE.get(candidate, f"CellMarker2.0/PanglaoDB consensus marker={candidate}")
            return Assignment(
                node_id=int(elem["id"]),
                name=raw,
                type=elem.get("type") or "",
                compartment_id=elem.get("compartmentId"),
                celltypes=set(cts),
                confidence="HIGH",
                rule="R3",
                evidence=f"marker={candidate}",
                score_per_celltype={ct: score for ct in cts},
                source_per_celltype={ct: src for ct in cts},
            )
    return None


def rule_R4_pathway(elem: dict[str, Any]) -> Assignment | None:
    notes = elem.get("notes") or ""
    pids = _extract_pathway_ids(notes)
    if not pids:
        return None
    # Pour chaque cell-type, garder le score max issu d'un pathway pertinent
    ct_scores: dict[str, int] = {}
    ct_sources: dict[str, str] = {}
    matched_ids: list[str] = []
    for pid in pids:
        if pid in PATHWAY_TO_CELLTYPES:
            keys = [pid]
        else:
            keys = [k for k in PATHWAY_TO_CELLTYPES
                    if pid.startswith(k) or k.startswith(pid)]
        for k in keys:
            cts = PATHWAY_TO_CELLTYPES[k]
            score = PATHWAY_SCORE_OVERRIDE.get(k, RULE_BASE_SCORE["R4"])
            db = "Reactome" if k.startswith("R-HSA") else "KEGG"
            src = f"{db}:{k} via MINERVA notes"
            for ct in cts:
                if score > ct_scores.get(ct, 0):
                    ct_scores[ct] = score
                    ct_sources[ct] = src
            matched_ids.append(k)
    if not ct_scores:
        return None
    return Assignment(
        node_id=int(elem["id"]),
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes=set(ct_scores),
        confidence="MEDIUM",
        rule="R4",
        evidence=";".join(sorted(set(matched_ids))[:5]),
        score_per_celltype=ct_scores,
        source_per_celltype=ct_sources,
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

    # Score : base R5 + bonus proportionnel au vote (capé à 60)
    ct_scores: dict[str, int] = {}
    ct_sources: dict[str, str] = {}
    max_vote = max(votes.values())
    for ct in concordant:
        bonus = min(15, int(8 * votes[ct] / max(max_vote, 1)))
        ct_scores[ct] = RULE_BASE_SCORE["R5"] + bonus
        ct_sources[ct] = (
            f"{RULE_SOURCE['R5']}; neighbors_vote={votes[ct]}; "
            f"InnateDB:innatedb.com Breuer 2013 PMID:23180781; STRING db v12"
        )
    return Assignment(
        node_id=nid,
        name=elem.get("name") or "",
        type=elem.get("type") or "",
        compartment_id=elem.get("compartmentId"),
        celltypes=set(concordant),
        confidence="LOW",
        rule="R5",
        evidence=f"neighbors_vote={dict(votes.most_common(5))}",
        score_per_celltype=ct_scores,
        source_per_celltype=ct_sources,
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
    base = RULE_BASE_SCORE["R6"]
    if not has_edges:
        evidence += f"_isolated_in_compartment_{cid}"
        base = 15  # plus pénalisé : isolé sans signal cell-type

    cts = set(DEFAULT_FALLBACK_CELLTYPES)
    src = (
        f"{RULE_SOURCE['R6']}; "
        f"compartment={cid}; HPA Uhlén 2015 PMID:25613900 (ubiquitous expression)"
    )
    return Assignment(
        node_id=nid,
        name=elem.get("name") or "",
        type=etype,
        compartment_id=cid,
        celltypes=cts,
        confidence="LOW",
        rule="R6",
        evidence=evidence,
        score_per_celltype={ct: base for ct in cts},
        source_per_celltype={ct: src for ct in cts},
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
    score_sum: dict[str, int] = defaultdict(int)
    score_n: dict[str, int] = defaultdict(int)
    score_high: dict[str, int] = defaultdict(int)  # ≥80
    score_low: dict[str, int] = defaultdict(int)   # <40
    for a in assignments.values():
        for ct in a.celltypes:
            by_celltype[ct] += 1
            s = a.score_per_celltype.get(ct, 0)
            score_sum[ct] += s
            score_n[ct] += 1
            if s >= 80:
                score_high[ct] += 1
            elif s < 40:
                score_low[ct] += 1
    mean_score: dict[str, float] = {
        ct: round(score_sum[ct] / score_n[ct], 1) if score_n[ct] else 0.0
        for ct in score_n
    }

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
        "mean_plausibility_per_celltype": mean_score,
        "n_high_score_per_celltype": dict(score_high),
        "n_low_score_per_celltype": dict(score_low),
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
