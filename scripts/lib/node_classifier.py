"""
node_classifier.py — Assignation des espèces de la SjD Map aux modules cellulaires.

Stratégie en 3 couches (appliquées dans l'ordre, première couche gagnante) :

    Couche 1 — Reactome pathway IDs (annotations objectives issues de la carte)
    Couche 2 — Gene symbol sets (connaissance biologique encodée manuellement)
    Couche 3 — Propagation de contexte sur le graphe de réactions (heuristique)

Modules cibles :
    SGEC     — Salivary Gland Epithelial Cells
    CD4      — CD4+ T cells (Th1/Th17)
    BCELL    — B cells / plasmablasts
    MACRO    — M1/M2 macrophages
    SHARED   — Espèces partagées délibérément entre modules (hubs signalétiques)
    EXTRACELLULAR — Ligands extracellulaires (compartiment Extracellular/Secreted)
    UNASSIGNED    — Non classifiées (à réviser manuellement)

Fonction publique principale :
    assign_modules(elements, reactions) → dict[int, str]

    Retourne {minerva_element_id → module_name}.

Ce fichier est conçu pour être édité itérativement :
    - Ajouter des entrées dans REACTOME_CELL_TYPE_MAP après exploration notebook
    - Compléter les *_MARKERS sets selon les résultats de l'audit
    - Ajouter des overrides manuels dans MANUAL_OVERRIDES
"""

from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Couche 1 — Reactome pathway IDs → module cellulaire
# ---------------------------------------------------------------------------
# Source : annotations Reactome dans le champ 'references' des éléments MINERVA.
# Les IDs sont extraits du lien (ex. "https://reactome.org/PathwayBrowser/#/R-HSA-909733")
# Format court attendu : "R-HSA-XXXXXX"

REACTOME_CELL_TYPE_MAP: dict[str, str] = {
    # ---------- SGEC (épithélium des glandes salivaires) ----------
    "R-HSA-909733":  "SGEC",  # Interferon alpha/beta signalling
    "R-HSA-877300":  "SGEC",  # Interferon alpha/beta signalling (alt)
    "R-HSA-936440":  "SGEC",  # Regulation of innate immune responses to SARS
    "R-HSA-168164":  "SGEC",  # Intrinsic pathway for apoptosis
    "R-HSA-109581":  "SGEC",  # Apoptosis
    "R-HSA-432047":  "SGEC",  # Aquaporin-mediated transport
    "R-HSA-5620971": "SGEC",  # Signalling by PTK6
    "R-HSA-5633007": "SGEC",  # Regulation of TP53 expression and degradation
    "R-HSA-5663213": "SGEC",  # RHO GTPases activate PAKs (epithelial)
    "R-HSA-2559586": "SGEC",  # DNA repair (epithelial cells)
    "R-HSA-5358747": "SGEC",  # Sensing of DNA double strand breaks
    "R-HSA-1169408": "SGEC",  # ISG15 antiviral mechanism
    "R-HSA-1169410": "SGEC",  # Antiviral mechanism by IFN-stimulated genes
    "R-HSA-168256":  "SGEC",  # Immune evasion by RSV (IFN pathway nodes)
    "R-HSA-6807505": "SGEC",  # RNA polymerase II transcription termination
    "R-HSA-9013148": "SGEC",  # CDC42 GTPase cycle (epithelial)

    # ---------- CD4+ T cells (Th1/Th17) ----------
    "R-HSA-202403":  "CD4",   # TCR signalling
    "R-HSA-202427":  "CD4",   # Downstream TCR signalling
    "R-HSA-202433":  "CD4",   # Generation of second messenger molecules
    "R-HSA-2454202": "CD4",   # Fc epsilon receptor (FCERI) — Th2 context
    "R-HSA-6804756": "CD4",   # Regulation of TP53 degradation (T cells)
    "R-HSA-5663220": "CD4",   # RHO GTPases activate PAKs (T cells)
    "R-HSA-389948":  "CD4",   # PD-1 signalling
    "R-HSA-202409":  "CD4",   # Costimulation by the CD28 family
    "R-HSA-6804757": "CD4",   # Regulation of TP53 expression (T cells)
    "R-HSA-5625740": "CD4",   # RHO GTPases activate KTN1 (T cells)
    "R-HSA-5218920": "CD4",   # VEGFA-VEGFR2 Pathway (Th17 angiogenesis)

    # ---------- B cells / plasmablastes ----------
    "R-HSA-983705":  "BCELL", # BCR signalling
    "R-HSA-5690714": "BCELL", # CD22 mediated BCR regulation
    "R-HSA-2029481": "BCELL", # FCGR activation (B cells)
    "R-HSA-5690716": "BCELL", # Activation of BIM and translocation to mitochondria
    "R-HSA-9013149": "BCELL", # CDC42 GTPase cycle (B cells)
    "R-HSA-2029482": "BCELL", # Regulation of actin dynamics for phagocytic cup
    "R-HSA-1257604": "BCELL", # PIP3 activates AKT signalling (BCR)
    "R-HSA-6806003": "BCELL", # Regulation of TP53 expression (B cells)
    "R-HSA-5625900": "BCELL", # RHO GTPases activate CIT (B cells)

    # ---------- Macrophages M1/M2 ----------
    "R-HSA-168249":  "MACRO", # Innate Immune System
    "R-HSA-168928":  "MACRO", # DDX58/IFIH1-mediated induction of interferon-alpha/beta
    "R-HSA-2029480": "MACRO", # FCGR-dependent phagocytosis
    "R-HSA-2029482": "MACRO", # Regulation of actin dynamics (phagocytic cup)
    "R-HSA-2871809": "MACRO", # FCERI-mediated Ca2+ mobilisation
    "R-HSA-5620612": "MACRO", # Signalling by FGFR4 (macrophage polarisation)
    "R-HSA-5625970": "MACRO", # RHO GTPases regulate CFTR (macrophages)
    "R-HSA-9013148": "MACRO", # CDC42 GTPase cycle (macrophages — partagé)
    "R-HSA-5663202": "MACRO", # RHO GTPases activate PAKs (macrophages)
    "R-HSA-168256":  "MACRO", # Innate immunity evasion
    "R-HSA-9660821": "MACRO", # SARS-CoV-2 innate immunity evasion
}

# ---------------------------------------------------------------------------
# Couche 2 — Gene symbol sets par module cellulaire
# ---------------------------------------------------------------------------

SGEC_MARKERS: frozenset[str] = frozenset({
    # Aquaporins / fonction sécrétoire
    "AQP5", "AQP3", "AQP1",
    # Mucines épithéliales
    "MUC1", "MUC5B", "MUC7",
    # BAFF / APRIL (produits par les SGEC dans la SjD)
    "TNFSF13", "TNFSF13B",
    # Récepteurs IFN de type I
    "IFNAR1", "IFNAR2",
    # MHC-I (présentation antigène épithéliale)
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
    # MHC-II (présentation aberrante dans SjD)
    "HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DMB",
    "CD74",  # chaîne invariante MHC-II
    # Voie JAK-STAT / IFN
    "STAT1", "STAT2", "JAK1", "TYK2",
    "IRF1", "IRF3", "IRF7", "IRF9",
    # ISGs (gènes stimulés par l'interféron — signature clé SjD)
    "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3",
    "IFI44", "IFI44L", "IFI6", "IFITM1", "IFITM3",
    "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
    # Apoptose / mort cellulaire épithéliale
    "CASP3", "CASP8", "CASP9", "CASP7",
    "BAX", "BAD", "BCL2", "BCL2L1", "BCL2L11",
    "CYCS", "APAF1",
    "FAS", "FASLG",
    # Voie NF-κB dans les épithélia
    "TNFRSF1A", "TNFRSF1B",  # récepteurs TNF
    # Sensing innée épithélial
    "TLR3", "TLR7", "TLR9",
    "MAVS", "STING1", "CGAS",
    "DDX58",  # RIG-I
    # Chimiokines produites par SGEC
    "CXCL10", "CXCL9", "CXCL11",  # CXCR3 ligands → recrutement T
    "CXCL13",  # recrutement B cells
    "CCL2", "CCL5",
    # Polarité épithéliale / jonctions
    "EPCAM", "CDH1",  # E-cadhérine
    # Stress du réticulum endoplasmique
    "XBP1",   # Partagé SGEC/BCELL — contexte dépendant
    "ATF6", "EIF2AK3",
})

CD4_MARKERS: frozenset[str] = frozenset({
    # Complexe TCR
    "CD3E", "CD3D", "CD3G", "CD247",  # CD3ζ
    "CD4",
    # Kinases proximal TCR
    "ZAP70", "LCK", "FYN",
    # Co-stimulation
    "CD28", "CD80", "CD86",
    "CD40LG",  # CD40L — exprimé sur les T activés
    "ICOS", "ICOSLG",
    # Checkpoints
    "CTLA4", "PDCD1",  # PD-1
    # Facteurs de transcription master Th1/Th17
    "TBX21",   # T-bet → Th1
    "RORC",    # RORγt → Th17
    "GATA3",   # Th2
    "FOXP3",   # Treg
    "BCL6",    # Tfh (partagé avec BCELL — contexte)
    # Cytokines T effectrices
    "IFNG",    # IFN-γ — signature Th1 majeure SjD
    "IL17A", "IL17F",  # Th17
    "IL21",    # Tfh / Th17
    "IL4", "IL5", "IL13",  # Th2
    "IL2", "IL2RA",   # croissance T
    "TGFB1",   # Treg / Th17 plasticity
    "IL10",    # Treg (partagé MACRO)
    "IL22",    # Th17 tissue
    # Voie PI3K/AKT dans T cells
    "PIK3CD", "AKT1",
    # Recrutement (récepteurs chimiokines)
    "CXCR3",   # ligand : CXCL10/9/11 produit par SGEC
    "CCR5", "CCR6", "CXCR5",
    # Signalisation en aval TCR
    "PLCG1", "PRKCA",
    "NFATC1", "NFATC2",
    "RELA",   # NF-κB (partagé — contexte T)
    # Activation / prolifération
    "MKI67",  # prolifération
    "CD69",   # activation précoce
    "CD44",
})

BCELL_MARKERS: frozenset[str] = frozenset({
    # Marqueurs de surface B
    "CD19", "CD20",   # MS4A1
    "CD22", "CD24",
    # Complexe BCR
    "CD79A", "CD79B",
    # BTK pathway
    "BTK", "BLNK",
    "PIK3CD", "PLCG2",
    "SYK",
    # Récepteurs BAFF/APRIL (produits par SGEC et MACRO)
    "TNFRSF13B",  # BAFF-R
    "TNFRSF13C",  # BCMA
    "TNFRSF17",   # TACI
    # Facteurs de transcription B / plasmablastes
    "PAX5",   # identité B cell
    "PRDM1",  # Blimp-1 → plasmablaste
    "IRF4",   # plasmablaste / MBC
    "BCL6",   # centre germinatif
    "AICDA",  # AID — commutation isotypique
    "XBP1",   # plasmablaste (stress RE)
    "MZB1",   # marginal zone B
    # Anticorps / immunoglobulines
    "IGHG1", "IGHG2", "IGHG3", "IGHG4",
    "IGHA1", "IGHA2",
    "IGHM",
    "IGLC1", "IGLC2",
    # Voie PI3K/AKT B cells
    "AKT1", "MTOR",
    # Lymphomagenèse (nœuds critiques SjD)
    "MYC",
    "BCL2",   # anti-apoptose B
    "CCND1",  # cycline D1
    # Signalisation CD40
    "CD40",   # récepteur → signal de survie B
    # Récepteurs chimiokines / recrutement
    "CXCR5",  # recrutement dans les follicules
    "CCR6",
    # Co-stimulation
    "CD86",   # exprimé sur B activés
    "PDCD1LG2",  # PD-L2
})

MACRO_MARKERS: frozenset[str] = frozenset({
    # Marqueurs macrophages
    "CD68", "CD163", "CD64",   # FCGR1A
    "MRC1",   # CD206 — mannose receptor, M2
    "MARCO",  # scavenger receptor
    "CD36",   # phagocytose/lipides
    # TLR signalling (innate sensing)
    "TLR1", "TLR2", "TLR4", "TLR5", "TLR6",
    "TLR7", "TLR8", "TLR9",
    "MYD88", "TIRAP",
    "TICAM1",  # TRIF
    # Voie NF-κB innée
    "IRAK1", "IRAK4",
    "TRAF6",
    "CHUK",   # IKKα
    "IKBKB",  # IKKβ
    "IKBKG",  # NEMO
    # Inflammasome
    "NLRP3",
    "PYCARD",  # ASC
    "CASP1",
    "IL1B",   # output inflammasome
    "IL18",
    "IL33",
    # Cytokines M1
    "TNF",
    "IL6",
    "IL12A", "IL12B",  # IL-12 → différenciation Th1
    "IL23A",           # IL-23 → Th17
    "CXCL8",           # IL-8
    # Cytokines M2 / résolution
    "IL10",
    "TGFB1",
    "VEGFA",
    # Enzymes M1/M2
    "NOS2",   # iNOS → M1
    "ARG1",   # arginase → M2
    # Récepteurs Fc (phagocytose complexes immuns)
    "FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A",
    # Facteurs de transcription M1/M2
    "IRF5",   # M1
    "PPARG",  # M2
    "STAT6",  # M2 (IL-4/IL-13)
    "STAT3",  # M2 (IL-10/IL-6)
    # Phagocytose / lysosome
    "LAMP1",
    "CTSD",   # cathepsine D
    # BAFF/APRIL production (macrophages aussi producteurs dans SjD)
    "TNFSF13B",  # BAFF — partagé avec SGEC
    "TNFSF13",   # APRIL — partagé avec SGEC
})

# ---------------------------------------------------------------------------
# Espèces extracellulaires / ligands sécrétés
# ---------------------------------------------------------------------------
# Ces noms apparaissent dans le compartiment Extracellular ou Secreted.
# Ils seront assignés à EXTRACELLULAR (ni intra-cellulaire, ni module spécifique).

EXTRACELLULAR_MARKERS: frozenset[str] = frozenset({
    # Cytokines extracellulaires génériques
    "IFNA1", "IFNA2", "IFNB1",
    "IFNG",   # extracellulaire (produit CD4, reçu SGEC)
    "IL17A", "IL17F",
    "IL21",
    "TNF",    # extracellulaire
    "IL6",
    "IL12",
    "IL10",
    "TGFB1",
    "IL4",
    "BAFF",   # = TNFSF13B extracellulaire
    "APRIL",  # = TNFSF13 extracellulaire
    "CXCL10", "CXCL9", "CXCL11", "CXCL13",
    "CCL2", "CCL5",
})

# ---------------------------------------------------------------------------
# Hubs partagés délibérément dupliqués entre plusieurs modules
# ---------------------------------------------------------------------------
# Ces nœuds sont upstream de ≥2 types cellulaires et doivent être présents
# dans chaque module avec un suffixe (_SGEC, _CD4, etc.)

SHARED_HUBS: frozenset[str] = frozenset({
    "STAT1",    # IFN → SGEC + MACRO
    "STAT3",    # IL-6 → CD4 + MACRO
    "RELA",     # NF-κB → SGEC + CD4 + MACRO
    "NFKB1",    # NF-κB → idem
    "NFKB2",
    "IRF9",     # ISGF3 → SGEC + MACRO
    "MAPK1",    # ERK2 → ubiquitaire
    "MAPK3",    # ERK1
    "MAPK8",    # JNK1
    "MAPK14",   # p38α
    "AKT1",     # PI3K → BCELL + CD4 + MACRO
    "PIK3CA",   # PI3K → idem
    "TP53",     # apoptose → ubiquitaire
    "BCL2",     # anti-apoptose → SGEC + BCELL
    "TGFB1",    # Treg/fibrosis → CD4 + MACRO + SGEC
    "IL10",     # résolution → CD4 + MACRO
    "TNF",      # M1 → MACRO + SGEC
    "CXCL10",   # SGEC produit, CD4 reçoit
    "IFNG",     # CD4 produit, SGEC reçoit
    "TNFSF13B", # BAFF → SGEC + MACRO produisent, BCELL reçoit
    "TNFSF13",  # APRIL → idem
    "BCL6",     # TFH/GC → CD4 + BCELL
    "XBP1",     # UPR → SGEC + BCELL
    "VEGFA",    # angiogenèse → MACRO + SGEC
})

# ---------------------------------------------------------------------------
# Overrides manuels (priorité maximale — couche 0)
# ---------------------------------------------------------------------------
# Format : {element_name_lower → module}
# À compléter itérativement après l'exploration notebook.

MANUAL_OVERRIDES: dict[str, str] = {
    # Nœuds phénotypiques → SHARED (ils apparaissent dans tous les modules)
    "inflammation":                   "SHARED",
    "apoptosis":                      "SHARED",
    "b_cell_activation/survival":     "SHARED",
    "mhc_class_1_activation":         "SHARED",
    "mhc_class_2_activation":         "SHARED",
    "t_cell_activation/differentiation": "SHARED",
    "lymphoid_organ_development":     "SHARED",
    "regulated_necrosis":             "SHARED",
    "phagocytosis":                   "SHARED",
    "chemotaxis/infiltration":        "SHARED",
    "fibrosis":                       "SHARED",
    "angiogenesis":                   "SHARED",
    "matrix_degradation":             "SHARED",
    "cell_proliferation/survival":    "SHARED",
}

# ---------------------------------------------------------------------------
# Fonctions utilitaires internes
# ---------------------------------------------------------------------------

def _normalize_name(name: str) -> str:
    """Normalise un nom pour la comparaison : minuscules, espaces → underscore."""
    return re.sub(r"[\s\-/]+", "_", (name or "").strip().lower())


def _extract_reactome_ids(elem: dict) -> list[str]:
    """
    Extrait les IDs Reactome depuis les références d'un élément MINERVA.

    Recherche dans :
        - elem['references'][*]['link'] → URL Reactome
        - elem['notes'] → texte libre avec mentions R-HSA-XXXXXX
    """
    ids: list[str] = []
    pattern = re.compile(r"R-HSA-\d+")

    for ref in elem.get("references", []) or []:
        ref_type = (ref.get("type") or "").upper()
        link = ref.get("link") or ref.get("resource") or ""
        if "REACTOME" in ref_type or "reactome.org" in link.lower():
            found = pattern.findall(link)
            ids.extend(found)

    notes = elem.get("notes") or ""
    ids.extend(pattern.findall(notes))

    return list(set(ids))


def _classify_by_reactome(elem: dict) -> str | None:
    """Couche 1 : classification par Reactome pathway IDs."""
    reactome_ids = _extract_reactome_ids(elem)
    votes: Counter[str] = Counter()
    for rid in reactome_ids:
        module = REACTOME_CELL_TYPE_MAP.get(rid)
        if module:
            votes[module] += 1
    if not votes:
        return None
    top_module, top_count = votes.most_common(1)[0]
    # Conflit = plusieurs modules avec le même score → SHARED
    if len(votes) > 1:
        second_count = votes.most_common(2)[1][1]
        if top_count == second_count:
            return "SHARED"
    return top_module


def _classify_by_symbol(elem: dict) -> str | None:
    """
    Couche 2 : classification par gene symbol.

    Cherche dans elem['name'] et dans les références de type HGNC_SYMBOL.
    """
    candidates: set[str] = set()

    name = (elem.get("name") or "").strip().upper()
    if name:
        candidates.add(name)

    for ref in elem.get("references", []) or []:
        if (ref.get("type") or "").upper() == "HGNC_SYMBOL":
            link = ref.get("link") or ""
            symbol = link.rstrip("/").split("/")[-1].upper()
            if symbol:
                candidates.add(symbol)

    matches: set[str] = set()

    for sym in candidates:
        if sym in SGEC_MARKERS:
            matches.add("SGEC")
        if sym in CD4_MARKERS:
            matches.add("CD4")
        if sym in BCELL_MARKERS:
            matches.add("BCELL")
        if sym in MACRO_MARKERS:
            matches.add("MACRO")
        if sym in EXTRACELLULAR_MARKERS:
            matches.add("EXTRACELLULAR")

    if not matches:
        return None
    if len(matches) == 1:
        return next(iter(matches))
    # Plusieurs modules → SHARED si c'est un hub connu, sinon UNASSIGNED
    if any(sym in SHARED_HUBS for sym in candidates):
        return "SHARED"
    return "UNASSIGNED"


def _classify_by_compartment(elem: dict) -> str | None:
    """
    Couche 3 (heuristique légère) : compartiment MINERVA comme indice.

    Les espèces dans Extracellular (21555) ou Secreted (21629) sans assignation
    sont placées dans EXTRACELLULAR.
    Les phénotypes (21540) sont SHARED.
    """
    cid = elem.get("compartmentId")
    if cid in (21555, 21629):
        return "EXTRACELLULAR"
    if cid == 21540:
        return "SHARED"
    return None


# ---------------------------------------------------------------------------
# Propagation de contexte sur le graphe (couche 3 complémentaire)
# ---------------------------------------------------------------------------

def _propagate_context(
    assignments: dict[int, str],
    reactions: list[dict],
    elements: list[dict],
    max_iterations: int = 3,
) -> dict[int, str]:
    """
    Propage les assignations connues aux nœuds UNASSIGNED adjacents dans le graphe.

    Pour chaque nœud UNASSIGNED :
        - Collecte les modules de tous ses voisins directs (réactants/produits/modifiers)
        - Si tous les voisins assignés appartiennent au même module → assigne ce module
        - Si conflit → reste UNASSIGNED (sera révisé manuellement)

    Args:
        assignments : Dict courant {element_id → module} (modifié in-place)
        reactions   : Liste des réactions MINERVA
        elements    : Liste des éléments MINERVA
        max_iterations : Nombre de passes de propagation

    Returns:
        Dict assignments mis à jour.
    """
    # Construction du graphe de voisinage
    neighbors: dict[int, set[int]] = defaultdict(set)
    for rxn in reactions:
        all_ids: list[int] = []
        for group_key in ("reactants", "products", "modifiers"):
            for entry in rxn.get(group_key) or []:
                elem_dict = entry.get("element") or {}
                eid = elem_dict.get("id")
                if eid is not None:
                    all_ids.append(int(eid))
        for i, eid_a in enumerate(all_ids):
            for eid_b in all_ids[i + 1:]:
                neighbors[eid_a].add(eid_b)
                neighbors[eid_b].add(eid_a)

    changed = True
    iteration = 0
    while changed and iteration < max_iterations:
        changed = False
        iteration += 1
        for eid, module in list(assignments.items()):
            if module != "UNASSIGNED":
                continue
            neighbor_modules = {
                assignments.get(nid)
                for nid in neighbors.get(eid, set())
                if assignments.get(nid) not in (None, "UNASSIGNED", "SHARED", "EXTRACELLULAR")
            }
            if len(neighbor_modules) == 1:
                new_module = next(iter(neighbor_modules))
                assignments[eid] = new_module
                changed = True
                logger.debug(
                    "Propagation iter%d : élément %d → %s",
                    iteration, eid, new_module,
                )

    return assignments


# ---------------------------------------------------------------------------
# Fonction publique principale
# ---------------------------------------------------------------------------

def assign_modules(
    elements: list[dict[str, Any]],
    reactions: list[dict[str, Any]],
    propagate: bool = True,
) -> dict[int, str]:
    """
    Assigne chaque espèce MINERVA à un module cellulaire.

    Stratégie en 4 couches (priorité décroissante) :
        0. Overrides manuels (MANUAL_OVERRIDES)
        1. Reactome pathway IDs
        2. Gene symbol sets
        3. Compartiment MINERVA (Extracellular/Phenotypes)
        + propagation de contexte sur le graphe (optionnel)

    Args:
        elements   : Liste des éléments MINERVA (résultat de get_all_elements)
        reactions  : Liste des réactions MINERVA (résultat de get_all_reactions)
        propagate  : Si True, propage les assignations aux voisins UNASSIGNED

    Returns:
        {minerva_element_id (int) → module_name (str)}
        module_name ∈ {"SGEC", "CD4", "BCELL", "MACRO",
                       "SHARED", "EXTRACELLULAR", "UNASSIGNED"}
    """
    assignments: dict[int, str] = {}

    for elem in elements:
        eid = elem.get("id")
        if eid is None:
            continue
        eid = int(eid)

        name_norm = _normalize_name(elem.get("name") or "")

        # Couche 0 — override manuel
        if name_norm in MANUAL_OVERRIDES:
            assignments[eid] = MANUAL_OVERRIDES[name_norm]
            continue

        # Couche 0b — nœud Phenotype → SHARED automatiquement
        if elem.get("type") == "Phenotype":
            assignments[eid] = "SHARED"
            continue

        # Couche 1 — Reactome
        module = _classify_by_reactome(elem)
        if module:
            assignments[eid] = module
            continue

        # Couche 2 — Gene symbols
        module = _classify_by_symbol(elem)
        if module:
            assignments[eid] = module
            continue

        # Couche 3 — Compartiment
        module = _classify_by_compartment(elem)
        if module:
            assignments[eid] = module
            continue

        # Non classifié
        assignments[eid] = "UNASSIGNED"

    # Propagation de contexte (couche 4)
    if propagate:
        assignments = _propagate_context(assignments, reactions, elements)

    return assignments


# ---------------------------------------------------------------------------
# Statistiques et rapport
# ---------------------------------------------------------------------------

def classification_report(
    assignments: dict[int, str],
    elements: list[dict],
) -> dict[str, Any]:
    """
    Retourne un rapport de classification pour QC.

    Returns:
        Dict avec :
            counts       : {module → count}
            pct_assigned : % d'éléments assignés (hors UNASSIGNED)
            unassigned   : liste des éléments non classifiés [{id, name, type}]
            shared       : liste des éléments partagés
    """
    id_to_elem = {int(e["id"]): e for e in elements if e.get("id") is not None}

    counts = Counter(assignments.values())
    total = len(assignments)
    n_unassigned = counts.get("UNASSIGNED", 0)
    pct_assigned = 100.0 * (total - n_unassigned) / total if total else 0.0

    unassigned_list = [
        {
            "id": eid,
            "name": id_to_elem.get(eid, {}).get("name", "?"),
            "type": id_to_elem.get(eid, {}).get("type", "?"),
            "compartmentId": id_to_elem.get(eid, {}).get("compartmentId"),
        }
        for eid, module in assignments.items()
        if module == "UNASSIGNED"
    ]

    shared_list = [
        {
            "id": eid,
            "name": id_to_elem.get(eid, {}).get("name", "?"),
        }
        for eid, module in assignments.items()
        if module == "SHARED"
    ]

    return {
        "counts": dict(counts.most_common()),
        "total": total,
        "pct_assigned": round(pct_assigned, 1),
        "unassigned": unassigned_list,
        "shared": shared_list,
    }


def get_module_elements(
    module: str,
    assignments: dict[int, str],
    elements: list[dict],
) -> list[dict]:
    """
    Retourne la liste des éléments MINERVA appartenant à un module donné.

    Args:
        module      : "SGEC", "CD4", "BCELL", "MACRO", "SHARED", "EXTRACELLULAR"
        assignments : Résultat de assign_modules()
        elements    : Liste complète des éléments

    Returns:
        Liste filtrée d'éléments.
    """
    module_ids = {eid for eid, m in assignments.items() if m == module}
    return [e for e in elements if e.get("id") is not None and int(e["id"]) in module_ids]


def check_required_markers(
    assignments: dict[int, str],
    elements: list[dict],
) -> dict[str, list[str]]:
    """
    Vérifie que les marqueurs essentiels de chaque module sont présents.

    Returns:
        {module → liste des marqueurs manquants}
        (dict vide = tout OK)
    """
    REQUIRED: dict[str, set[str]] = {
        "SGEC": {
            "AQP5", "TNFSF13B", "STAT1", "IRF7", "ISG15",
            "HLA-A", "HLA-DRA", "CASP3", "TLR3", "CXCL10",
        },
        "CD4": {
            "CD3E", "ZAP70", "IFNG", "IL17A", "TBX21",
            "RORC", "CD28", "CD40LG", "FOXP3", "CXCR3",
        },
        "BCELL": {
            "CD19", "CD79A", "BTK", "TNFRSF13B", "BCL6",
            "PRDM1", "AICDA", "CD40", "IRF4",
        },
        "MACRO": {
            "TLR4", "MYD88", "NLRP3", "IL1B", "TNF",
            "IL10", "CD68", "ARG1", "NOS2", "FCGR3A",
        },
    }

    # Index {name_upper → module}
    name_to_module: dict[str, str] = {}
    for e in elements:
        eid = e.get("id")
        if eid is None:
            continue
        module = assignments.get(int(eid), "UNASSIGNED")
        name = (e.get("name") or "").upper()
        if name:
            # Un nom peut exister dans plusieurs modules (SHARED) → on garde
            if name not in name_to_module or module != "UNASSIGNED":
                name_to_module[name] = module

    missing: dict[str, list[str]] = {}
    for module, required_set in REQUIRED.items():
        missing_markers = []
        for marker in required_set:
            assigned_module = name_to_module.get(marker.upper())
            if assigned_module is None:
                missing_markers.append(f"{marker} (absent de la carte)")
            elif assigned_module not in (module, "SHARED", "UNASSIGNED", "EXTRACELLULAR"):
                missing_markers.append(f"{marker} (assigné à {assigned_module})")
        if missing_markers:
            missing[module] = missing_markers

    return missing
