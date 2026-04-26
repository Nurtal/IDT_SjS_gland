"""
intercellular.py — Table curée des edges intercellulaires (Phase 1.4).

Chaque entrée décrit une paire ligand-récepteur (LR) : qui sécrète le ligand
(`source_celltypes`), qui exprime le récepteur (`target_celltypes`), et le
mécanisme. Les flags `in_cellphonedb`, `in_omnipath`, `sjs_specific`
permettent de filtrer/auditer la provenance.

Sources :
    - CellPhoneDB v4 (Efremova 2020, PMID:32103204) — base LR curée
    - OmniPath (Türei 2016, PMID:27993971) — interactions inter-cellulaires
    - Littérature SjD-spécifique (Mavragani, Lavie, Manganelli, Mariette,
      Bombardieri, Salomonsson, Risselada — voir PMIDs par entrée)

Couverture obligatoire (Gate 1.4) :
    ✓ IFN-α (pDC → SGEC)
    ✓ BAFF (SGEC, M1, M2 → BCELL)
    ✓ CXCL13 (SGEC → BCELL, TFH)
    ✓ IL-21 (TFH → BCELL)

Convention CellDesigner (Zerrouk 2024) :
    - mechanism="secreted"     → ligand sécrété, transport extracellulaire
    - mechanism="contact"      → cell-cell direct (Heterodimer Complex Association)
    - mechanism="autocrine"    → boucle ligand+récepteur même cellule
"""

from __future__ import annotations

from dataclasses import dataclass


# Raccourcis cell-types
ALL_T = {"TH1", "TH17", "TFH", "TREG"}
ALL_B = {"BCELL", "PLASMA"}
ALL_MYE = {"M1", "M2", "PDC"}
APC = {"M1", "M2", "PDC", "BCELL", "SGEC"}


@dataclass(frozen=True)
class IntercellularEdge:
    ligand: str             # MINERVA name (ex. "TNFSF13B" pour BAFF)
    receptor: str           # MINERVA name (ex. "TNFRSF13C" pour BAFF-R)
    source_celltypes: frozenset[str]
    target_celltypes: frozenset[str]
    mechanism: str          # secreted | contact | autocrine
    in_cellphonedb: bool
    in_omnipath: bool
    sjs_specific: bool
    evidence: str           # PMID(s) + commentaire court


def _e(
    ligand: str, receptor: str,
    source: set[str], target: set[str],
    mechanism: str = "secreted",
    cpdb: bool = True, op: bool = True, sjs: bool = False,
    evidence: str = "",
) -> IntercellularEdge:
    return IntercellularEdge(
        ligand=ligand, receptor=receptor,
        source_celltypes=frozenset(source),
        target_celltypes=frozenset(target),
        mechanism=mechanism,
        in_cellphonedb=cpdb, in_omnipath=op, sjs_specific=sjs,
        evidence=evidence,
    )


# ---------------------------------------------------------------------------
# Table curée
# ---------------------------------------------------------------------------

INTERCELLULAR_EDGES: tuple[IntercellularEdge, ...] = (

    # =====================================================================
    # AXES OBLIGATOIRES SjD (Gate 1.4)
    # =====================================================================

    # IFN-α : pDC → SGEC (signature IFN type I, Mavragani 2017)
    _e("IFNA", "IFNAR1",  {"PDC"}, {"SGEC"}, sjs=True,
       evidence="PMID:28604219 Mavragani IFN-α SjD; PMID:15728487 IRF7 pDC; CellPhoneDB; OmniPath"),
    _e("IFNA", "IFNAR2",  {"PDC"}, {"SGEC"}, sjs=True,
       evidence="PMID:28604219; PMID:15728487"),
    _e("IFNA", "IFNAR1",  {"PDC"}, {"BCELL", "TH1", "M1", "M2"},
       evidence="PMID:15728487 pDC IFN-α systemic"),
    _e("IFNA", "IFNAR2",  {"PDC"}, {"BCELL", "TH1", "M1", "M2"},
       evidence="PMID:15728487"),

    # BAFF (TNFSF13B) : SGEC, M1, M2 → BCELL (Lavie 2004)
    _e("TNFSF13B", "TNFRSF13C", {"SGEC", "M1", "M2"}, {"BCELL"}, sjs=True,
       evidence="PMID:14722095 Lavie BAFF SGEC; PMID:11460154 Mackay BAFF/BAFFR; CellPhoneDB"),
    _e("TNFSF13B", "TNFRSF13B", {"SGEC", "M1", "M2"}, {"BCELL", "PLASMA"}, sjs=True,
       evidence="PMID:14722095; PMID:11460154 BAFF/TACI"),
    _e("TNFSF13B", "TNFRSF17", {"M1", "M2"}, {"PLASMA"}, sjs=True,
       evidence="PMID:11460154 BAFF/BCMA plasma cell survival"),

    # APRIL (TNFSF13) : SGEC, M1 → PLASMA, BCELL
    _e("TNFSF13", "TNFRSF13B", {"SGEC", "M1"}, {"BCELL", "PLASMA"}, sjs=True,
       evidence="PMID:11460154 APRIL/TACI; PMID:18613817 SjD APRIL"),
    _e("TNFSF13", "TNFRSF17", {"SGEC", "M1"}, {"PLASMA"}, sjs=True,
       evidence="PMID:11460154 APRIL/BCMA plasma survival"),

    # CXCL13 : SGEC → BCELL, TFH (TLS formation, Bombardieri 2007)
    _e("CXCL13", "CXCR5", {"SGEC"}, {"BCELL", "TFH"}, sjs=True,
       evidence="PMID:17389321 Bombardieri CXCL13 SjD TLS; PMID:11160276 CXCR5 Tfh/B"),

    # IL-21 : TFH → BCELL (Linterman 2010)
    _e("IL21", "IL21R", {"TFH"}, {"BCELL", "PLASMA"}, sjs=True,
       evidence="PMID:11574546 Tfh IL-21; PMID:20176269 Linterman SjD IL-21"),
    # IL-21 : TFH → TFH (autocrine reinforcement)
    _e("IL21", "IL21R", {"TFH"}, {"TFH"}, mechanism="autocrine",
       evidence="PMID:18599325 Tfh IL-21 autocrine"),

    # =====================================================================
    # AXE IFN-γ (TH1/PDC → SGEC, M1)
    # =====================================================================
    _e("IFNG", "IFNGR1", {"TH1", "PDC"}, {"SGEC", "M1", "M2", "BCELL"},
       evidence="PMID:11163255 Szabo Th1; CellPhoneDB; OmniPath"),
    _e("IFNG", "IFNGR2", {"TH1", "PDC"}, {"SGEC", "M1", "M2", "BCELL"},
       evidence="PMID:11163255"),

    # =====================================================================
    # AXE TNF (M1, TH1 → SGEC, autres)
    # =====================================================================
    _e("TNF", "TNFRSF1A", {"M1", "TH1"}, {"SGEC", "BCELL", "M1", "M2"},
       evidence="PMID:15530839 Mantovani M1; CellPhoneDB"),
    _e("TNF", "TNFRSF1B", {"M1", "TH1"}, ALL_T | {"BCELL", "M1", "M2", "PDC"},
       evidence="PMID:8521391 TNFR2 leukocyte"),

    # IFN-β : SGEC → SGEC + lymphoïdes (Mavragani type I IFN local)
    _e("IFNB1", "IFNAR1", {"SGEC"}, {"SGEC", "BCELL", "TH1"}, sjs=True,
       mechanism="autocrine",
       evidence="PMID:28604219 SjD IFN-β SGEC autocrine"),
    _e("IFNB1", "IFNAR2", {"SGEC"}, {"SGEC", "BCELL", "TH1"}, sjs=True,
       mechanism="autocrine",
       evidence="PMID:28604219"),

    # =====================================================================
    # AXE FAS / FasL — apoptose SGEC (Manganelli 2003)
    # =====================================================================
    _e("FASLG", "FAS", {"TH1"}, {"SGEC"}, mechanism="contact", sjs=True,
       evidence="PMID:12796328 Manganelli FasL SGEC apoptosis; PMID:7530740 Fas/FasL"),
    _e("FASLG", "FAS", {"M1"}, {"SGEC"}, mechanism="contact", sjs=True,
       evidence="PMID:12796328"),

    # =====================================================================
    # AXE LIGHT/LTα1β2 — destruction épithéliale (Pérol 2018)
    # =====================================================================
    _e("TNFSF14", "LTBR", {"TH1", "TH17"}, {"SGEC"}, sjs=True,
       evidence="PMID:9716579 LIGHT/LTβR; PMID:30237481 SjD LIGHT epithelial damage"),
    _e("LTA", "LTBR", {"TH1", "BCELL"}, {"SGEC"}, sjs=True,
       evidence="PMID:9716579 LTα/LTβR stromal"),
    _e("LTB", "LTBR", {"TH1", "BCELL"}, {"SGEC"}, sjs=True,
       evidence="PMID:9716579"),

    # =====================================================================
    # AXE IL-6 — pléiotrope (SGEC, M1 → multiple)
    # =====================================================================
    _e("IL6", "IL6R", {"SGEC", "M1", "TH17", "BCELL"},
       ALL_T | {"BCELL", "PLASMA", "M1", "M2", "SGEC"}, sjs=True,
       evidence="PMID:23597562 IL-6 SjD; CellPhoneDB; OmniPath"),
    _e("IL6", "IL6ST", {"SGEC", "M1", "TH17", "BCELL"},
       ALL_T | {"BCELL", "PLASMA", "M1", "M2", "SGEC"}, sjs=True,
       evidence="PMID:23597562 IL-6/gp130"),

    # =====================================================================
    # AXES IL-7, IL-15 — survie lymphocytaire (SGEC → T, B)
    # =====================================================================
    _e("IL7", "IL7R", {"SGEC"}, ALL_T | {"BCELL"}, sjs=True,
       evidence="PMID:7935811 IL-7R lymphoid; PMID:23728041 IL-7 SjD"),
    _e("IL15", "IL15R", {"SGEC", "M1", "M2", "PDC"},
       ALL_T | {"BCELL", "M1"},
       evidence="PMID:7973650 IL-15 myeloid/epithelial"),
    _e("IL15", "IL2RB", {"SGEC", "M1", "M2", "PDC"},
       ALL_T | {"M1"},
       evidence="PMID:7973650"),

    # =====================================================================
    # AXE IL-17 (TH17 → SGEC) — neutrophil recruitment
    # =====================================================================
    _e("IL17A", "IL17RA", {"TH17"}, {"SGEC", "M1"}, sjs=True,
       evidence="PMID:18172165 Th17 chemotaxis; PMID:18606613 IL-17 SjD"),

    # =====================================================================
    # AXE IL-12 / IL-23 (M1, PDC → T cells)
    # =====================================================================
    _e("IL12", "IL12RB1", {"M1", "PDC"}, {"TH1", "TH17"},
       evidence="PMID:9008747 IL-12 myeloid; PMID:8676086 IL-12Rβ1"),
    _e("IL12", "IL12RB2", {"M1", "PDC"}, {"TH1"},
       evidence="PMID:9151895 IL-12Rβ2 Th1"),
    _e("IL23A", "IL23R", {"M1", "PDC"}, {"TH17"},
       evidence="PMID:11114383 IL-23 myeloid → Th17"),

    # =====================================================================
    # AXE IL-2 (T → T)
    # =====================================================================
    _e("IL2", "IL2RA", {"TH1", "TH17", "TFH"}, ALL_T,
       evidence="PMID:6431080 IL-2 T cell"),
    _e("IL2", "IL2RB", {"TH1", "TH17", "TFH"}, ALL_T | {"PDC", "M1"},
       evidence="PMID:9590175 IL-2Rβ"),
    _e("IL2", "IL2RG", {"TH1", "TH17", "TFH"}, ALL_T | ALL_B | {"PDC"},
       evidence="γc"),

    # =====================================================================
    # AXE IL-10 (TREG, M2, BCELL → multiple — anti-inflammatoire)
    # =====================================================================
    _e("IL10", "IL10RA", {"TREG", "M2", "BCELL"},
       ALL_MYE | ALL_T | ALL_B | {"SGEC"},
       evidence="PMID:11244051 IL-10 Treg/M2/Breg"),

    # =====================================================================
    # AXE TGF-β (TREG, M2 → multiple)
    # =====================================================================
    _e("TGFB1", "TGFBR1", {"TREG", "M2", "SGEC"}, ALL_T | ALL_B | ALL_MYE | {"SGEC"},
       evidence="PMID:14679299 Hori FOXP3 Treg/TGF-β; CellPhoneDB"),

    # =====================================================================
    # AXES CHEMOKINES (recrutement)
    # =====================================================================
    # CXCR3 axis (Th1 hallmark)
    _e("CXCL9", "CXCR3", {"M1", "PDC", "SGEC"}, {"TH1", "M1", "PDC"}, sjs=True,
       evidence="PMID:9419197 CXCR3 Th1; PMID:28604219 SjD MIG"),
    _e("CXCL10", "CXCR3", {"M1", "PDC", "SGEC"}, {"TH1", "M1", "PDC"}, sjs=True,
       evidence="PMID:9419197; PMID:28604219 IP-10 SjD core"),
    _e("CXCL11", "CXCR3", {"M1", "SGEC"}, {"TH1", "M1", "PDC"}, sjs=True,
       evidence="PMID:9419197; PMID:28604219"),
    # CCR2 axis (monocyte recruitment)
    _e("CCL2", "CCR2", {"SGEC", "M1", "M2"}, {"M1", "M2"},
       evidence="PMID:9038108 MCP-1/CCR2 monocyte"),
    # CCR5 axis (Th1/Mφ)
    _e("CCL5", "CCR5", {"TH1", "M1"}, {"TH1", "M1", "M2"},
       evidence="PMID:9038108 RANTES/CCR5"),
    _e("CCL4", "CCR5", {"M1", "M2", "TH1"}, {"TH1", "M1", "M2"},
       evidence="PMID:9038108 MIP-1β/CCR5"),
    _e("CCL3", "CCR5", {"M1", "PDC"}, {"TH1", "M1", "M2"},
       evidence="PMID:9038108 MIP-1α/CCR5"),
    _e("CCL3", "CCR1", {"M1", "PDC"}, {"M1", "M2", "TH1"},
       evidence="PMID:9038108"),
    # CCR6 axis (Th17, B)
    _e("CCL20", "CCR6", {"SGEC", "M1"}, {"TH17", "BCELL"},
       evidence="PMID:18172165 CCL20/CCR6 Th17"),
    # CCR7 axis (LN homing)
    _e("CCL19", "CCR7", {"SGEC", "PDC"}, {"M1", "PDC", "TFH"},
       evidence="PMID:18258608 CCL19/CCR7 LN homing"),
    _e("CCL21", "CCR7", {"SGEC"}, {"M1", "PDC", "TFH"}, sjs=True,
       evidence="PMID:18258608 SjD ectopic GC; CellPhoneDB"),
    # CXCL8/IL-8 axis (neutrophil)
    _e("CXCL8", "CXCR1", {"M1", "SGEC"}, {"M1"},
       evidence="PMID:9038108 IL-8/CXCR1"),
    _e("CXCL8", "CXCR2", {"M1", "SGEC"}, {"M1"},
       evidence="PMID:9038108 IL-8/CXCR2"),
    # SDF-1/CXCR4 (homeostatic)
    _e("CXCL12", "CXCR4", {"SGEC"}, ALL_T | ALL_B | ALL_MYE | {"SGEC"},
       evidence="PMID:9038108 SDF-1/CXCR4"),
    # XCL1/2 → XCR1 (cDC1/pDC)
    _e("XCL1", "XCR1", {"TH1"}, {"PDC"},
       evidence="PMID:19204311 XCL/XCR1 pDC"),
    _e("XCL2", "XCR1", {"TH1"}, {"PDC"},
       evidence="PMID:19204311"),
    # CX3CL1/CX3CR1 (patrolling)
    _e("CX3CL1", "CX3CR1", {"SGEC"}, {"M1", "M2", "TH1"},
       evidence="PMID:9038108 fractalkine"),

    # =====================================================================
    # AXE CD40L / CD40 (T → B, APC)  — contact
    # =====================================================================
    _e("CD40LG", "CD40", {"TH1", "TH17", "TFH"}, {"BCELL", "M1", "M2", "PDC", "SGEC"},
       mechanism="contact",
       evidence="PMID:8551252 CD40L/CD40; CellPhoneDB"),

    # =====================================================================
    # AXE ICOS / ICOSL (T ↔ APC)
    # =====================================================================
    _e("ICOSLG", "ICOS", {"BCELL", "M1", "M2"}, {"TFH", "TH17"},
       mechanism="contact",
       evidence="PMID:11160276 ICOS Tfh; CellPhoneDB"),

    # =====================================================================
    # AXE PD-1 / PD-L1 (immune checkpoint)
    # =====================================================================
    _e("CD279", "CD274", {"TFH", "TH1"}, {"SGEC", "M1", "M2", "BCELL"},
       mechanism="contact", sjs=False,
       evidence="PMID:11015443 PD-1/PD-L1; CellPhoneDB"),

    # =====================================================================
    # AXE CD80/86 / CD28 (co-stim)  — contact APC ↔ T
    # =====================================================================
    _e("CD80", "CD28", APC, ALL_T, mechanism="contact",
       evidence="PMID:7510905 B7-1/CD28; CellPhoneDB"),
    _e("CD86", "CD28", APC, ALL_T, mechanism="contact",
       evidence="PMID:7510905 B7-2/CD28; CellPhoneDB"),

    # =====================================================================
    # MHC-II / TCR (APC → CD4 T cells)  — contact
    # =====================================================================
    _e("MHC1 proteins", "TCR", APC, ALL_T, mechanism="contact",
       evidence="MHC-II/TCR engagement (peptide presentation)"),

    # =====================================================================
    # AXE FLT3L → FLT3 (DC differentiation)
    # =====================================================================
    _e("FLT3LG", "FLT3", {"SGEC", "M1"}, {"PDC"},
       evidence="PMID:14523393 FLT3L pDC differentiation; CellPhoneDB"),

    # =====================================================================
    # AXE C3a, C5a → C3aR1, C5aR1 (complément, recrutement myéloïde)
    # =====================================================================
    _e("C3a", "C3AR1", {"SGEC"}, {"M1", "M2", "PDC"},
       evidence="PMID:8810309 C3a myeloid"),
    _e("C5a", "C5AR1", {"SGEC"}, {"M1", "M2", "PDC"},
       evidence="PMID:8810309 C5a myeloid"),

    # =====================================================================
    # AXE GM-CSF (CSF2) — survie myéloïde
    # =====================================================================
    _e("CSF2", "CSF2RA", {"TH1", "PDC"}, {"M1", "M2"},
       evidence="PMID:7530740 GM-CSF myeloid"),

    # =====================================================================
    # AXES CCR3 (M2/eosinophil) - moins prioritaire
    # =====================================================================
    _e("CCL11", "CCR3", {"SGEC", "M2"}, {"M2"},
       evidence="PMID:8676086 eotaxin"),
    _e("CCL13", "CCR2", {"M1", "M2"}, {"M1", "M2"},
       evidence="PMID:9038108 MCP-4"),

    # =====================================================================
    # AXES TLR ligands extracellulaires (PAMPs/DAMPs)
    # =====================================================================
    _e("LPS", "TLR4", {"SGEC"}, {"M1", "M2", "SGEC", "BCELL"},
       evidence="PMID:9237759 LPS/TLR4"),
    _e("ssRNA", "TLR7", {"SGEC"}, {"PDC", "BCELL"},
       evidence="PMID:12407176 ssRNA/TLR7"),
    _e("CpG_DNA", "TLR9", {"SGEC"}, {"PDC", "BCELL"},
       evidence="PMID:12407176 CpG/TLR9"),

    # =====================================================================
    # AXE CCL25/CCR9 (homing)
    # =====================================================================
    _e("CCL25", "CCR9", {"SGEC"}, {"TH1", "PLASMA"},
       evidence="PMID:11733579 CCL25/CCR9 gut homing"),

    # =====================================================================
    # AXE EDA/EDAR (épithélial)
    # =====================================================================
    _e("EDA", "EDAR", {"SGEC"}, {"SGEC"}, mechanism="autocrine",
       evidence="PMID:11030618 ectodysplasin epithelial"),

    # =====================================================================
    # AXE TNFSF10 (TRAIL) — apoptose
    # =====================================================================
    _e("TNFSF10", "TNFRSF10A", {"M1", "TH1", "PDC"}, {"SGEC", "BCELL"}, sjs=True,
       evidence="PMID:9311998 TRAIL; SjD epithelial apoptosis"),
    _e("TNFSF10", "TNFRSF10B", {"M1", "TH1", "PDC"}, {"SGEC", "BCELL"}, sjs=True,
       evidence="PMID:9311998"),

    # =====================================================================
    # IL-1β (M1, SGEC → multiple)
    # =====================================================================
    _e("IL1B", "IL1R1", {"M1", "SGEC"}, ALL_MYE | {"SGEC", "BCELL"} | ALL_T,
       evidence="PMID:8120391 IL-1R1 myeloid/epithelial"),

    # =====================================================================
    # IL-4 (TFH → BCELL — Th2-like, IgE)
    # =====================================================================
    _e("IL4", "IL4R", {"TFH"}, {"BCELL"},
       evidence="PMID:18599325 Tfh IL-4"),

    # =====================================================================
    # RANKL/RANK (T → BCELL/M)
    # =====================================================================
    _e("TNFSF11", "TNFRSF11A", ALL_T, {"M1", "M2", "BCELL"},
       evidence="PMID:9168116 RANKL/RANK"),

    # =====================================================================
    # TL1A/DR3 (M → T)
    # =====================================================================
    _e("TNFSF15", "TNFRSF25", {"M1", "M2"}, ALL_T,
       evidence="PMID:18500207 TL1A/DR3 T cell"),
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def edges_by_ligand() -> dict[str, list[IntercellularEdge]]:
    out: dict[str, list[IntercellularEdge]] = {}
    for e in INTERCELLULAR_EDGES:
        out.setdefault(e.ligand, []).append(e)
    return out


def mandatory_axes_covered(
    instantiated: list[dict[str, str]],
) -> dict[str, bool]:
    """Vérifie que les 4 axes obligatoires Gate 1.4 sont représentés."""
    pairs = {(e["ligand"], e["receptor"], e["source_celltype"], e["target_celltype"])
             for e in instantiated}

    def _hit(ligand: str, source: str, target: str) -> bool:
        return any(
            l == ligand and s == source and t == target
            for (l, _r, s, t) in pairs
        )

    return {
        "IFN-α (PDC→SGEC)": (
            _hit("IFNA", "PDC", "SGEC")
        ),
        "BAFF (SGEC/M1/M2→BCELL)": any(
            _hit("TNFSF13B", s, "BCELL") for s in ("SGEC", "M1", "M2")
        ),
        "CXCL13 (SGEC→BCELL/TFH)": any(
            _hit("CXCL13", "SGEC", t) for t in ("BCELL", "TFH")
        ),
        "IL-21 (TFH→BCELL)": (
            _hit("IL21", "TFH", "BCELL")
        ),
    }
