"""
refinement.py — Raffinement expert des assignations R6 (fallback générique).

Surcharche par-nom les assignations R6 trop larges en s'appuyant sur :
    • PanglaoDB (Franzén 2019 PMID:30951143)
    • CellMarker 2.0 (Hu 2023 PMID:36300620)
    • Human Protein Atlas (HPA, Uhlén 2015 PMID:25613900)
    • InnateDB (Breuer 2013 PMID:23180781)
    • Reactome curated pathways (Jassal 2020 PMID:31691815)
    • ImmGen tissue expression (Heng 2008 PMID:18800157)
    • SjD-specific literature (Mavragani 2017, Manoussakis 2020, Verstappen 2021 review)

Sortie : règle `R6c` (curated), score 60–85 selon spécificité du knowledge.

Le matching se fait sur le `node_name` MINERVA en majuscules (alphanum-only).
Les complexes ("STAT1/STAT2/IRF9", "TGF/TGFBR1") sont décomposés et chaque
sous-composant est checké ; on prend l'intersection des cell-types curés.
"""

from __future__ import annotations

import logging
import re
from typing import Any

logger = logging.getLogger(__name__)

# Raccourcis lignées
ALL = {"SGEC", "TH1", "TH17", "TFH", "TREG", "BCELL", "PLASMA", "M1", "M2", "PDC"}
ALL_T = {"TH1", "TH17", "TFH", "TREG"}
ALL_B = {"BCELL", "PLASMA"}
ALL_LYMPH = ALL_T | ALL_B
ALL_MYE = {"M1", "M2", "PDC"}
ALL_INNATE = ALL_MYE | {"SGEC"}
APC = {"M1", "M2", "PDC", "BCELL", "SGEC"}  # antigen-presenting cells (MHC-II)


# ---------------------------------------------------------------------------
# Table curée : name → (celltypes, score, source)
# ---------------------------------------------------------------------------
#
# Source format : "DB1; DB2; key PMIDs"
# Score :
#   85 : cell-type très restreint, consensus fort
#   80 : restriction modérée, validée par 2+ DBs
#   75 : axe biologique défini (ex. ISG → TOUS)
#   70 : fallback informé (HPA ubiquitaire mais pertinence variable)

CURATED: dict[str, tuple[set[str], int, str]] = {

    # =====================================================================
    # IFN signaling axis (signature centrale SjD — Mavragani 2017 PMID:28604219)
    # =====================================================================
    # Récepteurs IFN type I : ubiquitaires, expression confirmée HPA
    "IFNAR1": (ALL, 75, "HPA Uhlén 2015 PMID:25613900; Reactome:R-HSA-909733; ubiquitous"),
    "IFNAR2": (ALL, 75, "HPA Uhlén 2015; Reactome:R-HSA-909733"),
    "IFNGR1": (ALL, 75, "HPA Uhlén 2015; Reactome:R-HSA-877312"),
    "IFNGR2": (ALL, 75, "HPA Uhlén 2015; Reactome:R-HSA-877312"),
    # Cytokines IFN — produites par sous-ensembles
    "IFNA":   ({"PDC", "SGEC"}, 82, "PMID:15728487 IRF7 pDC; Mavragani 2017 PMID:28604219 SGEC"),
    "IFNB1":  ({"SGEC", "PDC", "M1", "M2"}, 80, "Reactome:R-HSA-909733; PMID:28604219; SGEC IFN-β"),
    "IFNG":   ({"TH1", "PDC", "M1"}, 82, "PMID:11163255 Szabo TBX21; ImmGen Th1 signature"),
    # Effecteurs JAK-STAT (downstream IFN/IL-6)
    "STAT1": (ALL, 80, "Reactome:R-HSA-877300; KEGG:hsa04630; ubiquitous JAK-STAT"),
    "STAT2": (ALL, 78, "Reactome:R-HSA-909733; ubiquitous"),
    "STAT3": (ALL, 78, "KEGG:hsa04630; ubiquitous IL-6/IL-10"),
    "STAT4 HOMODIMER": ({"TH1", "PDC", "M1"}, 82, "PMID:8702611 Th1; HPA T cell"),
    "STAT5 HOMODIMER": (ALL_T | {"BCELL"}, 78, "KEGG:hsa04630; PMID:9590175 STAT5 lymphoid"),
    "STAT5A B": (ALL_T | {"BCELL"}, 78, "KEGG:hsa04630"),
    "JAK3": (ALL_LYMPH | {"PDC"}, 82, "PMID:7481803 JAK3 γc; lymphoid-restricted"),
    "TYK2": (ALL, 75, "Reactome:R-HSA-909733; ubiquitous"),
    "IRF3": (ALL, 75, "Reactome:R-HSA-3134963; ubiquitous innate"),
    "IRF9": (ALL, 75, "Reactome:R-HSA-877300; STAT1/2/IRF9 complex"),

    # ISGs — globalement induits dans toutes cellules sous IFN, pertinents partout
    "ADAR":   (ALL, 75, "PMID:30401873 ISG; Reactome:R-HSA-9013973"),
    "ADAR1 HOMODIMER": (ALL, 75, "PMID:30401873 ADAR1 ISG"),
    "EIF2AK2": (ALL, 75, "PKR antiviral, Reactome:R-HSA-1834941"),
    "EIF2AK2 HOMODIMER": (ALL, 75, "PKR antiviral"),
    "ISG15":  (ALL, 75, "Reactome:R-HSA-1169408 ISGylation"),
    "ISG20":  (ALL, 75, "ISG, antiviral exonuclease PMID:23394386"),
    "HERC5":  (ALL, 75, "ISG E3 ligase ISGylation PMID:16275642"),
    "IFI27":  (ALL, 75, "ISG; Mavragani 2017 PMID:28604219 SjD signature"),
    "IFI35":  (ALL, 75, "ISG; PMID:28604219"),
    "IFI6":   (ALL, 75, "ISG; PMID:28604219 SjD signature"),
    "IFIT1":  (ALL, 75, "ISG; Reactome:R-HSA-936964"),
    "IFIT2":  (ALL, 75, "ISG; Reactome:R-HSA-936964"),
    "IFIT3":  (ALL, 75, "ISG; PMID:28604219"),
    "IFIT5":  (ALL, 75, "ISG; PMID:28604219"),
    "IFITM1": (ALL, 75, "ISG; PMID:23354158"),
    "IFITM2": (ALL, 75, "ISG; PMID:23354158"),
    "IFITM3": (ALL, 75, "ISG; PMID:23354158"),
    "MX1":    (ALL, 75, "ISG GTPase antiviral; PMID:28604219"),
    "MX2":    (ALL, 75, "ISG GTPase"),
    "GBP2":   (ALL, 75, "ISG; PMID:28604219"),
    "BST2":   (ALL, 75, "Tetherin ISG; PMID:18342598"),
    "OAS1":   (ALL, 75, "ISG OAS family; Reactome:R-HSA-1834941"),
    "OAS2":   (ALL, 75, "ISG OAS family"),
    "OAS3":   (ALL, 75, "ISG OAS family"),
    "XAF1":   (ALL, 70, "ISG; pro-apoptotic"),
    "UBA7":   (ALL, 70, "UBE1L ISGylation E1; PMID:16275642"),
    "TARBP2": (ALL, 70, "PKR partner"),
    "PSME2":  (APC, 75, "Immunoproteasome; PMID:21134086 IFN-γ induced"),
    "NLRC5":  (ALL, 72, "MHC-I transactivator; PMID:21084438"),

    # =====================================================================
    # TLR/NLR/RLR — innate sensing (InnateDB curated)
    # =====================================================================
    "TLR1":  ({"M1", "M2", "SGEC", "PDC"}, 82, "InnateDB; PMID:11607032 TLR1 surface"),
    "TLR2":  ({"M1", "M2", "SGEC"}, 82, "InnateDB; PMID:9237759 TLR2 myeloid+epithelial"),
    "TLR3":  ({"SGEC", "M1", "M2", "PDC"}, 82, "InnateDB; PMID:14976262 TLR3 dsRNA epithelial+DC"),
    "TLR4":  ({"M1", "M2", "SGEC", "BCELL"}, 82, "InnateDB; PMID:9237759 LPS receptor"),
    "TLR5":  ({"SGEC", "M1", "M2"}, 80, "InnateDB; PMID:11323673 flagellin epithelial"),
    "TLR8":  ({"M1", "M2", "PDC"}, 80, "InnateDB; PMID:12407176 ssRNA endosomal"),
    "TLR10": ({"BCELL", "M1", "M2", "PDC"}, 75, "InnateDB; PMID:14751764 B-cell biased"),
    "LY96":  ({"M1", "M2", "SGEC"}, 80, "MD-2/TLR4 partner; PMID:10359581"),
    "CD14":  ({"M1", "M2"}, 85, "PanglaoDB monocyte/macrophage core; PMID:1698311"),
    "MYD88": (ALL, 75, "Universal TLR/IL-1R adaptor; Reactome:R-HSA-166058"),
    "TIRAP": (ALL_INNATE | ALL_B, 75, "TLR2/4 adaptor; PMID:11544529"),
    "TICAM1": (ALL_INNATE, 78, "TRIF; TLR3/4 adaptor; PMID:14530373"),
    "TICAM2": (ALL_INNATE, 78, "TRAM; TLR4 adaptor; PMID:14530373"),
    "IRAK1 2": (ALL, 72, "Reactome:R-HSA-9020702 IL-1/TLR signal"),
    "IRAK4":   (ALL, 72, "Reactome:R-HSA-9020702; PMID:12055258"),
    "IRAK4 IRAK1 2 IRAKM": (ALL, 72, "Myddosome; Reactome:R-HSA-975144"),
    "IRAKM":   ({"M1", "M2"}, 78, "PMID:12150927 IRAK-M myeloid"),
    "TBK1":   (ALL, 75, "Reactome:R-HSA-936964; ubiquitous innate"),
    "IKBKE":  (ALL, 75, "IKKε/IRF3 axis; Reactome:R-HSA-936964"),
    "IKBKE TBK1": (ALL, 75, "IKKε/TBK1 complex"),
    "MAVS":   (ALL_INNATE | ALL_B, 78, "RLR adaptor; PMID:16125763"),
    "RIGI":   (ALL_INNATE | ALL_B, 78, "DDX58 RLR; PMID:15208624"),
    "IFIH1":  (ALL_INNATE | ALL_B, 78, "MDA5 RLR; PMID:15208624"),
    # Co-receptors
    "CD36":   ({"M1", "M2", "SGEC"}, 78, "PanglaoDB myeloid; scavenger receptor"),
    "TLR1 TLR2 CD14":  ({"M1", "M2"}, 80, "TLR1/2 heterodimer myeloid"),
    "TLR2 TLR6 CD36":  ({"M1", "M2"}, 80, "TLR2/6/CD36 myeloid"),
    "TLR4 LY96":  ({"M1", "M2", "SGEC"}, 80, "TLR4/MD-2 LPS"),
    "TLR5 FLAGELLIN": ({"SGEC", "M1", "M2"}, 78, "TLR5 epithelial PMID:11323673"),
    "TLR8 SSRNA": ({"M1", "M2", "PDC"}, 78, "TLR8 endosomal PMID:12407176"),
    "TLR10L TLR10": ({"BCELL", "M1", "M2"}, 70, "TLR10 B/myeloid"),

    # =====================================================================
    # TNF receptor family / death signaling
    # =====================================================================
    "TNF":    ({"M1", "TH1", "PDC"}, 82, "Mantovani 2004 PMID:15530839 M1; ImmGen Th1"),
    "TNFRSF1A": (ALL, 75, "TNFR1 ubiquitous; HPA Uhlén 2015"),
    "TNFRSF1B": ({"M1", "M2", "PDC"} | ALL_LYMPH, 75, "TNFR2 leukocyte-biased; PMID:8521391"),
    "TNF TNFR1": (ALL, 75, "TNFR1 complex"),
    "TNF TNFR2": (ALL_LYMPH | ALL_MYE, 75, "TNFR2 leukocyte"),
    "TNFRSF10A": (ALL, 72, "TRAIL-R1; HPA"),
    "TNFRSF10A B": (ALL, 72, "TRAIL-R1/2"),
    "TRAIL TRAILR": (ALL, 72, "TRAIL/DR4-5 PMID:9311998"),
    "TNFRSF11A": ({"M1", "M2", "BCELL"}, 78, "RANK; PMID:9168116 osteoclast/myeloid"),
    "RANKL RANK": ({"M1", "M2", "BCELL"} | ALL_T, 75, "RANK/RANKL"),
    "TNFRSF25": ({"TH1", "TH17", "TFH", "TREG"}, 78, "DR3; PMID:18500207 T cell"),
    "TNFSF15 TNFRSF25": (ALL_T, 78, "TL1A/DR3 T-cell axis PMID:18500207"),
    "TNFSF13B": ({"M1", "M2", "SGEC", "PDC"}, 82, "BAFF source: myeloid+SGEC; PMID:11460154 Mackay; Lavie 2004 SGEC"),
    "FAS":    (ALL, 75, "Fas ubiquitous; Reactome:R-HSA-75157"),
    "FASL FASR": ({"TH1", "M1"} | {"SGEC"}, 78, "FasL T cell+M1; SGEC apoptosis Manganelli 2003 PMID:12796328"),
    # Death machinery
    "FADD": (ALL, 72, "Death adaptor; ubiquitous"),
    "TRADD": (ALL, 72, "TNFR1 adaptor"),
    "RIPK1": (ALL, 72, "Reactome:R-HSA-2562578"),
    "RIPK1 TRAF2 TRAF5": (ALL, 72, "TNFR1 complex I"),
    "RIPK3": (ALL, 70, "Necroptosis"),
    "CFLAR": (ALL, 70, "cFLIP anti-apoptotic"),
    "DAXX":  (ALL, 70, "Death-associated"),
    "CASP3": (ALL, 75, "Executioner caspase"),
    "CASP7": (ALL, 75, "Executioner caspase"),
    "CASP8": (ALL, 75, "Initiator caspase extrinsic"),
    "BCL2":  (ALL, 75, "Anti-apoptotic; KEGG:hsa04210"),
    "BCL2A1": ({"M1", "M2", "PDC"} | ALL_LYMPH, 75, "Lymphoid/myeloid HPA"),
    "BCL2L1": (ALL, 75, "Bcl-xL anti-apoptotic"),
    "BAD":   (ALL, 72, "Pro-apoptotic"),
    "BIRC2": (ALL, 70, "cIAP1"),
    "XIAP":  (ALL, 70, "X-IAP"),
    "TNFAIP3": (ALL, 75, "A20 deubiquitinase NFkB"),
    "TRAF3": (ALL, 72, "TNFR/TLR adaptor"),
    "TRAF5": (ALL, 72, "TNFR adaptor"),
    "TRAF6": (ALL, 72, "TLR/IL-1R adaptor"),
    "TRAF2 TRAF3": (ALL, 72, "TRAF complex"),

    # =====================================================================
    # NF-κB axis
    # =====================================================================
    "NFKB1": (ALL, 78, "p50; ubiquitous KEGG:hsa04064"),
    "NFKB2": (ALL, 78, "p52; non-canonical NF-kB"),
    "NFKB2 RELB": (ALL, 78, "Non-canonical NFkB; PMID:11283238"),
    "RELA":  (ALL, 78, "p65 ubiquitous"),
    "RELB":  (ALL, 78, "Non-canonical NFkB"),
    "NFKBIA": (ALL, 75, "IκBα ubiquitous"),
    "CHUK":  (ALL, 75, "IKKα; KEGG:hsa04064"),
    "CHUK CHUK": (ALL, 75, "IKKα homodimer"),
    "IKBKB": (ALL, 75, "IKKβ canonical NF-kB"),
    "IKBKG": (ALL, 75, "NEMO"),
    "IKBIP": (ALL, 70, "IKK interactor"),
    "MAP3K7": (ALL, 75, "TAK1 TLR/IL-1R; PMID:10882101"),
    "MAP3K14": (ALL, 75, "NIK non-canonical NFkB"),
    "TAB1":  (ALL, 70, "TAK1 partner"),
    "TAB2":  (ALL, 70, "TAK1 partner"),
    "TAB3":  (ALL, 70, "TAK1 partner"),
    "TAB1 TAB2 TAB3": (ALL, 70, "TAK1/TAB complex"),
    "CARD11": (ALL_LYMPH, 82, "CBM signalosome; PMID:14684827 lymphoid"),
    "BCL10":  (ALL_LYMPH | {"M1", "M2"}, 80, "CBM signalosome PMID:14684827"),
    "MALT1":  (ALL_LYMPH | {"M1", "M2"}, 80, "CBM signalosome"),
    "CARD11 BCL10 MALT1": (ALL_LYMPH, 82, "CBM complex lymphoid"),

    # =====================================================================
    # MAPK common
    # =====================================================================
    "MAPK3":  (ALL, 75, "ERK1; ubiquitous"),
    "MAPK8 9 10": (ALL, 75, "JNK1/2/3"),
    "MAPK11 12 13 14": (ALL, 75, "p38 MAPK"),
    "MAP2K1": (ALL, 75, "MEK1 RAS-RAF-MEK-ERK"),
    "MAP2K3": (ALL, 75, "MKK3 p38"),
    "MAP2K4": (ALL, 75, "MKK4 JNK"),
    "MAP2K6": (ALL, 75, "MKK6 p38"),
    "MAP2K7": (ALL, 75, "MKK7 JNK"),
    "MAP3K1": (ALL, 75, "MEKK1"),
    "MAP3K5": (ALL, 75, "ASK1"),
    "DUSP5": (ALL, 70, "MAPK phosphatase"),
    "DUSP6": (ALL, 70, "MAPK phosphatase"),
    "JUN":   (ALL, 72, "AP-1"),
    "FOS":   (ALL, 72, "AP-1"),
    "AP1":   (ALL, 72, "AP-1 transcription factor"),
    "ETS1":  (ALL_LYMPH | ALL_MYE, 75, "Lymphoid/myeloid TF; PMID:9252188"),
    "MAX":   (ALL, 70, "Myc partner"),

    # =====================================================================
    # T-cell receptor / co-stim
    # =====================================================================
    "TCR":   (ALL_T, 88, "T-cell receptor T-lineage only"),
    "MHC1 B2M TCR": (ALL_T, 80, "MHC-I/TCR engagement"),
    "LCK ZAP70": (ALL_T, 88, "TCR proximal; PMID:7585978"),
    "CD8":   ({"TH1"}, 70, "CD8 T cells (not modeled separately) - keep on TH1 cytotoxic"),
    "CD27 CD70": (ALL_T | ALL_B | {"PDC"}, 75, "CD27/CD70 axis lymphoid PMID:9596702"),
    "CD70 HOMOTRIMER": (ALL_T | {"BCELL"}, 75, "CD70 activated lymphocyte"),
    "CD80 86": (APC, 82, "B7-1/B7-2 co-stim APC; PMID:7510905"),
    # IL-2 axis (T cell)
    "IL2":   ({"TH1", "TH17", "TFH"}, 80, "PMID:6431080 activated T cells"),
    "IL2R":  (ALL_T, 80, "IL-2 receptor"),
    "IL2RB": (ALL_T | {"PDC", "M1"}, 78, "CD122; PMID:9590175"),
    "IL2RG": (ALL_LYMPH | {"PDC"}, 80, "γc subunit lymphoid"),
    # IL-12 (Th1)
    "IL12":   ({"M1", "PDC"}, 82, "PMID:9008747 IL-12 myeloid/DC"),
    "IL12R":  ({"TH1", "PDC"}, 82, "IL-12R Th1+pDC"),
    "IL12RB1": ({"TH1", "TH17", "PDC"}, 80, "IL-12Rβ1 PMID:8676086"),
    "IL12RB2": ({"TH1"}, 85, "IL-12Rβ2 Th1-specific PMID:9151895"),
    # IL-15 (NK/T memory)
    "IL15":    ({"M1", "M2", "PDC", "SGEC"}, 78, "IL-15 myeloid/epithelial; PMID:7973650"),
    "IL15R":   (ALL_T | {"BCELL", "M1"}, 75, "IL-15R"),
    "IL15RA":  (ALL_T | {"M1", "M2", "PDC"}, 75, "IL-15Rα"),
    # IL-4/5 (Th2 — non modélisé)
    "IL4":    ({"TFH"}, 78, "Tfh IL-4 PMID:18599325 also Th2"),
    "IL5":    (set(), 0, "Th2 not modeled — DROP"),
    # IL-6 family
    "IL6":    ({"M1", "SGEC", "TH17", "BCELL"}, 80, "PMID:23597562 IL-6 source; SjD SGEC"),
    "IL6R":   (ALL, 75, "IL-6R; KEGG:hsa04630"),
    "IL6RA":  (ALL, 75, "IL-6Rα"),
    "IL6ST":  (ALL, 75, "gp130 ubiquitous"),
    # IL-7 (lymphoid development)
    "IL7R":   (ALL_T | {"BCELL"}, 80, "IL-7R lymphoid PMID:7935811"),
    "IL7RA":  (ALL_T | {"BCELL"}, 80, "IL-7Rα"),
    # IL-10
    "IL10":   ({"M2", "TREG", "BCELL"}, 82, "PMID:11244051 IL-10 M2/Treg/Breg"),
    # IL-21 (Tfh primary)
    "IL21":   ({"TFH", "TH17"}, 85, "PMID:11574546 Tfh IL-21"),
    "IL21R":  ({"BCELL", "PLASMA", "TFH", "TH17"}, 80, "IL-21R lymphoid"),
    # IL-23
    "IL23A":  ({"M1", "PDC"}, 82, "PMID:11114383 IL-23 myeloid"),
    # IL-1 family
    "IL1R1":  (ALL_INNATE | {"BCELL"}, 78, "IL-1R1; PMID:8120391"),
    "IL1RA":  (ALL, 70, "IL-1Ra antagonist"),

    # CSF
    "CSF2":   ({"M1", "PDC"}, 78, "GM-CSF myeloid PMID:7530740"),
    # FLT3 (DC differentiation)
    "FLT3":   ({"PDC"}, 88, "FLT3 pDC differentiation; PMID:14523393"),
    "FLT3LG FLT3": ({"PDC"}, 85, "FLT3L/FLT3 DC PMID:14523393"),

    # =====================================================================
    # Chemokine receptors & ligands (cell-type-restricted by receptor)
    # =====================================================================
    "CCR1":   ({"M1", "M2", "TH1"}, 78, "PMID:9038108 myeloid/Th1"),
    "CCR2":   ({"M1", "M2"}, 82, "PanglaoDB monocyte; PMID:9038108"),
    "CCR3":   ({"M2"}, 75, "Eosinophil/Th2; PMID:8676086 — keep on M2 only"),
    "CCR5":   ({"TH1", "M1", "M2"}, 80, "Th1/Mφ PMID:9038108"),
    "CCR6":   ({"TH17", "BCELL"}, 85, "PMID:18172165 Th17 chemotaxis"),
    "CCR7":   ({"M1", "PDC"} | {"TFH"}, 80, "CCR7 lymph node homing"),
    "CCR9":   ({"TH1", "PLASMA"}, 75, "Gut homing PMID:11733579"),
    "CCRL2":  ({"M1", "M2", "PDC"}, 72, "Atypical chemokine receptor"),
    "CXCR1":  ({"M1"}, 80, "IL-8R PMID:9038108"),
    "CXCR2":  ({"M1"}, 80, "IL-8R; neutrophil/M1"),
    "CXCR3":  ({"TH1", "M1", "PDC"}, 85, "PMID:9419197 Th1 hallmark"),
    "CXCR4":  (ALL, 72, "Ubiquitous SDF-1R"),
    "CXCR5":  ({"TFH", "BCELL"}, 88, "PMID:11160276 Tfh+B-cell"),
    "CXCR6":  ({"TH1", "TH17"}, 78, "PMID:11460154 effector T"),
    "CX3CR1": ({"M1", "M2", "TH1"}, 78, "PMID:9038108 patrolling Mφ"),
    "FPR2":   ({"M1", "M2"}, 78, "Lipoxin receptor myeloid"),
    "XCR1":   ({"PDC"}, 80, "PMID:19204311 cDC1/pDC"),
    # Ligand-receptor complexes
    "CCL11 CCR3": ({"M2"}, 75, "Eotaxin/CCR3"),
    "CCL16 CCR1": ({"M1", "M2"}, 72, "Myeloid"),
    "CCL19 CCRL2": (ALL_MYE, 70, "Atypical receptor"),
    "CCL2 CCR2":  ({"M1", "M2"}, 80, "Monocyte"),
    "CCL20 CCR6": ({"TH17", "BCELL"}, 80, "Th17/B"),
    "CCL25 CCR9": ({"TH1", "PLASMA"}, 75, "Gut homing"),
    "CCL5 CCR5":  ({"TH1", "M1", "M2"}, 78, "RANTES"),
    "CX3CL1 CX3CR1": ({"M1", "M2", "TH1"}, 75, "Fractalkine"),
    "CXCL1 CXCR2": ({"M1"}, 78, "Neutrophil"),
    "CXCL12 CXCR4": (ALL, 70, "SDF-1/CXCR4"),
    "CXCL16 CXCR6": ({"TH1", "TH17"}, 75, "Effector T"),
    "CXCL5 CXCR2": ({"M1"}, 78, "Neutrophil"),
    "CXCL6 CXCR1": ({"M1"}, 78, "Neutrophil"),
    "CXCL6 CXCR2": ({"M1"}, 78, "Neutrophil"),
    "CXCL9 10 11 CXCR3": ({"TH1", "M1", "PDC"}, 82, "Th1 hallmark PMID:9419197"),
    "XCL1 XCL2 XCR1": ({"PDC"}, 80, "PMID:19204311"),
    "LXA4 CCL23 FPR2": ({"M1", "M2"}, 72, "Lipoxin myeloid"),
    # Free ligands (endogenous expression — not extracellular compartment)
    "CCL2":  ({"SGEC", "M1", "M2"}, 78, "MCP-1 SGEC+myeloid PMID:9038108"),
    "CCL3":  ({"M1", "PDC"}, 78, "MIP-1α myeloid"),
    "CCL4":  ({"M1", "M2", "TH1"}, 78, "MIP-1β"),
    "CCL5":  ({"TH1", "M1"}, 78, "RANTES"),
    "CCL11": ({"SGEC", "M2"}, 75, "Eotaxin epithelial"),
    "CCL13": ({"M1", "M2"}, 72, "MCP-4 myeloid"),
    "CCL19": ({"SGEC", "PDC"}, 75, "Stromal+DC"),
    "CCL21": ({"SGEC"}, 78, "PMID:18258608 SjD ectopic GC SGEC source"),
    "CXCL1": ({"M1", "SGEC"}, 75, "GROα"),
    "CXCL2": ({"M1", "SGEC"}, 75, "GROβ"),
    "CXCL3": ({"M1"}, 72, "GROγ"),
    "CXCL5": ({"SGEC", "M1"}, 72, "ENA-78 epithelial"),
    "CXCL8": ({"M1", "SGEC"}, 78, "IL-8 source"),
    "CXCL9": ({"M1", "PDC", "SGEC"}, 80, "MIG IFN-γ-induced PMID:28604219 SjD"),
    "CXCL10": ({"M1", "PDC", "SGEC"}, 82, "IP-10 IFN-γ-induced; SjD core PMID:28604219"),
    "CXCL11": ({"M1", "SGEC"}, 78, "I-TAC IFN-induced"),
    "CXCL12": ({"SGEC"}, 75, "SDF-1 stromal"),

    # =====================================================================
    # Complement / Fc receptors
    # =====================================================================
    "C3AR1": ({"M1", "M2", "PDC"}, 78, "C3a receptor myeloid PMID:8810309"),
    "C3A CA3R1": ({"M1", "M2"}, 75, "C3a/C3aR myeloid"),
    "C5AR1": ({"M1", "M2", "PDC"}, 78, "C5a receptor myeloid"),
    "C5A C5AR1": ({"M1", "M2"}, 75, "C5a/C5aR myeloid"),
    "FCER1A": ({"PDC", "M1"}, 75, "FcεRI α — pDC; PMID:11160276"),
    "FCER1G": ({"M1", "M2", "PDC"}, 78, "FcεRI γ-chain shared myeloid"),
    "MS4A2":  ({"M1"}, 70, "FcεRI β — mast cell primary; minor M1"),
    "IGG":    ({"PLASMA"}, 88, "Plasma cell IgG product"),

    # =====================================================================
    # MHC / antigen presentation
    # =====================================================================
    "HLA-A": (ALL, 78, "MHC-I ubiquitous nucleated"),
    "HLA-B": (ALL, 78, "MHC-I"),
    "HLA-C": (ALL, 78, "MHC-I"),
    "HLA-G": ({"SGEC", "M2"}, 75, "Tolerogenic MHC-Ib PMID:11181016"),
    "HLA-H": (ALL, 70, "MHC-I pseudo"),
    # Adhesion
    "ICAM1": (ALL_INNATE | ALL_LYMPH, 78, "Adhesion induced by IFN-γ; PMID:2563381"),
    "VCAM1": ({"SGEC", "M1", "M2"}, 75, "Endothelial/epithelial adhesion"),
    "VCAM ITGA4 ITGB1": (ALL_LYMPH, 75, "VLA-4/VCAM lymphocyte trafficking"),
    "ITGA4": (ALL_LYMPH | ALL_MYE, 78, "VLA-4 leukocyte"),
    "ITGB1": (ALL, 72, "β1 integrin ubiquitous"),
    "ITGA4 ITGB1": (ALL_LYMPH | ALL_MYE, 78, "VLA-4"),

    # =====================================================================
    # Lymphocyte / B-cell-specific
    # =====================================================================
    "CD22":  ({"BCELL"}, 90, "PanglaoDB B cell core; CellMarker2.0"),
    "CD40":  (APC, 82, "B+APC; PMID:8551252"),
    "CD72":  ({"BCELL"}, 88, "PanglaoDB B cell"),
    "EDA":   ({"SGEC"}, 78, "Ectodysplasin epithelial PMID:11030618"),
    "EDAR":  ({"SGEC"}, 78, "EDA receptor epithelial"),
    "EDA EDAR": ({"SGEC"}, 78, "Ectodysplasin axis epithelial"),
    "EDARADD": ({"SGEC"}, 78, "EDA adaptor"),
    "LTBR":  ({"SGEC", "BCELL"}, 75, "LTβR stromal/epithelial PMID:9716579"),
    "LTB-TNFSF14 LTBR": ({"SGEC"}, 75, "LTα1β2/LIGHT to LTβR"),

    # =====================================================================
    # TGF-β axis (TREG, M2 emphasis)
    # =====================================================================
    "TGF TGFBR1": (ALL, 70, "TGF-β/TGFBR1"),
    "TGFBR1": (ALL, 75, "Ubiquitous"),
    "SMAD2": (ALL, 72, "TGF-β downstream"),
    "SMAD3": (ALL, 72, "TGF-β downstream"),
    "SMAD4": (ALL, 72, "TGF-β/BMP downstream"),
    "SMAD2 SMAD3": (ALL, 72, "R-Smads"),
    "SMAD2 SMAD3 SMAD4": (ALL, 72, "Smad complex"),

    # =====================================================================
    # Tumor suppressors / cell cycle / mTOR / housekeeping
    # =====================================================================
    "TP53":   (ALL, 70, "Ubiquitous tumor suppressor"),
    "MDM2":   (ALL, 70, "p53 regulator"),
    "CDKN1A": (ALL, 70, "p21 ubiquitous"),
    "GADD45B": (ALL, 70, "Stress response"),
    "MTOR":   (ALL, 72, "mTOR ubiquitous"),
    "MTORC1": (ALL, 72, "mTOR complex 1"),
    "MLST8":  (ALL, 70, "mTOR partner"),
    "RPTOR":  (ALL, 70, "mTORC1 partner"),
    "RPLP0":  (ALL, 70, "Ribosomal P0"),
    "LMNB1":  (ALL, 70, "Lamin B1 housekeeping"),
    "HNRNPA2B1": (ALL, 70, "RNA binding housekeeping"),
    "KPNB1":  (ALL, 70, "Importin β"),
    "VIM":    ({"SGEC", "M1", "M2", "PDC"}, 75, "Mesenchymal/myeloid HPA"),
    "ARF1":   (ALL, 70, "Trafficking ubiquitous"),
    "GORAB":  (ALL, 70, "Golgi"),
    "GRB2":   (ALL, 72, "Adaptor ubiquitous"),
    "SHC1":   (ALL, 72, "Adaptor ubiquitous"),
    "SOS2":   (ALL, 72, "RAS GEF"),
    "RAC1":   (ALL, 72, "Rho GTPase"),
    "RRAS":   (ALL, 70, "Ras family"),
    "RAF":    (ALL, 72, "MAPK kinase"),
    "SRC":    (ALL, 72, "Tyrosine kinase ubiquitous"),
    "PTK2":   (ALL, 70, "FAK ubiquitous"),
    "PTK2B":  (ALL_LYMPH | ALL_MYE, 75, "Pyk2 hematopoietic"),
    "PI3K":   (ALL, 72, "PI3K ubiquitous"),
    "PI(3,4,5)P3": (ALL, 70, "Lipid signal"),
    "PLCB":   (ALL, 70, "PLC-β"),
    "PPP2R1A": (ALL, 70, "PP2A scaffold"),
    "CALM1":  (ALL, 70, "Calmodulin"),
    "ITPR3":  (ALL, 70, "IP3R"),
    "CA2":    (ALL, 65, "Calcium signal"),
    "GLUT1":  (ALL, 70, "SLC2A1 ubiquitous"),
    "KDR":    (set(), 0, "VEGFR-2 endothelial — not modeled DROP"),
    "RRM2":   (ALL, 65, "Ribonucleotide reductase"),
    "MMP9":   ({"M1", "M2", "SGEC"}, 75, "MMP-9 myeloid+epithelial"),
    "PTGS2":  ({"M1", "SGEC"}, 75, "COX-2 inflammatory"),
    "MUC1":   ({"SGEC"}, 88, "PMID:9237759 epithelial mucin"),
    "CMTM6":  (ALL, 70, "PD-L1 stabilizer"),
    "CMKLR1": ({"M1", "M2", "PDC"}, 75, "ChemR23 myeloid"),
    "GNAI":   (ALL, 70, "Gαi ubiquitous"),
    "AP1":    (ALL, 70, "AP-1 TF"),
    # Phosphatases
    "PTPN2":  (ALL, 75, "TC-PTP ubiquitous PMID:21646515"),
    "PTPN6":  (ALL_LYMPH | ALL_MYE, 80, "SHP-1 hematopoietic PMID:7479798"),
    "PTPRC":  (ALL_LYMPH | ALL_MYE, 88, "CD45 pan-leukocyte; CellMarker2.0"),
    # Ubiquitin/SUMO
    "TANK":   (ALL, 70, "TBK1 partner"),
    "SOCS1":  (ALL, 75, "JAK-STAT inhibitor PMID:9590175"),
    "SOCS2":  (ALL, 72, "JAK-STAT inhibitor"),
    "SOCS3":  (ALL, 75, "JAK-STAT inhibitor"),
    "CISH":   (ALL, 72, "CIS JAK-STAT inhibitor"),
    "CMTM6":  (ALL, 70, "PD-L1 stabilizer PMID:28813417"),
    "HDAC3":  (ALL, 70, "Histone deacetylase"),
    "MARCHF8": (APC, 75, "MHC-II ubiquitin ligase PMID:18794341"),
    "LAMP3":   ({"PDC", "SGEC"}, 75, "DC-LAMP PMID:11160276; SGEC-DC like Mavragani 2017"),

    # =====================================================================
    # Lineage TFs
    # =====================================================================
    "TOX":    ({"TFH", "TH17"}, 78, "PMID:31375812 T cell exhaustion/Tfh"),
    "SNAI1":  ({"SGEC"}, 75, "EMT epithelial PMID:21326285"),

    # =====================================================================
    # SGEC-specific innate
    # =====================================================================
    "SIGLEC1": ({"M1", "M2", "PDC"}, 80, "CD169 macrophage; PMID:18024188"),
    "SA SIGLEC1": ({"M1", "M2"}, 75, "Sialic acid/Siglec-1"),

    # Drop / unmodeled
    "IP6K2":   (ALL, 65, "IP6 kinase"),
    "SESN1":   (ALL, 65, "Sestrin stress"),
    "SHISA5":  (ALL, 65, "Apoptosis regulator"),
    "SOCS2":   (ALL, 70, "JAK-STAT inhibitor"),
}


# ---------------------------------------------------------------------------
# Application sur un dict d'assignations
# ---------------------------------------------------------------------------


def _normalise(name: str) -> str:
    """Match key : majuscules + espaces collapse + ponctuation→space."""
    s = (name or "").upper().strip()
    s = re.sub(r"[^A-Z0-9/\- ]", " ", s)
    s = re.sub(r"[/\-]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def lookup(name: str) -> tuple[set[str], int, str] | None:
    """Lookup curated avec normalisation."""
    if not name:
        return None
    key = _normalise(name)
    if key in CURATED:
        return CURATED[key]
    # try uppercase raw
    raw = name.upper().strip()
    if raw in CURATED:
        return CURATED[raw]
    return None


def refine_assignments(assignments: dict[int, Any]) -> dict[str, int]:
    """
    Modifie en place les Assignment R6 dont le nom matche CURATED.

    Returns:
        stats {n_refined, n_dropped, n_no_match, n_unchanged_non_r6}
    """
    n_refined = 0
    n_dropped = 0
    n_no_match = 0
    n_non_r6 = 0
    for nid, a in assignments.items():
        if a.rule != "R6":
            n_non_r6 += 1
            continue
        hit = lookup(a.name)
        if hit is None:
            n_no_match += 1
            continue
        cts, score, src = hit
        if not cts:
            # Marker explicite "DROP" — basculer en R7
            a.celltypes = set()
            a.rule = "R7"
            a.confidence = ""
            a.evidence = "expert_curated_drop"
            a.score_per_celltype = {}
            a.source_per_celltype = {}
            n_dropped += 1
            continue
        a.celltypes = set(cts)
        a.rule = "R6c"
        a.confidence = "MEDIUM" if score >= 75 else "LOW"
        a.evidence = f"curated_from_R6;orig_evidence={a.evidence}"
        a.score_per_celltype = {ct: score for ct in cts}
        a.source_per_celltype = {ct: src for ct in cts}
        n_refined += 1
    logger.info(
        "Raffinement : %d refinés, %d droppés (R7), %d sans match (R6 conservé), "
        "%d non-R6 inchangés",
        n_refined, n_dropped, n_no_match, n_non_r6,
    )
    return {
        "n_refined": n_refined,
        "n_dropped": n_dropped,
        "n_no_match": n_no_match,
        "n_unchanged_non_r6": n_non_r6,
        "curated_entries": len(CURATED),
    }
