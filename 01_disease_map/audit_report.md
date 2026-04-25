# Audit topologique & sémantique — SjD Map

> **Phase 1.1 ROADMAP** — caractérisation de la SjD Map avant dissociation cell-type (Phase 1.2).

## 1. Topologie

| Métrique | Valeur |
|---|---|
| Nombre de nœuds | 1157 |
| Nombre d'edges | 1151 |
| Composantes faibles | 316 |
| Plus grosse composante faible | 842 |
| Composantes fortes | 1121 |
| Plus grosse composante forte | 23 |
| Nœuds isolés (degree 0) | 315 |
| Nœuds sources (in=0) | 371 |
| Nœuds puits (out=0) | 36 |
| Hubs (degree ≥ 20) | 6 |

### Top 20 hubs

| id | name | type | in | out | total |
|---|---|---|---|---|---|
| 21541 | Inflammation | Phenotype | 44 | 0 | 44 |
| 21240 | STAT1 homodimer | Protein | 1 | 42 | 43 |
| 21485 | STAT1/STAT2/IRF9 | Complex | 2 | 40 | 42 |
| 21272 | RELA/NFKB1 | Complex | 1 | 27 | 28 |
| 20805 | GNAI | Protein | 20 | 2 | 22 |
| 21544 | Chemotaxis/Infiltration | Phenotype | 22 | 0 | 22 |

## 2. Couverture des annotations

**Éléments** : 1157

| Annotation | Nb éléments |
|---|---|
| GO | 1010 |
| Entrez | 890 |
| UniProt | 890 |
| HGNC | 890 |
| Ensembl | 888 |
| KEGG | 369 |
| DOI | 281 |
| ChEBI | 17 |

- Avec ≥1 annotation externe : **1090** (94.21%)
- Sans aucune annotation : **67**
- Avec PubMed ID : **508**

**Réactions** : 598
- Avec PubMed ID : **100** (16.72%)

### Gate 1.1

- ≥80% éléments avec annotation : **PASS** (94.21%)
- ≥50% réactions avec PMID : **FAIL** (16.72%)

## 3. Compartiments

| compartmentId | n_species | top types |
|---|---|---|
| 20513 | 624 | Protein, Complex, Simple molecule, Compartment |
| 21231 | 308 | RNA, Gene, Protein, Complex |
| None | 93 | Protein, Simple molecule, Compartment |
| 21555 | 73 | Protein, Simple molecule, Complex |
| 21629 | 40 | Protein |
| 21540 | 14 | Phenotype |
| 20730 | 5 | Protein, Complex, Simple molecule, Ion |

## 4. Marqueurs cell-type détectés

Pour chaque cell-type cible, nombre de marqueurs canoniques (cf. CONVENTIONS.md) trouvés par correspondance de nom :

| Cell-type | trouvés | attendus | % | exemples |
|---|---|---|---|---|
| SGEC | 0 | 19 | 0.0% |  |
| TH1 | 15 | 7 | 214.3% | IFNG, CD3G, IFNG, CD3E, IFNG, IFNG |
| TH17 | 10 | 8 | 125.0% | IL21, IL17A, IL17A, IL17A, STAT3, IL21 |
| TFH | 6 | 4 | 150.0% | PDCD1, CXCR5, ICOS, CXCR5, ICOS, PDCD1 |
| TREG | 7 | 4 | 175.0% | IL2RA, IL2RA, CTLA4, IL2RA, IL2RA, CTLA4 |
| BCELL | 10 | 11 | 90.9% | CD79A, BTK, TNFRSF13B, TNFRSF13C, TNFRSF13C, BTK |
| PLASMA | 2 | 7 | 28.6% | TNFRSF17, TNFRSF17 |
| M1 | 24 | 8 | 300.0% | IL1B, FCGR1A, FCGR1A, TNF, CXCL10, CXCL10 |
| M2 | 5 | 7 | 71.4% | TGFB1, TGFB1, IL10, IL10, IL10 |
| PDC | 9 | 6 | 150.0% | IRF7, IRF7, TLR7, IRF7, TLR7, IRF7 |

## 5. Pivots intercellulaires (cytokines/chemokines clés)

**92 / 36** cytokines d'intérêt SjD trouvées dans la carte.

| id | name | type | matched |
|---|---|---|---|
| 21657 | CCL21 | Protein | CCL21 |
| 21667 | TNFSF13B | Protein | TNFSF13B |
| 21666 | IL1B | Protein | IL1B |
| 21290 | CXCL9 | Gene | CXCL9 |
| 21209 | IFNG | Protein | IFNG |
| 21633 | IL23A | Protein | IL23A |
| 21575 | LTA | Protein | LTA |
| 21615 | IFNB1 | Protein | IFNB1 |
| 21329 | IL21 | RNA | IL21 |
| 21310 | IL15 | RNA | IL15 |
| 21662 | IFNG | Protein | IFNG |
| 21571 | TNF | Protein | TNF |
| 20744 | LTA | Protein | LTA |
| 21262 | CXCL10 | RNA | CXCL10 |
| 21658 | CXCL10 | Protein | CXCL10 |
| 21254 | IL17A | RNA | IL17A |
| 21265 | IFNB1 | Gene | IFNB1 |
| 21562 | IL15 | Protein | IL15 |
| 21266 | IL6 | Gene | IL6 |
| 21645 | CXCL9 | Protein | CXCL9 |
| 21567 | CCL19 | Protein | CCL19 |
| 21637 | CCL4 | Protein | CCL4 |
| 21580 | CXCL13 | Protein | CXCL13 |
| 21341 | CCL3 | RNA | CCL3 |
| 21263 | IFNB1 | RNA | IFNB1 |
| 21622 | CXCL11 | Protein | CXCL11 |
| 21649 | IL17A | Protein | IL17A |
| 21375 | IL23A | Gene | IL23A |
| 21373 | CCL5 | Gene | CCL5 |
| 20680 | TNFSF13B | Protein | TNFSF13B |
| 21023 | IL6 | Protein | IL6 |
| 21277 | CCL4 | Gene | CCL4 |
| 21289 | IL17A | Gene | IL17A |
| 21490 | CXCL10 | Gene | CXCL10 |
| 21565 | TNFSF14 | Protein | TNFSF14 |
| 21642 | TNF | Protein | TNF |
| 21253 | IFNG | Gene | IFNG |
| 21476 | IFNG | RNA | IFNG |
| 20985 | CCL19 | Protein | CCL19 |
| 21625 | IL6 | Protein | IL6 |
| 21647 | IL15 | Protein | IL15 |
| 21267 | CXCL11 | RNA | CXCL11 |
| 21659 | CCL19 | Protein | CCL19 |
| 21579 | CXCL9 | Protein | CXCL9 |
| 21338 | CCL3 | Gene | CCL3 |
| 21246 | CXCL11 | Gene | CXCL11 |
| 21623 | LTB | Protein | LTB |
| 21651 | CCL5 | Protein | CCL5 |
| 21237 | IL1B | Gene | IL1B |
| 21445 | CCL2 | Gene | CCL2 |
| 21344 | IL6 | RNA | IL6 |
| 21089 | CCL2 | Protein | CCL2 |
| 21199 | CXCL13 | Protein | CXCL13 |
| 21193 | IL15 | Protein | IL15 |
| 21609 | CCL21 | Protein | CCL21 |
| 21443 | CCL19 | Gene | CCL19 |
| 21635 | IFNB1 | Protein | IFNB1 |
| 21632 | CXCL11 | Protein | CXCL11 |
| 21669 | CCL3 | Protein | CCL3 |
| 21036 | IL21 | Protein | IL21 |
| 21337 | IL15 | Gene | IL15 |
| 21636 | IL21 | Protein | IL21 |
| 21120 | CCL5 | Protein | CCL5 |
| 21241 | CCL4 | RNA | CCL4 |
| 21655 | IL6 | Protein | IL6 |
| 20688 | TNF | Protein | TNF |
| 21380 | CCL21 | Gene | CCL21 |
| 21535 | CCL19 | RNA | CCL19 |
| 21239 | CCL5 | RNA | CCL5 |
| 21648 | IL10 | Protein | IL10 |
| 21442 | TNF | Gene | TNF |
| 21480 | TNFSF13B | RNA | TNFSF13B |
| 21446 | IL1B | RNA | IL1B |
| 21400 | IL21 | Gene | IL21 |
| 21469 | IL10 | RNA | IL10 |
| 21311 | CCL2 | RNA | CCL2 |
| 21577 | IL21 | Protein | IL21 |
| 21559 | CCL5 | Protein | CCL5 |
| 21275 | TNFSF13B | Gene | TNFSF13B |
| 21572 | CCL2 | Protein | CCL2 |
| 21538 | IL10 | Gene | IL10 |
| 21410 | CCL21 | RNA | CCL21 |
| 21641 | CCL2 | Protein | CCL2 |
| 21368 | IL23A | RNA | IL23A |
| 21318 | CXCL9 | RNA | CXCL9 |
| 21435 | TNF | RNA | TNF |
| 21578 | CXCL10 | Protein | CXCL10 |
| 21223 | IFNB1 | Protein | IFNB1 |
| 21581 | TNFSF13 | Protein | TNFSF13 |
| 21585 | IFNG | Protein | IFNG |
| 21147 | TNF | Protein | TNF |
| 21563 | TNFSF13B | Protein | TNFSF13B |

**Cytokines clés manquantes** (11) : CCL22, IFNA1, IFNA13, IFNA2, IFNL1, IFNL2, IFNL3, IL12A, IL12B, IL17F, IL22

## 6. Indices de pathway extraits des notes

**235** éléments avec ≥1 identifiant Reactome / KEGG / GO dans leurs notes.

## 7. Recommandations pour Phase 1.2 (dissociation)

- ✅ **6 compartiments** présents → approche compartiment-driven viable comme première passe.
- ✅ Marqueurs bien représentés : TH1, TH17, TFH, TREG, BCELL, M1, M2, PDC
- ⚠️ Marqueurs peu représentés : SGEC → envisager extension littérature ou retrait du cell-type.

---

_Généré par `scripts/02_audit_map.py`_