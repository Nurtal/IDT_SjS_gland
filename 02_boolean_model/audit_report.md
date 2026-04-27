# Phase 2.1 — Audit du modèle Booléen CaSQ

- Fichier `.bnet` : `02_boolean_model/casq_output/SjS_boolean.bnet`
- Fichier transitions : `02_boolean_model/casq_output/SjS_boolean_Transitions.csv`

## Compteurs globaux

| Métrique | Valeur |
|---|---|
| Nœuds totaux | 5015 |
| Nœuds avec règle non triviale | 2435 (48.6 %) |
| Nœuds inputs (X = X) | 2580 (51.4 %) |
| Lignes Transitions.csv | 2345 |
| Opérateurs AND `&` | 665 |
| Opérateurs OR `\|` | 1191 |
| Opérateurs NOT `!` | 304 |

## Répartition par cell-type

| Cell-type | Total nœuds | Non triviaux | % non triviaux |
|---|---:|---:|---:|
| SGEC | 440 | 197 | 44.8 |
| TH1 | 469 | 226 | 48.2 |
| TH17 | 455 | 204 | 44.8 |
| TFH | 443 | 201 | 45.4 |
| TREG | 441 | 197 | 44.7 |
| BCELL | 494 | 280 | 56.7 |
| PLASMA | 458 | 246 | 53.7 |
| M1 | 530 | 241 | 45.5 |
| M2 | 501 | 218 | 43.5 |
| PDC | 464 | 212 | 45.7 |
| EXTRA | 49 | 13 | 26.5 |
| OTHER | 271 | 200 | 73.8 |

## Top 20 nœuds (in-degree)

| Rang | Nœud | Régulateurs |
|---:|---|---:|
| 1 | Inflammation_phenotype | 220 |
| 2 | Apoptosis_phenotype | 60 |
| 3 | MHC_Class_1_Activation_phenotype | 48 |
| 4 | MHC_Class_2_Activation_phenotype | 40 |
| 5 | T_Cell_Activation/Differentiation_phenotype | 20 |
| 6 | CHUK/IKBKB/IKBKG_complex_BCELL Cytoplasm | 16 |
| 7 | CHUK/IKBKB/IKBKG_complex_PLASMA Cytoplasm | 16 |
| 8 | JAK1_BCELL Cytoplasm | 13 |
| 9 | JAK1_M1 Cytoplasm | 13 |
| 10 | JAK1_M2 Cytoplasm | 13 |
| 11 | JAK1_PDC Cytoplasm | 13 |
| 12 | JAK1_PLASMA Cytoplasm | 13 |
| 13 | JAK1_SGEC Cytoplasm | 13 |
| 14 | JAK1_TFH Cytoplasm | 13 |
| 15 | JAK1_TH17 Cytoplasm | 13 |
| 16 | JAK1_TH1 Cytoplasm | 13 |
| 17 | JAK1_TREG Cytoplasm | 13 |
| 18 | Phagocytosis_phenotype | 12 |
| 19 | B_Cell_Activation/Survival_phenotype | 12 |
| 20 | Regulated_Necrosis_phenotype | 10 |

## Critères de Gate Phase 2.1

| Critère | Statut |
|---|---|
| ≥ 40 % nœuds non triviaux | PASS |
| ≥ 1 règle AND | PASS |
| ≥ 1 règle OR | PASS |
| ≥ 1 règle NOT | PASS |
| Couverture ≥ 8 / 10 cell-types | PASS |

**Gate Phase 2.1 : PASS**
