# CONVENTIONS — SjS-DigitalTwin

> Conventions de nommage, organisation des fichiers, règles de codage et de commit.
> **Statut** : VERROUILLÉ — toute modification implique mise à jour synchronisée des scripts dépendants.

---

## 1. Nommage des espèces (nodes)

### 1.1 Identifiants externes (priorité décroissante)

Chaque nœud DOIT être annoté avec ≥1 identifiant externe lorsqu'il existe :

| Type biologique | Source canonique | Exemple |
|---|---|---|
| Gène / protéine humaine | **HGNC symbol** (officiel) | `IFNG`, `BAFF`, `BTK`, `STAT1` |
| Métabolite | **ChEBI** ID | `CHEBI:15422` (ATP) |
| Phénotype | **GO** ID + libellé interne | `GO:0006915` → `phen_apoptosis` |
| Complexe protéique | Concaténation HGNC `_` séparé | `IFNAR1_IFNAR2` |
| Famille (groupes) | Préfixe `fam_` + nom court | `fam_ISGs`, `fam_TLRs` |

### 1.2 Suffixes cell-type

Quand un nœud existe dans plusieurs types cellulaires, on **clone** avec suffixe :

| Suffixe | Signification |
|---|---|
| `_SGEC` | Salivary gland epithelial cell (acinaire/ductale) |
| `_TH1` | CD4+ T helper 1 |
| `_TH17` | CD4+ T helper 17 |
| `_TFH` | T follicular helper |
| `_TREG` | T regulatory |
| `_BCELL` | B cell naïve / mémoire / GC |
| `_PLASMA` | Plasmablaste / plasmocyte |
| `_M1` | Macrophage M1 (pro-inflammatoire) |
| `_M2` | Macrophage M2 (réparation / résolution) |
| `_PDC` | Plasmacytoid dendritic cell |
| `_EXTRA` | Espèce extracellulaire partagée (ligands sécrétés, immunoglobulines) |

**Exemples** : `STAT1_SGEC`, `STAT1_TH17`, `TNFSF13B_EXTRA` (BAFF sécrété), `CD40_BCELL` vs `CD40LG_TH17`.

**Règle** : un ligand sécrété est UNIQUE et appartient à `_EXTRA` ; ses récepteurs sont cell-type-spécifiques.

### 1.3 Phénotypes (nœuds de sortie)

Toujours préfixés `phen_` :

```
phen_apoptosis_SGEC
phen_secretory_loss
phen_inflammation
phen_lymphocytic_infiltration
phen_tls_formation
phen_b_cell_activation
phen_t_cell_activation
phen_lymphomagenesis
phen_fibrosis
phen_mhc1_activation
phen_mhc2_activation
phen_phagocytosis
phen_regulated_necrosis
phen_angiogenesis
phen_matrix_degradation
phen_proliferation
phen_chemotaxis
```

(14 phénotypes principaux issus du papier SjD Map + ajouts SjS-spécifiques `secretory_loss`, `tls_formation`, `lymphocytic_infiltration`).

---

## 2. Compartiments

| ID compartiment | Signification |
|---|---|
| `extracellular` | Milieu extracellulaire commun à tous les cell-types |
| `cytosol_<celltype>` | Cytoplasme d'un cell-type donné |
| `nucleus_<celltype>` | Noyau |
| `membrane_<celltype>` | Membrane plasmique |
| `endosome_<celltype>` | Compartiment endosomal (TLR7/9, MHC-II loading) |
| `mitochondrion_<celltype>` | Mitochondrie (apoptose intrinsèque) |
| `granule_<celltype>` | Granules (cytotoxicité, sécrétion) |

Le compartiment `Default` produit par CaSQ est **interdit** dans la carte multi-cellulaire finale — toute espèce doit être assignée.

---

## 3. Identifiants SBML / CellDesigner

### 3.1 Espèces

Format : `s_<celltype>_<symbol>` ou `s_extra_<symbol>`.

```
s_sgec_STAT1
s_th17_RORC
s_extra_BAFF
s_bcell_CD20
```

- Préfixe `s_` obligatoire (CellDesigner exige ID alphanumérique commençant par lettre)
- Casse : tout en minuscules sauf le symbole HGNC (gardé en MAJUSCULE)
- Underscore comme séparateur unique (pas de tiret, pas de point)

### 3.2 Réactions

Format : `r_<type>_<celltype>_<n>` où `<type>` ∈ {act, inh, trans, complex, phen, transport}.

```
r_act_th17_001       # Activation Th17 réaction #001
r_inh_sgec_042       # Inhibition SGEC réaction #042
r_transport_023      # Transport extracellulaire réaction #023
r_complex_bcell_007  # Formation de complexe B cell réaction #007
r_phen_017           # Réaction de phénotype #017
```

### 3.3 Compartiments

Format : ID = `c_<name>` ; `name` exact selon §2.

```
c_extracellular
c_cytosol_sgec
c_nucleus_th17
```

---

## 4. Organisation des fichiers

### 4.1 Code source

- `scripts/` : scripts d'orchestration numérotés `NN_<verbe>_<objet>.py` exécutables en CLI
- `scripts/lib/` : modules réutilisables (pas exécutables directement)
- `tests/` : tests pytest, miroir de l'arborescence `scripts/`
- `notebooks/` : notebooks Jupyter de QC / exploration (suffixés `.ipynb`)

### 4.2 Données

- `01_disease_map/` à `06_results/` selon ROADMAP
- `01_disease_map/cache/` : caches JSON API MINERVA (gitignored sauf petits fichiers de fixtures)
- `data/` : données brutes externes (transcriptomes GEO) — **gitignored**, récupérées via script
- Artéfacts > 1 Mo versionnés via Git LFS (carte XML, modèle Booléen)

### 4.3 Documentation

- `README.md` : présentation projet
- `ROADMAP.md` : plan détaillé + suivi avancement
- `CONVENTIONS.md` : ce fichier
- `journal.md` : log chronologique des actions
- `CLAUDE.md` : instructions assistant

---

## 5. Règles de codage Python

- **Style** : PEP 8 + ruff (`ruff check` + `ruff format`)
- **Type hints** : obligatoires sur toute fonction publique de `scripts/lib/`, `mypy --strict`
- **Docstrings** : Google style, en français pour cohérence avec le reste du projet
- **Imports** : `from __future__ import annotations` en tête de chaque module
- **Logging** : `logging.getLogger(__name__)`, pas de `print()` dans lib/
- **CLI** : `argparse` avec sous-commandes si pertinent ; `--verbose` pour DEBUG
- **Pas de globals mutables** : utiliser `dataclasses` ou paramètres explicites
- **Chemins** : `pathlib.Path` jamais `os.path`

---

## 6. Règles de codage R

- **Style** : tidyverse style guide
- **Pipe** : `|>` (R ≥ 4.1) plutôt que `%>%`
- **Scripts numérotés** dans `04_validation/scripts/NN_<verbe>.R`
- **Sessioninfo** : chaque script termine par `sessionInfo()` capturé dans un log

---

## 7. Git & commits

- Branches : `main` protégée ; développement sur `feat/<topic>` ou `phase-X.Y/<sub>`
- Commits : conventional commits (`feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`)
- Tags : `v0.1` après gate Phase 0, `v1.0` à soumission manuscrit
- Messages en anglais (compatibilité GitHub) ; commentaires de code en français OK
- Pas de commit de secrets, données patient, fichiers binaires non-LFS

---

## 8. Reproductibilité

- **Environnement scellé** : `envs/environment.yml` + `Dockerfile` + tag d'image
- **Seeds** : tout RNG (numpy, R `set.seed`, MaBoSS) reçoit un seed explicite ; default `42`
- **Caches** : tout cache externe (MINERVA, GEO) versionné via empreinte SHA256 dans le code
- **Versions des artéfacts** : naming explicite `<artifact>_v<N>_<YYYYMMDD>.xml` quand pertinent

---

## 9. Validation scientifique

Toute affirmation biologique dans le code, doc ou commit DOIT être :
- Soit traçable à un PMID
- Soit annotée comme `# HYPOTHESIS` (revue à venir)
- Soit annotée comme `# COMPUTATIONAL` (déduction du modèle, à valider)

---

*Dernière mise à jour : 2026-04-25*
