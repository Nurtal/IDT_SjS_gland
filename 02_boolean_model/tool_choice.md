# Phase 2.2 — Choix d'outil pour l'énumération d'attracteurs

> Comparatif réalisé sur le modèle Booléen `SjS_boolean.bnet`
> (5 015 nœuds, 2 435 règles non triviales — Phase 2.1).

## Résultats du POC

| Outil | Statut | Temps | Sortie | Échelle observée |
|---|---|---:|---|---|
| **mpbn** (4.3.2) | ✅ OK | 210 s | 5 000 trap-spaces minimaux (cap atteint) | 5 015 nœuds tractable |
| **pyboolnet** (3.0.16) | ❌ Intractable | >180 s sur la seule conversion `bnet2primes` | Aucune | `BNetToPrime` ne scale pas au-delà de ~ 1 000 nœuds |
| **MaBoSS** | Non testé Phase 2.2 | n/a | n/a | À évaluer Phase 2.5 si la dynamique temporelle devient nécessaire |
| **BMA** | Non testé Phase 2.2 | n/a | n/a | Méthode Zerrouk 2024 — env. Java/HPC requis (reporté) |
| **bioLQM** | Non testé Phase 2.2 | n/a | n/a | Réduction formelle — pertinent pour rendre pyboolnet/BMA tractables |

## Détail mpbn

- Charge le `.bnet` normalisé en 0.3 s.
- Énumère 5 000 trap-spaces minimaux en 3.5 minutes (cap arbitraire — l'énumération continue au-delà).
- **Tous les trap-spaces sont des point fixes complets** (`free_dim=0`) — pas de wildcards `*`.
  - Implication forte : le réseau présente un grand nombre de point fixes asynchrones, ce qui est cohérent avec un modèle multi-cellulaire où chaque cell-type a sa propre dynamique stable, et où les nœuds-inputs (51 % du modèle, voir Phase 2.1) figent une grande partie de l'espace d'états.

## Détail pyboolnet

- `pyboolnet.file_exchange.bnet2primes` ne termine pas après 180 s.
  - La conversion repose sur `BNetToPrime` (binaire externe) dont la complexité explose en présence de gros disjonctions OR.
  - Test indépendant confirmé : `bnet2primes` seul (sans `compute_steady_states`) hang aussi.
- Cross-check des point fixes mpbn ↔ pyboolnet **non réalisable à cette échelle** sans réduction préalable.

## Décision

- **mpbn retenu pour Phase 2.3** (énumération + filtrage par phénotype).
- Cross-check pyboolnet **reporté à Phase 2.1bis**, après une étape de **réduction bioLQM** (qui devrait ramener le modèle à <1 500 nœuds en éliminant les nœuds tampons et les inputs non régulants).
- BMA reporté à Phase 2.5 / discussion finale — non bloquant pour le pipeline scientifique.
- MaBoSS sera utilisé Phase 2.5 si l'on a besoin de simulations stochastiques (perturbations populationnelles, robustesse).

## Conséquences pour Gate 2.2

| Critère Gate 2.2 | Statut |
|---|---|
| ≥ 1 outil produit un set d'attracteurs en < 24 h | ✅ mpbn (3.5 min pour 5 000 trap-spaces) |
| Cohérence ≥ 2 outils sur les point fixes | ⚠️ Reportée à 2.1bis (pyboolnet intractable sans réduction) |
| Document `tool_choice.md` justifiant l'outil retenu | ✅ Ce fichier |

**Gate 2.2 : PASS partiel** — mpbn validé, cross-check inter-outils reporté.

## Livrables Phase 2.2

- `02_boolean_model/poc_results/mpbn/SjS_boolean_normalized.bnet` — `.bnet` réécrit avec identifiants safe (espaces et caractères spéciaux normalisés).
- `02_boolean_model/poc_results/mpbn/name_map.tsv` — mapping `nom CaSQ` → `nom safe`.
- `02_boolean_model/poc_results/mpbn/trap_spaces.tsv` — 5 000 lignes × 5 015 colonnes.
- `02_boolean_model/poc_results/mpbn/summary.json`.
- `02_boolean_model/poc_results/pyboolnet/summary.json` — preuve du timeout.
- `scripts/08_mpbn_attractors.py`, `scripts/09_pyboolnet_fixedpoints.py`.
