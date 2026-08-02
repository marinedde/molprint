# MolPrint

Projet in silico de prédiction de réponse aux médicaments, construit en complément d'[OncoPrint](https://github.com/marinedde/cdsd-certification/tree/main/bloc6-direction-projet/oncoprint) (classification moléculaire des sous-types de cancer du sein, TCGA-BRCA, XGBoost + SHAP).

Objectif : couvrir les compétences bioinformatique/chémoinformatique/biologie des systèmes demandées sur les postes Data Scientist Préclinique en R&D pharmaceutique (cible marché : Suisse, Bâle).

## Roadmap

| Phase | Sujet | Outils | Statut |
| --- | --- | --- | --- |
| 1 | Séquences ARN/ADN | Biopython | **terminée** |
| 2 | Chémoinformatique (screening in silico, QSAR) | RDKit, PubChem BioAssay, XGBoost | **terminée** |
| 3 | Biologie des systèmes (modèle ODE d'une voie de signalisation) | Tellurium / PySB | à venir |
| 4 | Base de données interne (patientes ↔ gènes ↔ molécules ↔ voies) | SQLite | à venir |
| 5 | Assemblage (dashboard Streamlit) + veille bibliographique | - | à venir |

## Phase 2 — Chémoinformatique

Cible thérapeutique : **HER2/ERBB2** (choisie car elle relie directement au sous-type "HER2-enrichi" d'OncoPrint — le trastuzumab est un traitement ciblé de référence).

1. `notebooks/01_bioactivity_data_acquisition.ipynb` — récupération des molécules testées sur ERBB2 via l'API PubChem BioAssay (table "bioactivity concise"). ChEMBL était visée initialement mais son API était en panne (erreurs 500) au moment de l'écriture.
2. `notebooks/02_rdkit_descriptors_qsar.ipynb` — calcul des descripteurs moléculaires (RDKit, règles de Lipinski) et modèle QSAR (XGBoost + SHAP) prédisant l'activité d'une molécule à partir de sa structure
3. `notebooks/03_virtual_screening_optimization.ipynb` — criblage virtuel et optimisation in silico : validation du modèle sur 5 médicaments anti-HER2 approuvés (probabilité prédite cohérente avec leur statut d'inhibiteurs connus), génération de nouveaux candidats par recombinaison de fragments (BRICS) à partir des molécules les plus actives, filtre de diversité structurale (clustering Butina sur empreintes de Morgan) et sélection d'un top 15
4. `notebooks/04_fingerprint_qsar.ipynb` — comparaison en validation croisée de trois représentations moléculaires (descripteurs seuls, empreintes de Morgan seules, combinaison des deux) ; le modèle combiné (ROC AUC 0,98) devient le modèle final, nettement plus fiable sur les médicaments de référence que le modèle à descripteurs seuls

## Phase 1 — Séquences ARN/ADN

Gènes choisis : les 5 features les plus importantes en SHAP dans OncoPrint, un par sous-type de cancer du sein (ESR1, ERBB2, FOXA1, AR, GATA3) — ancrage direct sur les résultats déjà produits plutôt que des gènes choisis arbitrairement.

`notebooks/05_gene_sequence_analysis.ipynb` (fonctions dans `src/sequence_analysis.py`) :
- Récupération des séquences RefSeq canoniques (NCBI, `Bio.Entrez`), composition GC, traduction ADN → protéine, descripteurs physico-chimiques de la protéine (poids moléculaire, indice d'instabilité, point isoélectrique)
- Vérification, directement dans la séquence traduite, que la position 537 d'ESR1 porte bien une tyrosine — le site exact de la mutation de résistance à l'hormonothérapie **Y537S**, l'une des mieux documentées dans le cancer du sein métastatique — puis simulation de la mutation et comparaison des propriétés physico-chimiques avant/après

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Structure

```
molprint/
├── data/
│   ├── raw/          # données brutes (PubChem, ChEMBL...) — non versionnées
│   └── processed/    # jeux de données nettoyés
├── notebooks/         # notebooks d'exploration, un par étape
├── src/               # fonctions réutilisables (extraites des notebooks)
├── models/            # modèles entraînés (non versionnés)
└── docs/              # fiches de veille bibliographique
```
