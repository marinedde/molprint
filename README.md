# MolPrint

Projet in silico de prédiction de réponse aux médicaments, construit en complément d'[OncoPrint](https://github.com/marinedde/cdsd-certification/tree/main/bloc6-direction-projet/oncoprint) (classification moléculaire des sous-types de cancer du sein, TCGA-BRCA, XGBoost + SHAP).

Objectif : couvrir les compétences bioinformatique/chémoinformatique/biologie des systèmes demandées sur les postes Data Scientist Préclinique en R&D pharmaceutique (cible marché : Suisse, Bâle).

## Roadmap

| Phase | Sujet | Outils | Statut |
| --- | --- | --- | --- |
| 1 | Séquences ARN/ADN | Biopython | à venir |
| 2 | Chémoinformatique (screening in silico, QSAR) | RDKit, PubChem BioAssay, XGBoost | **en cours** |
| 3 | Biologie des systèmes (modèle ODE d'une voie de signalisation) | Tellurium / PySB | à venir |
| 4 | Base de données interne (patientes ↔ gènes ↔ molécules ↔ voies) | SQLite | à venir |
| 5 | Assemblage (dashboard Streamlit) + veille bibliographique | - | à venir |

## Phase 2 — Chémoinformatique

Cible thérapeutique : **HER2/ERBB2** (choisie car elle relie directement au sous-type "HER2-enrichi" d'OncoPrint — le trastuzumab est un traitement ciblé de référence).

1. `notebooks/01_bioactivity_data_acquisition.ipynb` — récupération des molécules testées sur ERBB2 via l'API PubChem BioAssay (table "bioactivity concise"). ChEMBL était visée initialement mais son API était en panne (erreurs 500) au moment de l'écriture.
2. `notebooks/02_rdkit_descriptors_qsar.ipynb` — calcul des descripteurs moléculaires (RDKit, règles de Lipinski) et modèle QSAR (XGBoost + SHAP) prédisant l'activité d'une molécule à partir de sa structure

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
