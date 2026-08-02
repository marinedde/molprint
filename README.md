# MolPrint

**[🧪 Dashboard live](https://huggingface.co/spaces/marinedde/molprint-dashboard)** — cribler une molécule, explorer les candidats sélectionnés, voir les gènes analysés et la validation méthodologique.

**[📖 Guide complet](docs/GUIDE_COMPLET.md)** — tout le projet expliqué pas à pas (analogies, code, résultats, SQL/Python/ML pour les nuls), plus un audit technique/architecture complet et une vérification data leakage/overfitting.

Projet in silico de prédiction de réponse aux médicaments, construit en complément d'[OncoPrint](https://github.com/marinedde/cdsd-certification/tree/main/bloc6-direction-projet/oncoprint) (classification moléculaire des sous-types de cancer du sein, TCGA-BRCA, XGBoost + SHAP).

Objectif : couvrir les compétences bioinformatique/chémoinformatique/biologie des systèmes demandées sur les postes Data Scientist Préclinique en R&D pharmaceutique (cible marché : Suisse, Bâle).

## Roadmap

| Phase | Sujet | Outils | Statut |
| --- | --- | --- | --- |
| 1 | Séquences ARN/ADN | Biopython | **terminée** |
| 2 | Chémoinformatique (screening in silico, QSAR) | RDKit, PubChem BioAssay, XGBoost | **terminée** |
| 3 | Biologie des systèmes (modèle ODE d'une voie de signalisation) | Tellurium | **terminée** |
| 4 | Base de données interne (gènes ↔ molécules ↔ voies ↔ sous-types) | SQLite, GDSC | **terminée** |
| 5 | Assemblage (dashboard Streamlit) + veille bibliographique | Streamlit, HuggingFace Spaces | **dashboard en ligne** |

## Phase 2 — Chémoinformatique

Cible thérapeutique : **HER2/ERBB2** (choisie car elle relie directement au sous-type "HER2-enrichi" d'OncoPrint — le trastuzumab est un traitement ciblé de référence).

1. `notebooks/01_bioactivity_data_acquisition.ipynb` — récupération des molécules testées sur ERBB2 via l'API PubChem BioAssay (table "bioactivity concise"). ChEMBL était visée initialement mais son API était en panne (erreurs 500) au moment de l'écriture.
2. `notebooks/02_rdkit_descriptors_qsar.ipynb` — calcul des descripteurs moléculaires (RDKit, règles de Lipinski) et modèle QSAR (XGBoost + SHAP) prédisant l'activité d'une molécule à partir de sa structure
3. `notebooks/03_virtual_screening_optimization.ipynb` — criblage virtuel et optimisation in silico : validation du modèle sur 5 médicaments anti-HER2 approuvés (probabilité prédite cohérente avec leur statut d'inhibiteurs connus), génération de nouveaux candidats par recombinaison de fragments (BRICS) à partir des molécules les plus actives, filtre de diversité structurale (clustering Butina sur empreintes de Morgan) et sélection d'un top 15
4. `notebooks/04_fingerprint_qsar.ipynb` — comparaison en validation croisée de trois représentations moléculaires (descripteurs seuls, empreintes de Morgan seules, combinaison des deux) ; le modèle combiné (ROC AUC 0,98) devient le modèle final, nettement plus fiable sur les médicaments de référence que le modèle à descripteurs seuls
5. `notebooks/06_scaffold_validation_applicability_domain.ipynb` — **rigueur méthodologique** : le split aléatoire (ROC AUC 0,98) est optimiste car le jeu de données contient beaucoup d'analogues proches (499 squelettes chimiques pour 1070 molécules). Un split par squelette de Bemis-Murcko (méthode standard du domaine, type MoleculeNet) donne un ROC AUC honnête de 0,91 sur des squelettes jamais vus (rappel actif 0,59) — et un contrôle de domaine d'applicabilité (similarité de Tanimoto au train) montre que les candidats générés en notebook 03 restent dans une zone "modérément nouvelle", donc leurs scores sont des pistes à explorer plutôt que des prédictions fiables

## Phase 4 — Base de données interne

Relie dans une seule base SQLite (`data/processed/molprint.db`) les briques construites séparément : gènes (Phase 1), molécules candidates (Phase 2), voie de signalisation (Phase 3), et une source externe indépendante — **GDSC** (Genomics of Drug Sensitivity in Cancer, Sanger Institute), qui mesure la sensibilité réelle en laboratoire de 51 lignées cellulaires de cancer du sein à des centaines de médicaments.

`notebooks/08_internal_database.ipynb` :
- 7 tables : `genes`, `subtypes`, `pathways`, `pathway_genes`, `molecules`, `cell_lines`, `drug_response`
- Annotation des 51 lignées GDSC par sous-type moléculaire (classification manuelle issue de la littérature — Neve et al. 2006, Kao et al. 2009 — volontairement incomplète : les lignées non documentées restent "Non classé" plutôt que d'être devinées)
- Requête roadmap : *"pour un sous-type donné, quels gènes clés, quelles molécules candidates, quelle voie de signalisation ?"*
- **Validation empirique indépendante** : les données GDSC (réelles, sans lien avec nos modèles) confirment que le lapatinib est ~10x plus puissant sur les lignées HER2-enrichi (IC50 moyen ≈ 1,8 µM) que sur Luminal A ou Triple Négatif (≈ 17-21 µM) — cohérent avec les hypothèses des Phases 2 et 3

## Phase 3 — Biologie des systèmes

Voie de signalisation **HER2 → PI3K → AKT** (le cœur de la prolifération dans le sous-type "HER2-enrichi"), modélisée par équations différentielles ordinaires (ODE) avec Tellurium — cascade de kinases à 3 étages, chaque protéine passant d'un état inactif à actif.

`notebooks/07_systems_biology_her2_pathway.ipynb` :
- Simulation de la cascade sans traitement, puis avec blocage de HER2 (type lapatinib) — le signal s'éteint sur les trois étages
- **Mécanisme de résistance** : une mutation activatrice de *PIK3CA* (mécanisme de résistance au trastuzumab documenté cliniquement) maintient PI3K actif indépendamment de HER2 — même avec HER2 bloqué à 95%, AKT reste presque aussi actif que sans traitement
- Simulateur générique pour choisir le point de blocage (HER2, PI3K ou AKT) et comparer au signal non traité
- Courbe dose-réponse : intensité du blocage de HER2 vs niveau d'AKT actif à l'état stationnaire

Constantes cinétiques illustratives (non fittées sur données expérimentales) — l'objectif est de démontrer la méthode de modélisation, pas de produire un résultat biologique quantitatif validé.

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
