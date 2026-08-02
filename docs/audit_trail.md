# Piste d'audit — décisions et exclusions

Registre des décisions de traitement de données et de modélisation prises dans MolPrint, avec le raisonnement et l'impact quantifié quand c'est possible. Objectif : dans six mois, pouvoir répondre à "pourquoi cette ligne a disparu / pourquoi ce choix" sans avoir à reconstituer le fil depuis zéro.

Format : Décision — Raisonnement — Impact mesuré — Notebook/commit source.

---

### Source de données Phase 2 : PubChem plutôt que ChEMBL
**Décision** : basculer de ChEMBL (visé initialement) vers PubChem BioAssay.
**Raisonnement** : l'API ChEMBL retournait des erreurs 500 au moment de l'écriture — vérifié directement en `curl`, pas supposé.
**Impact** : aucune perte de rigueur — PubChem est la deuxième référence publique du domaine.
**Source** : `notebooks/01_bioactivity_data_acquisition.ipynb`, intro.

### Agrégation des essais multiples par molécule
**Décision** : quand une molécule (CID) apparaît dans plusieurs essais biologiques, garder le résultat majoritaire (`Activity Outcome`) et la meilleure valeur observée (IC50 le plus bas).
**Raisonnement** : réduire chaque molécule à une seule ligne, condition nécessaire pour éviter des doublons entre train et test.
**Impact** : 1178 lignes brutes → 1070 molécules uniques.
**Source** : `notebooks/01_bioactivity_data_acquisition.ipynb`, section 3.

### Bug `dropna()` trop agressif (notebook 02)
**Décision** : ne faire `dropna()` que sur les colonnes de descripteurs, pas sur `pactivity` (absente pour les essais qualitatifs sans valeur chiffrée).
**Raisonnement** : la première version droppait silencieusement toutes les lignes sans valeur IC50 chiffrée, biaisant le jeu vers ~95% actif.
**Impact mesuré** : jeu de test passé d'un déséquilibre extrême à un split équilibré (105/109) après correction.
**Source** : `notebooks/02_rdkit_descriptors_qsar.ipynb`, commentaire dans la cellule de dropna ; trouvé et corrigé en cours de session.

### Choix du modèle final : combiné plutôt que descripteurs seuls
**Décision** : le modèle de production est celui entraîné sur descripteurs + empreintes de Morgan, pas sur descripteurs seuls.
**Raisonnement** : comparé en validation croisée à 5 plis, pas au jugé.
**Impact mesuré** : ROC AUC 0,9545 (descripteurs) → 0,9747 (empreintes) → 0,9803 (combiné).
**Source** : `notebooks/04_fingerprint_qsar.ipynb`.

### Bug de cohérence : notebook 03 chargeait l'ancien modèle
**Décision** : corriger `notebooks/03_virtual_screening_optimization.ipynb` pour qu'il charge le modèle combiné (notebook 04), pas le modèle descripteurs-seuls (notebook 02).
**Raisonnement** : trouvé en relisant le code pendant l'audit data leakage/overfitting — le criblage virtuel, le domaine d'applicabilité (notebook 06) et la table `molecules` de la base (Phase 4) héritaient tous du modèle le moins rigoureux.
**Impact mesuré** : le shortlist de candidats a changé après correction (nouveaux SMILES, nouveaux scores) ; notebooks 03 → 06 → 08 rejoués dans l'ordre, dashboard et base redéployés.
**Source** : commit "Corriger un bug de cohérence : notebook 03 utilisait l'ancien modèle".

### Molécule outlier (CID 136027061) : conservée, pas exclue
**Décision** : ne pas retirer cette molécule (14 donneurs H, TPSA 355, très atypique) du jeu de données.
**Raisonnement** : c'est une vraie molécule, pas une erreur de saisie ; elle est étiquetée inactive donc ne biaise pas la classe active.
**Impact mesuré** : test de sensibilité — ROC AUC en validation croisée avec (0,9803) vs sans (0,9799) l'outlier — différence dans le bruit statistique, retrait non justifié.
**Source** : `notebooks/09_exploratory_data_analysis.ipynb`, section 4 ; vérification complémentaire faite en session.

### Split par squelette plutôt qu'aléatoire pour évaluer la généralisation
**Décision** : citer le ROC AUC par squelette (0,91) comme chiffre honnête de généralisation, pas le ROC AUC aléatoire (0,98).
**Raisonnement** : 1070 molécules pour seulement 499 squelettes uniques (un squelette a 66 dérivés) — un split aléatoire laisse fuiter des quasi-doublons entre train et test.
**Impact mesuré** : ROC AUC hold-out 0,9812 (aléatoire) → 0,9119 (squelette), rappel classe active 0,59 sur squelettes inédits.
**Source** : `notebooks/06_scaffold_validation_applicability_domain.ipynb`.

### Seuils du domaine d'applicabilité
**Décision** : similarité de Tanimoto ≥ 0,6 = "fiable", 0,35-0,6 = "modérément nouveau", < 0,35 = "extrapolation (prudence)".
**Raisonnement** : seuils usuels du domaine QSAR (pas de calibration statistique formelle sur ce projet — limite assumée).
**Impact** : 1 des 15 candidats du shortlist final classé "extrapolation", le reste "modérément nouveau" — aucun n'est "proche du train" malgré des scores élevés, d'où l'importance de l'avertissement dans le dashboard.
**Source** : `notebooks/06_scaffold_validation_applicability_domain.ipynb`, section 3 ; `dashboard/streamlit_app.py`.

### Classification des lignées cellulaires GDSC par sous-type : volontairement incomplète
**Décision** : classer seulement 33 des 51 lignées de cancer du sein GDSC par sous-type moléculaire ; les 18 restantes marquées "Non classé".
**Raisonnement** : classification manuelle basée sur la littérature (Neve et al. 2006, Kao et al. 2009) — seules les lignées avec un consensus suffisant dans ces sources sont classées. Deviner aurait introduit une erreur silencieuse.
**Impact** : les 18 lignées non classées sont exclues des analyses par sous-type (ex. comparaison de sensibilité au lapatinib), mais restent dans la base pour d'autres usages.
**Source** : `notebooks/08_internal_database.ipynb`, section 2 ; `docs/veille/phase4_base_donnees.md`.

### Pas de rééquilibrage de classe (SMOTE ou autre)
**Décision** : aucune technique de rééquilibrage appliquée sur le jeu de données ERBB2.
**Raisonnement** : le jeu est déjà équilibré (51%/49%, vérifié en EDA) — un rééquilibrage artificiel n'aurait fait qu'ajouter du risque (ex. la fuite SMOTE-avant-split rencontrée sur un autre projet, OncoPrint) sans bénéfice.
**Impact** : aucun — décision de ne rien faire, documentée pour que ce ne soit pas relu comme un oubli.
**Source** : `notebooks/09_exploratory_data_analysis.ipynb`, section 1 ; `src/data_quality.py` (le contrôle qualité alerte si le déséquilibre dépassait 25%/75%).

### Filtre de poids moléculaire sur les candidats générés (BRICS)
**Décision** : exclure les candidats générés par recombinaison BRICS dont le poids moléculaire dépasse 700 g/mol.
**Raisonnement** : la recombinaison combinatoire peut chaîner des fragments indéfiniment et produire des molécules chimiquement peu plausibles comme candidats-médicaments.
**Impact** : réduit le nombre de candidats bruts avant filtrage de diversité (voir notebook 03 pour les volumes exacts à chaque étape).
**Source** : `notebooks/03_virtual_screening_optimization.ipynb`, section 2.

---

## Contrôles automatiques en place

Depuis l'ajout de `src/data_quality.py`, `notebooks/01_bioactivity_data_acquisition.ipynb` vérifie automatiquement, avant de sauvegarder le fichier que tous les autres notebooks consomment :
- présence des colonnes requises
- nombre de lignes dans une fourchette attendue (500–5000)
- pas de SMILES ou d'étiquette `active` manquants
- pas de CID ni de SMILES dupliqués
- tous les SMILES parsables par RDKit
- déséquilibre de classe dans une fourchette raisonnable (25%–75% actif)

Toute violation lève `DataQualityError` et **arrête le notebook** plutôt que de sauvegarder un fichier douteux en silence.
