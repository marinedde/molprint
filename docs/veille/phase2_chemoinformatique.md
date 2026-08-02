# Veille bibliographique — Phase 2 (chémoinformatique / QSAR)

*Note de méthode : ces fiches pointent vers des sources à consulter, comme le ferait une vraie veille — elles ne remplacent pas la lecture de l'article original avant citation formelle dans un document officiel.*

## Sujet : la règle des 5 de Lipinski (druggabilité)

**Source** :

1. Lipinski, C.A., Lombardo, F., Dominy, B.W., Feeney, P.J. (1997/2001). *Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings*. Advanced Drug Delivery Reviews, 23(1-3), 3–25.

**Synthèse** : à partir d'une analyse rétrospective de médicaments oraux approuvés, Lipinski propose une heuristique simple (poids moléculaire ≤ 500, LogP ≤ 5, ≤ 5 donneurs de liaison H, ≤ 10 accepteurs) prédictive de la biodisponibilité orale. Ce n'est pas une loi physique — c'est une observation statistique avec de nombreuses exceptions connues (médicaments injectables, produits naturels) — d'où son usage dans MolPrint comme indicateur informatif (`lipinski_violations`), jamais comme filtre d'exclusion strict.

## Sujet : les empreintes moléculaires de Morgan / ECFP

**Source** :

2. Rogers, D., Hahn, M. (2010). *Extended-Connectivity Fingerprints*. Journal of Chemical Information and Modeling, 50(5), 742–754.

**Synthèse** : formalise les empreintes circulaires (ECFP, dérivées de l'algorithme de Morgan) comme représentation standard en QSAR et recherche de similarité — chaque bit encode la présence d'un environnement atomique local (rayon fixé, ici 2 liaisons dans MolPrint). Papier fondateur qui justifie le choix du notebook 04 de basculer des descripteurs globaux vers cette représentation, plus fine.

## Sujet : BRICS — fragmentation et recombinaison moléculaire

**Source** :

3. Degen, J., Wegscheid-Gerlach, C., Zaliani, A., Rarey, M. (2008). *On the Art of Compiling and Using 'Drug-Like' Chemical Fragment Spaces*. ChemMedChem, 3(10), 1503–1507.

**Synthèse** : introduit BRICS (*Breaking Retrosynthetically Interesting Chemical Substructures*) — une méthode de fragmentation qui coupe les molécules à des liaisons chimiquement sensées (celles qu'un chimiste romprait pour une synthèse), contrairement à une fragmentation aléatoire. C'est l'algorithme utilisé tel quel (implémentation RDKit) dans `notebooks/03_virtual_screening_optimization.ipynb` pour générer de nouveaux candidats.

## Sujet : squelettes de Bemis-Murcko et split par squelette

**Sources** :

4. Bemis, G.W., Murcko, M.A. (1996). *The Properties of Known Drugs. 1. Molecular Frameworks*. Journal of Medicinal Chemistry, 39(15), 2887–2893.
5. Wu, Z. et al. (2018). *MoleculeNet: A Benchmark for Molecular Machine Learning*. Chemical Science, 9(2), 513–530.

**Synthèse** : Bemis & Murcko définissent la notion de "squelette moléculaire" (l'ossature d'une molécule, décorations périphériques retirées) — la méthode utilisée dans `notebooks/06_scaffold_validation_applicability_domain.ipynb` pour regrouper les molécules par famille structurelle. MoleculeNet (Wu et al.) établit la pratique désormais standard du domaine : évaluer un modèle QSAR avec un split train/test par squelette plutôt qu'aléatoire, précisément parce qu'un split aléatoire surestime la performance sur une chimie réellement nouvelle — le résultat quantifié dans ce même notebook (ROC AUC 0,98 → 0,91).

## Sujet : résistance au trastuzumab et mutations PIK3CA

**Source** :

6. Berns, K. et al. (2007). *A functional genetic approach identifies the PI3K pathway as a major determinant of trastuzumab resistance in breast cancer*. Cancer Cell, 12(4), 395–402.

**Synthèse** : identifie par criblage fonctionnel que l'activation de la voie PI3K (notamment via des mutations de PIK3CA) est un déterminant majeur de résistance au trastuzumab, indépendamment du statut HER2 lui-même. Base scientifique directe du scénario de résistance simulé dans `notebooks/07_systems_biology_her2_pathway.ipynb`.
