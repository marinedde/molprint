# Veille bibliographique — Phase 4 (base de données interne / GDSC)

*Note de méthode : ces fiches pointent vers des sources à consulter, comme le ferait une vraie veille — elles ne remplacent pas la lecture de l'article original avant citation formelle dans un document officiel.*

## Sujet : GDSC — la base de sensibilité aux médicaments

**Sources** :

1. Yang, W. et al. (2013). *Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells*. Nucleic Acids Research, 41(D1), D955–D961.
2. Iorio, F. et al. (2016). *A Landscape of Pharmacogenomic Interactions in Cancer*. Cell, 166(3), 740–754.

**Synthèse** : ces deux publications documentent la construction et la méthodologie de GDSC (Sanger Institute) — le criblage systématique de centaines de lignées cellulaires de cancer contre des centaines de composés, avec courbes dose-réponse ajustées (IC50, AUC). Iorio et al. (2016) relient en plus ces réponses aux profils génomiques des lignées, posant les bases de la logique "biomarqueur génomique → réponse au traitement" que `notebooks/08_internal_database.ipynb` reproduit à petite échelle (sous-type moléculaire → réponse au lapatinib).

## Sujet : classification moléculaire des lignées cellulaires de cancer du sein

**Sources** :

3. Neve, R.M. et al. (2006). *A collection of breast cancer cell lines for the study of functionally distinct cancer subtypes*. Cancer Cell, 10(6), 515–527.
4. Kao, J. et al. (2009). *Molecular Profiling of Breast Cancer Cell Lines Defines Relevant Tumor Models and Provides a Resource for Cancer Gene Discovery*. PLoS ONE, 4(7), e6146.

**Synthèse** : ces deux études caractérisent, par profil d'expression génique, un large panel de lignées cellulaires de cancer du sein et les classent selon les mêmes catégories moléculaires que celles utilisées en clinique (Luminal, HER2-enrichi, Basal-like/Triple Négatif). C'est la source directe de la table de correspondance manuelle `SUBTYPE_MAP` codée dans `notebooks/08_internal_database.ipynb` — limitée aux lignées où ces deux études (et la littérature qui s'en inspire) offrent un consensus suffisant.

## Sujet : sensibilité au lapatinib et biomarqueurs HER2

**Source** :

5. Konecny, G.E. et al. (2006). *Activity of the dual kinase inhibitor lapatinib (GW572016) against HER-2-overexpressing and trastuzumab-treated breast cancer cells*. Cancer Research, 66(3), 1630–1639.

**Synthèse** : caractérise en laboratoire la sensibilité différentielle des lignées de cancer du sein au lapatinib selon leur statut HER2 — cohérent avec, et antérieur à, l'observation empirique reproduite dans MolPrint à partir des données GDSC (lapatinib ~10x plus puissant sur les lignées HER2-enrichi que sur les autres sous-types).
