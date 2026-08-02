# Veille bibliographique — Phase 1 (séquences ARN/ADN)

*Note de méthode : ces fiches pointent vers des sources à consulter, comme le ferait une vraie veille — elles ne remplacent pas la lecture de l'article original avant citation formelle dans un document officiel.*

## Sujet : la mutation ESR1 Y537S et la résistance à l'hormonothérapie

**Sources** :

1. Toy, W. et al. (2013). *ESR1 ligand-binding domain mutations in hormone-resistant breast cancer*. Nature Genetics, 45(12), 1439–1445.
2. Robinson, D.R. et al. (2013). *Activating ESR1 mutations in hormone-resistant metastatic breast cancer*. Nature Genetics, 45(12), 1446–1451.

**Synthèse** : ces deux études, publiées quasi simultanément fin 2013, ont établi que des mutations récurrentes dans le domaine de liaison au ligand (LBD) du récepteur aux œstrogènes ESR1 — dont Y537S — apparaissent spécifiquement dans les tumeurs métastatiques résistantes à l'hormonothérapie, alors qu'elles sont quasi absentes des tumeurs primaires non traitées. Le mécanisme : ces mutations stabilisent une conformation active du récepteur *indépendamment* de la présence d'œstrogène, rendant les traitements anti-œstrogéniques classiques (inhibiteurs de l'aromatase notamment) inefficaces. C'est le fondement clinique direct du résultat de `notebooks/05_gene_sequence_analysis.ipynb` : vérifier que la position 537 de la séquence RefSeq porte bien une tyrosine avant mutation.

**Lien avec MolPrint** : ce résultat motive directement l'idée (non implémentée, piste d'évolution) d'étendre la Phase 2 à ESR1 comme deuxième cible thérapeutique, avec les mutants de résistance comme cas d'usage naturel pour tester un futur modèle de prédiction d'impact de mutation.

## Sujet : NCBI RefSeq comme référentiel de séquences

**Source** :

3. O'Leary, N.A. et al. (2016). *Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation*. Nucleic Acids Research, 44(D1), D733–D745.

**Synthèse** : décrit la méthodologie de curation de RefSeq — la distinction entre séquences "prédites" (préfixe `XM_`/`XP_`, générées automatiquement) et séquences "curées" (préfixe `NM_`/`NP_`, validées manuellement par des biocurateurs). C'est la raison pour laquelle `src/sequence_analysis.py` filtre explicitement sur les accessions `NM_` plutôt que d'accepter le premier résultat de recherche.
