# Veille bibliographique — Phase 3 (biologie des systèmes)

*Note de méthode : ces fiches pointent vers des sources à consulter, comme le ferait une vraie veille — elles ne remplacent pas la lecture de l'article original avant citation formelle dans un document officiel.*

## Sujet : la voie HER2/PI3K/AKT dans le cancer du sein

**Source** :

1. Baselga, J., Swain, S.M. (2009). *Novel anticancer targets: revisiting ERBB2 and discovering ERBB3*. Nature Reviews Cancer, 9(7), 463–475.

**Synthèse** : revue de référence sur la biologie de la famille de récepteurs ErbB (dont HER2/ERBB2) et leur signalisation en aval via PI3K/AKT et MAPK. Sert de base biologique générale à la topologie du modèle ODE de `notebooks/07_systems_biology_her2_pathway.ipynb` (cascade HER2 → PI3K → AKT) — un modèle volontairement simplifié par rapport à la complexité réelle décrite dans cette revue (qui inclut aussi ERBB1/3/4, les hétérodimères, et de multiples voies en aval).

## Sujet : modélisation par équations différentielles ordinaires (ODE) en biologie des systèmes

**Source** :

2. Kholodenko, B.N. (2006). *Cell-signalling dynamics in time and space*. Nature Reviews Molecular Cell Biology, 7(3), 165–176.

**Synthèse** : revue méthodologique sur la modélisation dynamique des voies de signalisation cellulaire par systèmes d'équations différentielles — la classe de modèle utilisée par Tellurium. Aborde notamment pourquoi la cinétique (pas seulement la présence/absence d'une interaction) est nécessaire pour comprendre des phénomènes comme les boucles de rétroaction, absentes du modèle simplifié de MolPrint (qui n'a volontairement aucune rétroaction, une simplification explicitement signalée dans le notebook).

## Sujet : outils logiciels — Tellurium / Antimony / SBML

**Source** :

3. Choi, K. et al. (2018). *Tellurium: An extensible python-based modeling environment for systems and synthetic biology*. Biosystems, 171, 74–79.

**Synthèse** : article de présentation de l'outil utilisé directement dans MolPrint. Antimony (le langage texte utilisé pour écrire le modèle) et SBML (*Systems Biology Markup Language*, le format d'échange standard sous-jacent) sont les briques qui permettent d'écrire un modèle de réseau de réactions de façon lisible plutôt qu'en équations différentielles brutes.
