---
title: MolPrint Dashboard
emoji: 🧪
colorFrom: purple
colorTo: blue
sdk: streamlit
sdk_version: "1.38.0"
app_file: streamlit_app.py
pinned: false
license: mit
---

# 🧪 MolPrint — Dashboard

**Dashboard Streamlit — Marine Deldicque**

Criblage virtuel et prédiction de réponse aux médicaments in silico, cible ERBB2/HER2. Complément d'[OncoPrint](https://huggingface.co/spaces/marinedde/oncoprint-dashboard) côté chémoinformatique.

## Pages

| Page | Description |
|------|-------------|
| 🏠 Accueil | Vue d'ensemble et métriques clés |
| 🧪 Cribler une molécule | Prédiction d'activité à partir d'un SMILES |
| 🎯 Candidats sélectionnés | Top 15 candidats du criblage virtuel (BRICS) |
| 🧬 Séquences géniques | Analyse Biopython des 5 gènes clés d'OncoPrint |
| 📊 Méthodologie & rigueur | Validation par squelette chimique, domaine d'applicabilité |
| ℹ️ À propos | Contexte et stack technique |

## Stack

`Python` · `RDKit` · `XGBoost` · `Streamlit` · `HuggingFace Spaces`

## Code source

[github.com/marinedde/molprint](https://github.com/marinedde/molprint)
