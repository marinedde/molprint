import io
import sqlite3

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Draw, Lipinski, rdFingerprintGenerator

st.set_page_config(page_title="MolPrint", page_icon="🧪", layout="wide")

DESCRIPTOR_COLS = [
    "molecular_weight", "logp", "h_bond_donors", "h_bond_acceptors",
    "rotatable_bonds", "tpsa", "aromatic_rings", "ring_count",
    "heavy_atoms", "molar_refractivity",
]
MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


@st.cache_resource
def load_model():
    return joblib.load("models/qsar_erbb2_combined_xgboost.pkl")


@st.cache_data
def load_training_data():
    return pd.read_csv("data/erbb2_activities.csv")


@st.cache_resource
def load_training_fingerprints(_train_df):
    return [MORGAN_GEN.GetFingerprint(Chem.MolFromSmiles(s)) for s in _train_df["canonical_smiles"]]


@st.cache_data
def load_shortlist():
    return pd.read_csv("data/erbb2_candidate_shortlist_with_domain.csv")


@st.cache_data
def load_gene_table():
    return pd.read_csv("data/gene_sequence_analysis.csv")


@st.cache_resource
def get_db_connection():
    return sqlite3.connect("data/molprint.db", check_same_thread=False)


def compute_descriptors(mol):
    return [
        Descriptors.MolWt(mol), Descriptors.MolLogP(mol), Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol), Descriptors.NumRotatableBonds(mol), Descriptors.TPSA(mol),
        Descriptors.NumAromaticRings(mol), Descriptors.RingCount(mol), Descriptors.HeavyAtomCount(mol),
        Descriptors.MolMR(mol),
    ]


def lipinski_violations(desc):
    mw, logp, hbd, hba = desc[0], desc[1], desc[2], desc[3]
    return int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)


def applicability_domain(similarity):
    if similarity >= 0.6:
        return "proche du train (fiable)", "🟢"
    if similarity >= 0.35:
        return "modérément nouveau", "🟡"
    return "extrapolation (prudence)", "🔴"


def screen_smiles(smiles, model, train_fps):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    desc = compute_descriptors(mol)
    fp = MORGAN_GEN.GetFingerprint(mol)
    features = np.array(desc + list(fp)).reshape(1, -1)
    proba = float(model.predict_proba(features)[0, 1])
    similarity = max(DataStructs.BulkTanimotoSimilarity(fp, train_fps))
    domain_label, domain_icon = applicability_domain(similarity)
    return {
        "mol": mol,
        "activity_probability": proba,
        "lipinski_violations": lipinski_violations(desc),
        "molecular_weight": desc[0],
        "logp": desc[1],
        "nn_similarity": similarity,
        "domain_label": domain_label,
        "domain_icon": domain_icon,
    }


def mol_image(mol, size=(300, 300)):
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


PAGES = [
    "🏠 Accueil",
    "🧪 Cribler une molécule",
    "🎯 Candidats sélectionnés",
    "🧬 Séquences géniques",
    "🗄️ Base de données",
    "📊 Méthodologie & rigueur",
    "ℹ️ À propos",
]

st.sidebar.title("🧪 MolPrint")
page = st.sidebar.radio("Navigation", PAGES, label_visibility="collapsed")

model = load_model()
train_df = load_training_data()
train_fps = load_training_fingerprints(train_df)

if page == "🏠 Accueil":
    st.title("🧪 MolPrint")
    st.caption("Prédicteur de réponse aux médicaments in silico — cible ERBB2/HER2")
    st.markdown(
        "Complément d'[OncoPrint](https://huggingface.co/spaces/marinedde/oncoprint-dashboard) "
        "(classification moléculaire du cancer du sein) côté chémoinformatique : criblage virtuel "
        "et prédiction d'activité de molécules sur la cible HER2."
    )

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Molécules", len(train_df), help="Nombre de molécules dans le jeu d'entraînement")
    col2.metric("AUC aléatoire", "0.98", help="ROC AUC — split train/test aléatoire (optimiste, voir Méthodologie)")
    col3.metric("AUC squelette", "0.91", help="ROC AUC — split par squelette chimique, estimation honnête sur des squelettes inédits (voir l'onglet Méthodologie)")
    col4.metric("Candidats", len(load_shortlist()), help="Candidats retenus après criblage virtuel")

    st.markdown("---")
    st.markdown(
        "**Comment lire ce dashboard** : commence par *Cribler une molécule* pour tester "
        "n'importe quel SMILES, regarde les *Candidats sélectionnés* issus de l'optimisation "
        "in silico, ou va directement à *Méthodologie* pour voir les limites honnêtes du modèle."
    )

elif page == "🧪 Cribler une molécule":
    st.title("🧪 Cribler une molécule")
    st.markdown("Entre un SMILES pour prédire son activité sur la cible ERBB2/HER2.")

    examples = {
        "Lapatinib (approuvé)": "CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
        "Tucatinib (approuvé)": "CC1=C(C=CC(=C1)NC2=NC=NC3=C2C=C(C=C3)NC4=NC(CO4)(C)C)OC5=CC6=NC=NN6C=C5",
        "Aspirine (négatif attendu)": "CC(=O)OC1=CC=CC=C1C(=O)O",
    }
    choice = st.selectbox("Exemple rapide (optionnel)", ["—"] + list(examples.keys()))
    default_smiles = examples.get(choice, "")
    smiles = st.text_input("SMILES", value=default_smiles, placeholder="ex : CC(=O)OC1=CC=CC=C1C(=O)O")

    if smiles:
        result = screen_smiles(smiles, model, train_fps)
        if result is None:
            st.error("SMILES invalide — impossible de le parser avec RDKit.")
        else:
            col1, col2 = st.columns([1, 2])
            with col1:
                st.image(mol_image(result["mol"]), caption=smiles)
            with col2:
                st.metric("Probabilité d'activité (ERBB2/HER2)", f"{result['activity_probability']:.1%}")
                st.metric(
                    "Domaine d'applicabilité",
                    f"{result['domain_icon']} {result['domain_label']}",
                    help=f"Similarité de Tanimoto à la molécule la plus proche du jeu d'entraînement : {result['nn_similarity']:.2f}",
                )
                sub1, sub2, sub3 = st.columns(3)
                sub1.metric("Poids moléculaire", f"{result['molecular_weight']:.0f} g/mol")
                sub2.metric("LogP", f"{result['logp']:.2f}")
                sub3.metric("Violations Lipinski", result["lipinski_violations"])
                if result["nn_similarity"] < 0.35:
                    st.warning(
                        "Ce candidat est chimiquement éloigné du jeu d'entraînement : le score de "
                        "probabilité doit être pris avec prudence (voir onglet Méthodologie)."
                    )

elif page == "🎯 Candidats sélectionnés":
    st.title("🎯 Candidats sélectionnés par criblage virtuel")
    st.markdown(
        "Top 15 candidats générés par recombinaison de fragments (BRICS) à partir des molécules "
        "les plus actives, filtrés par diversité structurale (clustering Butina) — "
        "voir `notebooks/03_virtual_screening_optimization.ipynb`."
    )
    shortlist = load_shortlist()
    display_cols = [
        "canonical_smiles", "activity_probability", "nn_similarity_to_train",
        "applicability_domain", "molecular_weight", "logp", "lipinski_violations",
    ]
    st.dataframe(
        shortlist[display_cols].rename(columns={
            "canonical_smiles": "SMILES",
            "activity_probability": "Probabilité d'activité",
            "nn_similarity_to_train": "Similarité au train",
            "applicability_domain": "Domaine d'applicabilité",
            "molecular_weight": "PM (g/mol)",
            "logp": "LogP",
            "lipinski_violations": "Violations Lipinski",
        }),
        use_container_width=True,
    )

    st.markdown("### Structures des 6 meilleurs candidats")
    top6 = shortlist.head(6)
    cols = st.columns(3)
    for i, row in enumerate(top6.itertuples()):
        mol = Chem.MolFromSmiles(row.canonical_smiles)
        with cols[i % 3]:
            st.image(mol_image(mol, size=(220, 220)))
            st.caption(f"p={row.activity_probability:.2f} · {row.applicability_domain}")

elif page == "🧬 Séquences géniques":
    st.title("🧬 Séquences géniques (Phase 1 — Biopython)")
    st.markdown(
        "Les 5 gènes les plus importants en SHAP dans OncoPrint, un par sous-type de cancer du "
        "sein : **ESR1, ERBB2, FOXA1, AR, GATA3**. Séquences RefSeq (NCBI) récupérées et "
        "traduites via `src/sequence_analysis.py`."
    )
    genes = load_gene_table()
    st.dataframe(
        genes.rename(columns={
            "gene": "Gène", "accession": "Accession RefSeq", "mrna_length": "Longueur ARNm (nt)",
            "mrna_gc_percent": "GC%", "protein_length_aa": "Longueur protéine (aa)",
            "protein_mw_kda": "PM protéine (kDa)", "instability_index": "Indice d'instabilité",
            "isoelectric_point": "Point isoélectrique",
        }).drop(columns=["description", "cds_length_nt", "aromaticity"], errors="ignore"),
        use_container_width=True,
    )
    st.info(
        "**Fait notable** : la position 537 de la protéine ESR1 traduite depuis ce RefSeq est une "
        "tyrosine (Y) — le site exact de la mutation de résistance à l'hormonothérapie **Y537S**, "
        "l'une des mieux documentées dans le cancer du sein métastatique. "
        "Détail dans `notebooks/05_gene_sequence_analysis.ipynb`."
    )

elif page == "🗄️ Base de données":
    st.title("🗄️ Base de données interne (Phase 4)")
    st.markdown(
        "Relie gènes (Phase 1), molécules candidates (Phase 2), voie de signalisation (Phase 3) "
        "et sous-types moléculaires — enrichie par **GDSC** (Sanger), des données réelles de "
        "sensibilité de 51 lignées cellulaires de cancer du sein à des thérapies ciblées. "
        "Voir `notebooks/08_internal_database.ipynb`."
    )
    conn = get_db_connection()

    st.markdown("### Explorer par sous-type")
    subtype_choice = st.selectbox("Sous-type moléculaire", ["Luminal A", "HER2-enrichi", "Triple Negatif"])

    key_genes_df = pd.read_sql(
        "SELECT key_gene, shap_rank FROM subtypes WHERE subtype = ? ORDER BY shap_rank",
        conn, params=(subtype_choice,),
    )
    key_genes = key_genes_df["key_gene"].tolist()

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Gènes clés (SHAP OncoPrint)**")
        st.dataframe(
            key_genes_df.rename(columns={"key_gene": "Gène", "shap_rank": "Rang SHAP"}),
            use_container_width=True, hide_index=True,
        )
    with col2:
        n_lines = pd.read_sql(
            "SELECT COUNT(*) as n FROM cell_lines WHERE subtype = ?", conn, params=(subtype_choice,)
        )["n"][0]
        st.metric("Lignées cellulaires GDSC de ce sous-type", n_lines)

    if "ERBB2" in key_genes:
        st.markdown("**Molécules candidates ciblant ERBB2** (Phase 2)")
        mol_q = pd.read_sql(
            "SELECT canonical_smiles, activity_probability, applicability_domain "
            "FROM molecules WHERE target_gene = 'ERBB2' ORDER BY activity_probability DESC LIMIT 5",
            conn,
        )
        st.dataframe(mol_q, use_container_width=True, hide_index=True)
        st.markdown("**Voie de signalisation modélisée (Phase 3)** : HER2 → PI3K → AKT")
    else:
        st.info(
            "Aucune molécule candidate ni voie modélisée pour ce sous-type dans MolPrint pour "
            "l'instant — seul ERBB2/HER2 a été couvert en Phase 2/3."
        )

    st.markdown("### Validation empirique GDSC : sensibilité au lapatinib par sous-type")
    lap_df = pd.read_sql(
        '''
        SELECT c.subtype, r.LN_IC50
        FROM drug_response r
        JOIN cell_lines c ON r.COSMIC_ID = c.COSMIC_ID
        WHERE r.DRUG_NAME = 'Lapatinib' AND c.subtype != 'Non classé'
        ''',
        conn,
    )
    fig, ax = plt.subplots(figsize=(6, 3.5))
    lap_df.boxplot(column="LN_IC50", by="subtype", ax=ax)
    ax.set_ylabel("ln(IC50) — plus bas = plus sensible")
    ax.set_title("")
    plt.suptitle("")
    st.pyplot(fig)
    st.caption("Données réelles de laboratoire (GDSC/Sanger) — indépendantes de nos propres modèles.")

    st.markdown("### Requêter la base toi-même (SQL)")
    st.caption("Lecture seule — uniquement des requêtes SELECT.")
    example_queries = {
        "Tables disponibles": "SELECT name FROM sqlite_master WHERE type='table'",
        "Tous les gènes": "SELECT * FROM genes",
        "Médicaments testés sur les lignées HER2-enrichi (triés par efficacité)": (
            "SELECT r.DRUG_NAME, AVG(r.LN_IC50) as ic50_moyen, COUNT(*) as n "
            "FROM drug_response r JOIN cell_lines c ON r.COSMIC_ID = c.COSMIC_ID "
            "WHERE c.subtype = 'HER2-enrichi' GROUP BY r.DRUG_NAME ORDER BY ic50_moyen"
        ),
    }
    choice_q = st.selectbox("Exemple de requête (optionnel)", ["—"] + list(example_queries.keys()))
    default_query = example_queries.get(choice_q, "SELECT * FROM genes LIMIT 5")
    query = st.text_area("Requête SQL", value=default_query, height=100)

    if st.button("Exécuter"):
        if not query.strip().lower().startswith("select"):
            st.error("Seules les requêtes SELECT sont autorisées.")
        else:
            try:
                st.dataframe(pd.read_sql(query, conn), use_container_width=True)
            except Exception as e:
                st.error(f"Erreur SQL : {e}")

elif page == "📊 Méthodologie & rigueur":
    st.title("📊 Méthodologie & rigueur")
    st.markdown(
        "Un split train/test aléatoire est optimiste en QSAR : le jeu de données contient "
        "beaucoup d'analogues proches (**499 squelettes chimiques uniques pour 1070 molécules**), "
        "donc des quasi-doublons peuvent se retrouver de part et d'autre du split."
    )

    col1, col2 = st.columns(2)
    with col1:
        st.metric("ROC AUC — split aléatoire", "0.981", help="Notebook 04 — optimiste")
    with col2:
        st.metric("ROC AUC — split par squelette", "0.912", help="Notebook 06 — squelettes jamais vus en entraînement, estimation honnête")

    st.markdown(
        "Sur des squelettes chimiques réellement inédits, le rappel sur la classe active tombe "
        "à **0.59** : le modèle est fiable pour prioriser des analogues de séries chimiques "
        "connues, beaucoup moins pour explorer une chimie franchement nouvelle.\n\n"
        "C'est pour ça que chaque prédiction de l'onglet *Cribler une molécule* est accompagnée "
        "d'un **domaine d'applicabilité** (similarité de Tanimoto à la molécule la plus proche du "
        "jeu d'entraînement) plutôt que d'un score brut présenté seul."
    )
    st.markdown("Détail complet : `notebooks/06_scaffold_validation_applicability_domain.ipynb`.")

    st.markdown("### \"0,98, ce n'est pas trop beau ?\"")
    st.markdown(
        "Test de permutation : réentraîner le modèle avec les labels actif/inactif mélangés "
        "au hasard fait retomber le ROC AUC à **~0,5** (0,567 / 0,504 / 0,469 sur 3 essais) — "
        "pas de fuite de données cachée.\n\n"
        "Mais le bit d'empreinte le plus important (35% du poids du modèle) correspond au motif "
        "**aminopyrimidine/quinazoline**, le point d'ancrage classique des inhibiteurs de kinases "
        "— présent dans 86% des actifs contre 5% des inactifs. Le modèle discrimine surtout "
        "*\"type inhibiteur de kinase\"* vs *\"non\"*, pas finement la sélectivité pour HER2 "
        "spécifiquement. Détail dans `docs/GUIDE_COMPLET.md`."
    )

elif page == "ℹ️ À propos":
    st.title("ℹ️ À propos de MolPrint")
    st.markdown(
        """
Projet in silico de prédiction de réponse aux médicaments, construit en complément d'
[OncoPrint](https://github.com/marinedde/cdsd-certification/tree/main/bloc6-direction-projet/oncoprint)
pour couvrir les compétences bioinformatique / chémoinformatique demandées sur les postes
Data Scientist Préclinique en R&D pharmaceutique.

**Stack** : `Python` · `RDKit` · `Biopython` · `XGBoost` · `SHAP` · `scikit-learn` · `Streamlit` · `HuggingFace Spaces`

**Code source** : [github.com/marinedde/molprint](https://github.com/marinedde/molprint)

**Roadmap** :
- ✅ Phase 1 — Séquences ARN/ADN (Biopython)
- ✅ Phase 2 — Chémoinformatique (RDKit, QSAR, criblage virtuel, validation par squelette)
- ✅ Phase 3 — Biologie des systèmes (modèle ODE HER2/PI3K/AKT, Tellurium)
- ✅ Phase 4 — Base de données interne (SQLite + GDSC)
        """
    )
