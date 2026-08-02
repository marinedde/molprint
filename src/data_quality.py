"""Assertions de qualité sur le jeu de données ERBB2.

Principe : alerter (lever une exception claire) plutôt que corriger en
silence. Utilisé dans notebooks/01_bioactivity_data_acquisition.ipynb
juste avant de sauvegarder le CSV final, et rejouable indépendamment
sur n'importe quelle version du fichier.
"""
import pandas as pd
from rdkit import Chem


class DataQualityError(Exception):
    """Levée quand une assertion de qualité échoue — le pipeline doit s'arrêter, pas continuer silencieusement."""


def check_erbb2_activities(df, min_rows=500, max_rows=5000, min_active_share=0.25, max_active_share=0.75):
    """Vérifie le jeu de données ERBB2 nettoyé. Lève DataQualityError au premier problème trouvé.

    Retourne un dict de statistiques si tout passe (pour log/affichage).
    """
    errors = []

    required_cols = {"cid", "canonical_smiles", "active"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        errors.append(f"Colonnes manquantes : {missing_cols}")
        raise DataQualityError("; ".join(errors))

    n = len(df)
    if not (min_rows <= n <= max_rows):
        errors.append(f"Nombre de lignes hors bornes attendues [{min_rows}, {max_rows}] : {n}")

    n_null_smiles = df["canonical_smiles"].isnull().sum()
    if n_null_smiles > 0:
        errors.append(f"{n_null_smiles} lignes avec canonical_smiles manquant")

    n_null_active = df["active"].isnull().sum()
    if n_null_active > 0:
        errors.append(f"{n_null_active} lignes avec 'active' manquant")

    n_dup_cid = df["cid"].duplicated().sum()
    if n_dup_cid > 0:
        errors.append(f"{n_dup_cid} CID dupliqués (une molécule ne doit apparaître qu'une fois)")

    n_dup_smiles = df["canonical_smiles"].duplicated().sum()
    if n_dup_smiles > 0:
        errors.append(f"{n_dup_smiles} canonical_smiles dupliqués (risque de fuite train/test)")

    n_invalid_smiles = sum(1 for s in df["canonical_smiles"] if Chem.MolFromSmiles(s) is None)
    if n_invalid_smiles > 0:
        errors.append(f"{n_invalid_smiles} SMILES non parsables par RDKit")

    active_share = df["active"].mean()
    if not (min_active_share <= active_share <= max_active_share):
        errors.append(
            f"Déséquilibre de classe hors bornes attendues [{min_active_share:.0%}, {max_active_share:.0%}] : "
            f"{active_share:.1%} actif"
        )

    if errors:
        raise DataQualityError("Contrôle qualité échoué :\n- " + "\n- ".join(errors))

    return {
        "n_rows": n,
        "active_share": round(active_share, 4),
        "n_unique_cid": df["cid"].nunique(),
        "n_unique_smiles": df["canonical_smiles"].nunique(),
    }
