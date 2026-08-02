"""Récupération et analyse de séquences géniques (NCBI RefSeq) via Biopython.

Utilisé par notebooks/01_gene_sequence_analysis.ipynb pour les gènes
les plus importants en SHAP dans OncoPrint (ESR1, ERBB2, FOXA1, AR, GATA3).
"""
import os
import time

import certifi
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

os.environ.setdefault("SSL_CERT_FILE", certifi.where())

Entrez.email = "marine.deldicque@gmail.com"


def fetch_canonical_mrna_accession(gene_symbol, organism="Homo sapiens"):
    """Renvoie le numéro d'accession RefSeq curé (NM_...) le plus pertinent pour ce gène."""
    term = f"{gene_symbol}[Gene Name] AND {organism}[Organism] AND biomol_mrna[PROP] AND RefSeq[Filter]"
    ids = Entrez.read(Entrez.esearch(db="nucleotide", term=term, retmax=50))["IdList"]
    if not ids:
        return None
    summaries = Entrez.read(Entrez.esummary(db="nucleotide", id=",".join(ids)))
    nm_hits = sorted(s["Caption"] for s in summaries if s["Caption"].startswith("NM_"))
    return nm_hits[0] if nm_hits else None


def fetch_gene_record(accession):
    """Télécharge l'enregistrement GenBank complet (séquence + annotations) pour une accession."""
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def gc_content(seq):
    seq = str(seq).upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def analyze_gene(gene_symbol, pause=0.4):
    """Pipeline complet : accession -> séquence -> traduction -> descripteurs protéiques."""
    accession = fetch_canonical_mrna_accession(gene_symbol)
    if accession is None:
        return {"gene": gene_symbol, "error": "aucune séquence RefSeq trouvée"}
    time.sleep(pause)

    record = fetch_gene_record(accession)
    time.sleep(pause)

    cds_features = [f for f in record.features if f.type == "CDS"]
    if not cds_features:
        return {"gene": gene_symbol, "accession": accession, "error": "pas de CDS annotée"}
    cds = cds_features[0]
    cds_seq = cds.extract(record.seq)
    protein_seq = str(cds_seq.translate(to_stop=True))

    analysis = ProteinAnalysis(protein_seq)

    return {
        "gene": gene_symbol,
        "accession": accession,
        "description": record.description,
        "mrna_length": len(record.seq),
        "mrna_gc_percent": round(gc_content(record.seq), 2),
        "cds_length_nt": len(cds_seq),
        "protein_length_aa": len(protein_seq),
        "protein_mw_kda": round(analysis.molecular_weight() / 1000, 2),
        "instability_index": round(analysis.instability_index(), 2),
        "aromaticity": round(analysis.aromaticity(), 3),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "protein_sequence": protein_seq,
    }
