from Bio import SeqIO
from io import StringIO
import requests

from src.config import UNIPROT_DIR


def fetch_uniprot_sequence(uniprot_id, save_fasta=False):
    """
    Retrieve a UniProt protein sequence as a BioPython SeqRecord.

    Args:
        uniprot_id (str): UniProt accession ID
        save_fasta (bool): optionally save FASTA to data/uniprot/

    Returns:
        SeqRecord
    """

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    r = requests.get(url)
    r.raise_for_status()

    fasta_text = r.text
    fasta_stream = StringIO(fasta_text)

    record = SeqIO.read(fasta_stream, "fasta")

    if save_fasta:
        out_file = UNIPROT_DIR / f"{uniprot_id}.fasta"
        with open(out_file, "w") as f:
            f.write(fasta_text)

    return record


if __name__ == "__main__":
    seq = fetch_uniprot_sequence("Q1HGV3", save_fasta=True)
    print(seq.id)
    print(len(seq.seq))