from Bio import SeqIO
from io import StringIO
import requests


def fetch_uniprot_sequence(uniprot_id):
    """
    CREATE FASTA OBJECT FROM UNIPROT

    Args:
        uniprot_id (str): Uniprot Identifier
    Returns: 
        SeqRecord object & metadata dic
    """
    #API CALL TO UNIPROT
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url)
    r.raise_for_status() #flag unsuccessful responses as errors
    fasta = r.text

    #ITERATE FASTA RESPONSE & SAVE
    fasta_io = StringIO(fasta) # create in memory file object
    record = SeqIO.read(fasta_io, "fasta") # iterate object and return SeqRecord Object

    return record

# EXAMPLE RUN
if __name__ == "__main__":
    print(fetch_uniprot_sequence('P22303'))

