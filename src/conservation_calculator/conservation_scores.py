import os
import subprocess
from Bio import SeqIO, Entrez, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from data_fetch import uniprot_fetch
import math
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
import json
project_root = '/Users/ethanlucsik/Desktop/python_work/cancer_mutation_project'
output_dir = os.path.join(project_root, "results")
os.makedirs(output_dir, exist_ok=True)


def blast_uniprot_remote(uniprot_id, out_xmlfile="blast_results.xml", out_fastafile="msa_input.fasta"):
    """
    Run BLASTp remotely using NCBI servers and save XML results.
    
    uniprot_id
        Biopython SeqRecord object (protein sequence)
    output_file : str
        Path to save the BLAST XML results
    """
    # 
    seq_record = uniprot_fetch.fetch_uniprot_sequence(uniprot_id)

    # Run remote BLAST
    result_handle = NCBIWWW.qblast(
        program="blastp",       # Protein BLAST
        database="swissprot",   # SwissProt database
        sequence=seq_record.seq,
        hitlist_size=500        # Max number of hits
    )
    
    # Save results to XML
    with open(f"{output_dir}/{out_xmlfile}", "w") as f:
        f.write(result_handle.read())
    result_handle.close()
    
    # Parse BLAST XML
    with open(f"{output_dir}/{out_xmlfile}") as f:
        blast_records = NCBIXML.parse(f)
        blast_record = next(blast_records)

        seq_records = []
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                seq_records.append(SeqRecord(Seq(hsp.sbjct), id=alignment.hit_def, description=""))
                break
            
        seq_records.append(SeqRecord(seq_record.seq, id=seq_record.id, description=""))
        SeqIO.write(seq_records, f"{output_dir}/{out_fastafile}", "fasta")

    return out_fastafile


def run_msa(out_fastafile, out_msafile="msa_output.fasta", verbose=False, force=False):
    out_path = f"{output_dir}/{out_msafile}"
    if os.path.exists(out_path) and not force:
        print(f"File {out_msafile} already exists. Skipping MSA.")
        alignment = AlignIO.read(out_path, "fasta")
        return alignment
    else:
        mafft_cline = MafftCommandline(input=f"{output_dir}/{out_fastafile}")
        if verbose:
            mafft_cline.set_parameter("quiet", False)
        
        stdout, stderr = mafft_cline()
        if stderr:
            print(f"Std Err:\n{stderr}")
        
        # Write alignment to file
        with open(f"{output_dir}/{out_msafile}", "w") as f:
            f.write(stdout)

            alignment = AlignIO.read(f"{output_dir}/{out_msafile}", "fasta")
        
        return alignment


def compute_per_residue_entropy(alignment, out_entropyfile=True):
    scores = []
    for i in range(alignment.get_alignment_length()):
        sequences = [sequ.seq[i] for sequ in alignment]
        freq = Counter(sequences)
        if '-' in freq: del freq['-']
        total = sum(freq.values())
        if total == 0:
            scores.append(0.0)
        else:
            entropy = -sum((f/total) * math.log2(f/total) for f in freq.values())
            scores.append(entropy)
    #normalize
    conservation = [1 - (s / math.log2(20)) for s in scores]
    
    if out_entropyfile:
        with open(f"{output_dir}/res_entropy_results.txt", "w") as f:
            json.dump(conservation, f, indent=4)
    
    return conservation


if __name__ == "__main__":
    #out_fastafile = blast_uniprot_remote('P22303')
    alignment = run_msa("msa_input.fasta")
    compute_per_residue_entropy(alignment, True)

