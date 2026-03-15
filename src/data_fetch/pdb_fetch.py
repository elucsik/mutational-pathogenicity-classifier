import requests
import os


def fetch_pdb_sequence(pdb_id, out_dir=".",fmt = "pdb"):
    """
    DOWNLOAD PDB/CIF STRUCTURE FILE FROM RSCB

    Args:
        pdb_id (str): PDB Identifier
        out_dir (str): Output directory where PDB/CIF will be saved
        fmt (str): file format ('pdb' or 'cif')

    Returns: 
        Download: pdb/cif file
        Path to downloaded file (str)
    """    
    #API CALL TO RSCB
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.{fmt}"
    r = requests.get(url)
    r.raise_for_status() #flag unsuccessful responses as errors

    #MAKE DIRECTORY FOR PDB/CIF FILE
    os.makedirs(out_dir, exist_ok=True)
    file_path = os.path.join(out_dir, f"{pdb_id}.{fmt}") #build file path

    #WRITE FILE
    with open(file_path, 'w') as f:
        f.write(r.text)
    
    return file_path

if __name__ == "__main__":
    print(fetch_pdb_sequence('4PQE', 'data_fetch', 'cif'))