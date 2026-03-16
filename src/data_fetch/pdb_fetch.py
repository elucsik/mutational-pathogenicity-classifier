import requests
from src.config import PDB_DIR


def fetch_pdb_structure(pdb_id, fmt="cif", save_file=True):
    """
    Download a structure file from the RCSB PDB.

    Args:
        pdb_id (str): PDB identifier
        fmt (str): file format ("pdb" or "cif")
        save_file (bool): optionally save the structure file

    Returns:
        str: path to downloaded file if saved
    """

    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.{fmt}"

    r = requests.get(url)
    r.raise_for_status()

    if save_file:
        file_path = PDB_DIR / f"{pdb_id}.{fmt}"

        with open(file_path, "w") as f:
            f.write(r.text)

        return file_path

    return r.text


if __name__ == "__main__":
    path = fetch_pdb_structure("1TUP")
    print(path)