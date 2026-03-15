from pathlib import Path

# project root
ROOT = Path(__file__).resolve().parents[1]

# main directories
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
NOTEBOOKS_DIR = ROOT / "notebooks"

# specific data sources
CLINVAR_DIR = DATA_DIR / "clinvar"
UNIPROT_DIR = DATA_DIR / "uniprot"
PDB_DIR = DATA_DIR / "pdb"

# output directories
STRUCT_FEATURE_DIR = RESULTS_DIR / "structural_features"
CONSERVATION_DIR = RESULTS_DIR / "conservation"

# ensure directories exist
for d in [
    DATA_DIR,
    RESULTS_DIR,
    CLINVAR_DIR,
    UNIPROT_DIR,
    PDB_DIR,
    STRUCT_FEATURE_DIR,
    CONSERVATION_DIR,
]:
    d.mkdir(parents=True, exist_ok=True)
