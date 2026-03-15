from Bio.PDB import PDBParser, PPBuilder, MMCIFParser, DSSP, PDBIO
import py3Dmol
import sys
import os
import webbrowser
import json
project_root = os.path.abspath('/Users/ethanlucsik/Desktop/python_work/cancer_mutation_project/cancer_mutation_project')
sys.path.append(project_root)
from data_fetch import pdb_fetch

output_dir = os.path.join(project_root, "structural_features")
os.makedirs(output_dir, exist_ok=True)


def parse_pdb_file(pdb_id, file_type, data_dir): #data_dir must point to PDB file in data_fetch
    """
    Fetch/parse PDB/CIF File
    Returns Structure, pdb_id, and file_path
    """
    #SET ID & FILETYPE
    file_path = pdb_fetch.fetch_pdb_sequence(pdb_id, data_dir, file_type)

    #PARSE BASED ON FILETYPE
    if '.cif' in file_path:
        parser = MMCIFParser(QUIET=True)
    elif '.pdb' in file_path:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")
    
    structure = parser.get_structure(pdb_id, file_path)
    return structure, pdb_id, file_path


def return_meta_data(structure, pdb_id, file_path, output_json = True):
    """
    Return dictionary of meta_data (models, chains, residues, atoms) for given pdb_id
    """
    #CREATE META_DATA DIC
    meta_data = {
        "pdb_id": pdb_id,
        "num_models": len(structure),
        "models": {}
    }
    # Models
    for model in structure:
        model_dict = {"chains":{}}
        # Chains
        meta_data["models"][model.id] = model_dict
        for chain in model:
            chain_dict = {"residues": []}
            model_dict["chains"][chain.id] = chain_dict
            # Residues
            for residue in chain:
                residue_dict = {
                    "resname": residue.resname,
                    "id":residue.id,
                    # Atoms
                    "atoms": [(atom.name, atom.coord.tolist()) for atom in residue]

                }
                chain_dict["residues"].append(residue_dict)
                # Atoms

    #META_DATA JSON IN ./structural_features DIRECTORY
    if output_json:
        out_file = os.path.join(output_dir, f"{pdb_id}_meta_data.json") #creating a new file
        with open(out_file, 'w') as f:
            json.dump(meta_data, f, indent=4)

    return meta_data


def return_secondary_structure_data(structure, pdb_id, file_path, output_json = True):
    """
    Return dictionary of secondary structure data (
    chain_id, res_id, residue_name, ss_type, asa)
    """
    #SECONDARY STRUCTURE DETAILS
    #NOTE: DSSP Tuple Indece
    # - 0 --> Residue aa (one letter code)
    # - 1 --> Sec struc data (table below)
    # - 2 --> Relative Accessible SA
    # - 3 --> Phi Torsion Ang
    # - 4 --> Psi Torsion Ang
    # - 5 --> H Bond info
    model = structure[0]
    d = DSSP(model, file_path)
    ss_dict = {
        "pdb_id": pdb_id,
        "secondary_structure": {}
    }
    for key, val in d.property_dict.items():
        chain_id = key[0]
        res_id = key[1][1]
        residue_name = val[0]
        ss_type = val[1]
        asa = val[2]
            
        ss_entry = {
            "residue_id": res_id,
            "residue_name": residue_name,
            "secondary_structure": ss_type,
            "accessible_sa": asa,
        }
        if chain_id not in ss_dict["secondary_structure"]:
            ss_dict["secondary_structure"][chain_id] = []
        ss_dict["secondary_structure"][chain_id].append(ss_entry)

    """
    DSSP Codes Assignment
    _____________________________________
    Code        Sec Struc Type
    _____________________________________
    H           α-helix (4-turn helix)
    B           Isolated β-bridge residue
    E           Extended strand, participates in β-ladder
    G           3-helix (310 helix)
    I           5-helix (π-helix)
    P           poly-proline helix (κ-helix)
    T           Turn
    S           Bend
    -           Loop/irregular/coil
    More at https://pdb-redo.eu/dssp/about
    """

    if output_json:
        out_file = os.path.join(output_dir, f"{pdb_id}_secondary_structure.json")#creating a new file
        with open(out_file, 'w') as f:
            json.dump(ss_dict, f, indent=4)

    return ss_dict


def return_aa_sequence(structure, pdb_id, file_path, output_json=True):
    """
    Return single letter aa sequences per chain for structure
    """
    ppb = PPBuilder()
    chain_sequences = {}
    for pp in ppb.build_peptides(structure):
        seq = str(pp.get_sequence())
        chain_id = pp[0].get_parent().id # pp[0] = first residue in peptide, .getparent() chain object that contains residue pp[0], .id is parent chain id
        chain_sequences[chain_id] = seq
    
    if output_json:
        out_file = os.path.join(output_dir, f"{pdb_id}_chain_sequences.json")#creating a new file
        with open(out_file, 'w') as f:
            json.dump(chain_sequences, f, indent=4)

    return chain_sequences


def return_pdb_with_asa(structure, pdb_id, file_path, output_asa_pdb=True):
    model = structure[0]
    dssp = DSSP(model, file_path)
    for (chain_id, res_id), dssp_vals in dssp.property_dict.items():
        asa = dssp_vals[3]*100 # ASA value
        if asa is None or isinstance(asa, str):
            try:
                asa = float(asa)
            except Exception:
                asa = 0.0
        
        if chain_id in structure[0]:
            chain = structure[0][chain_id]
            if res_id in chain:
                for atom in chain[res_id]: #replace bfactors
                    atom.bfactor = float(asa) # asa value = to atom.bfactor for pymol mapping

    if output_asa_pdb:
        io = PDBIO()
        io.set_structure(structure)
        output_path = f"{output_dir}/{pdb_id}_pdb_asa.pdb"
        io.save(output_path)
        """
        load 1HZH_with_ASA.pdb
        spectrum b, red_white_blue
        """
        return output_path


def extract_features(pdb_id, file_type, data_dir, output_json=True, output_asa_pdb=True):
    structure, pdb_id, file_path = parse_pdb_file(pdb_id, file_type, data_dir)
    meta = return_meta_data(structure, pdb_id, file_path, output_json)
    ss = return_secondary_structure_data(structure, pdb_id, file_path, output_json)
    aa = return_aa_sequence(structure, pdb_id, file_path, output_json)
    return_pdb_with_asa(structure, pdb_id, file_path, output_asa_pdb)
    return {"meta_data": meta, "secondary_structure": ss, "aa_sequence": aa}


"""
NOTE: ASA File for PyMol can be manually exported if Pymol is installed or can be ran through script
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-qc'])

cmd.load(f"{pdb_id}_pdb_asa.pdb")
cmd.spectrum("b", "blue_white_red")
cmd.png("asa_colored.png", 1200, 900, dpi=300, ray=1)
cmd.save(f"{output_dir}/1HZH_colored.pdb")
"""


def return_pymol_html(output_path, html_out="pdb_view.html", auto_open = True):
    with open(output_path, "r") as f:
        pdb_data = f.read()
    
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_data, "pdb")
    view.setStyle(
        {"cartoon": {
            "color": "spectrum",
            "colorscheme": {"prop": "bfactor", "gradient": "roygb"}
        }}
    )
    # bfactors = []
    # with open(output_path) as f:
    #     for line in f:
    #         if line.startswith("ATOM") or line.startswith("HETATM"):
    #             bfactors.append(float(line[60:66]))
    # min_b, max_b = min(bfactors), max(bfactors)
    # view.addLegend({
    #     "gradient": "roygb",
    #     "min": min_b, "max": max_b,
    #     "title": "ASA (Å²)"
    # })
  
    # color by ASA stored in B-factors
    view.zoomTo()
    
    with open(html_out, "w") as out:
        out.write(view._make_html())
    
    print(f"3D viewer saved → {html_out}")

    if auto_open:
        webbrowser.open(f"file://{os.path.abspath(html_out)}")


def view_pdb_asa_html(pdb_id, file_type="pdb", data_dir="structural_features"):
    file_path = pdb_fetch.fetch_pdb_sequence(pdb_id, data_dir, file_type)

    if file_type == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file_path)
    pdb_asa_file = return_pdb_with_asa(structure, pdb_id, file_path, output_dir)
    return_pymol_html(pdb_asa_file, html_out=f"{output_dir}/{pdb_id}_asa_view.html")


if __name__ == "__main__":
    #view_pdb_asa_html("1HZH")
    extract_features('1HZH', 'pdb','structural_features',False, False)