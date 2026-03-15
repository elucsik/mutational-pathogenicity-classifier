import requests
import xml.etree.ElementTree as ET #for parsing the xml format
#NOTE: ET always evaluates to True for OR statements --> Use mutliple if statements when .find()
import os
import json
import re
project_root = os.path.abspath('/Users/ethanlucsik/Desktop/python_work/cancer_mutation_project')
output_dir = os.path.join(project_root, "data_fetch")
os.makedirs(output_dir, exist_ok=True)

def return_clinvar_esearch_uid(gene_symbol, retmax=100):
    """
    USE ESEARCH TO FIND VARIATION ID's FOR CLINVAR GENE

    Args:
        gene_symbol (str): Gene symbol and c. or p. (coding DNA or protein), e.g. PTEN c.34A>C
        retmax (str): max returned results

    Returns: 
        Variation ID's for provided gene (str)
    """    
    #API CALL TO NCBI FOR ID's
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_params = {
        "db": "clinvar",
        "term": f"{gene_symbol}[gene]",
        "retmax": retmax,
        "retmode": "xml", #xml file for xml elementtree
    }

    #CONVERT & RETURN ID RESPONSE
    r = requests.get(url, params=search_params)
    r.raise_for_status()

    record = ET.fromstring(r.text)
    id_list = [elem.text for elem in record.findall(".//Id")]
    
    return id_list, gene_symbol
    

def return_clinvar_uid_to_vcv(uid_list, gene_symbol):
    """
    Convert a list of UIDs from ClinVar esearch into VCV accessions

    Args:
        uid_list (list[str]): List of ClinVar UIDs

    Returns: 
        Returns:
        list[str]: List of VCV accessions
    """    
    #API CALL TO NCBI
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    search_params = {
        "db": "clinvar",
        "id": ",".join(uid_list),
        "retmode": "xml",
    }

    #CONVERT & RETURN XML RECORD RESPONSE
    r = requests.get(url, params=search_params)
    r.raise_for_status()
    record = ET.fromstring(r.text)

    vcv_list = []
    for elem in record.findall(".//DocumentSummary"):
            accession = elem.find(".//accession")
            if accession is not None:
                vcv_list.append(accession.text)

    return vcv_list, gene_symbol


def return_clinvar_efetch_vcv(vcv_list, gene_symbol):
    """
    Fetch ClinVar VariationArchive XML records for a list of VCV accessions
    """
    xml_list = []
    for vcv in vcv_list:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "clinvar",
            "id": vcv,            # VCV accession
            "rettype": "vcv",     # <-- must be 'vcv'
            "retmode": "xml"
        }
        response = requests.get(url, params=params)
        if response.ok:
            xml_list.append(response.text)
            print(f"{vcv} → {response.text[:200]}")
        else:
            print(f"Error fetching {vcv}: {response.status_code}")
    
    return xml_list, gene_symbol


def parse_clinvar_vcv(vcv_xml, gene_symbol="GENE"):
    all_variants = []

    for xml_str in vcv_xml:
        try:
            parsed_xml = ET.fromstring(xml_str)
        except ET.ParseError:
            print(f"Could not parse XML: {xml_str[:100]}")
            continue

        # loop through all VariationArchive entries in this XML
        for record in parsed_xml.findall(".//VariationArchive"):
            cs = record.find(".//Classifications/GermlineClassification/Description")
            clinical_significance = cs.text.strip() if cs is not None else None

            rs = record.find(".//Classifications/GermlineClassification/ReviewStatus")
            review_status = rs.text.strip() if rs is not None else None

            hgvs_categorized = {
                "protein_missense": [],
                "protein_frameshift": [],
                "protein_other": [],
                "coding": [],
                "genomic": [],
                "other": []
            }
            for h in record.findall(".//HGVSlist/HGVS"):
                expr = h.find("ProteinExpression/Expression")
                if expr is None:
                    expr = h.find("NucleotideExpression/Expression")
                if expr is None:
                    expr = h.find("Expression")
                if expr is not None and expr.text:
                    line = expr.text.strip()
                    line = line.split(":")[1]

                    if "p." in line:
                        if "fs" in line:
                            hgvs_categorized["protein_frameshift"].append(line)
                        elif re.match(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", line):
                            hgvs_categorized["protein_missense"].append(line)
                        else:
                            hgvs_categorized["protein_other"].append(line)
                    elif "c." in line:
                        hgvs_categorized["coding"].append(line)
                    elif "n." in line:
                        hgvs_categorized["coding"].append(line)
                    elif "g." in line:
                        hgvs_categorized["genomic"].append(line)
                    else:
                        hgvs_categorized["other"].append(line)

            all_variants.append({
                "clinical_significance": clinical_significance,
                "review_status": review_status,
                "hgvs": hgvs_categorized
            })
    
    print(f"✅ Retrieved {len(all_variants)} variants.")
    return all_variants


def clinvar_get_protein_variants(gene_symbol, output_json=True):
    uid_list, gene_symbol = return_clinvar_esearch_uid(gene_symbol)
    vcv_list, gene_symbol = return_clinvar_uid_to_vcv(uid_list, gene_symbol)
    xml_response_list, _ = return_clinvar_efetch_vcv(vcv_list, gene_symbol)
    all_variants = parse_clinvar_vcv(xml_response_list, gene_symbol)

    if output_json:
        if all_variants:
            out_file = os.path.join(output_dir, f"{gene_symbol}_clinvar.json")#creating a new file
            with open(out_file, 'w') as f:
                json.dump(all_variants, f, indent=4)
        else:
            print(f"No protien variants found for {gene_symbol}")
    
    return all_variants


if __name__ == "__main__":
    clinvar_get_protein_variants("TP53")

    

    