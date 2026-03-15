import requests
import xml.etree.ElementTree as ET
from pathlib import Path
import json
import re
from src.config import CLINVAR_DIR

def esearch_clinvar(gene, retmax=100):
    """Return ClinVar UID list for a gene."""
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "clinvar",
        "term": f"{gene}[gene]",
        "retmax": retmax,
        "retmode": "xml",
    }

    r = requests.get(url, params=params)
    r.raise_for_status()

    root = ET.fromstring(r.text)
    return [elem.text for elem in root.findall(".//Id")]


def uid_to_vcv(uid_list):
    """Convert ClinVar UIDs to VCV accessions."""

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    params = {
        "db": "clinvar",
        "id": ",".join(uid_list),
        "retmode": "xml",
    }

    r = requests.get(url, params=params)
    r.raise_for_status()

    root = ET.fromstring(r.text)

    vcv_list = []
    for elem in root.findall(".//DocumentSummary"):
        accession = elem.find(".//accession")
        if accession is not None:
            vcv_list.append(accession.text)

    return vcv_list


def fetch_vcv_xml(vcv_list):
    """Fetch XML records for VCV accessions."""

    xml_records = []

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    for vcv in vcv_list:

        params = {
            "db": "clinvar",
            "id": vcv,
            "rettype": "vcv",
            "retmode": "xml",
        }

        r = requests.get(url, params=params)

        if r.ok:
            xml_records.append(r.text)

    return xml_records


def parse_vcv(xml_records):
    """Parse ClinVar XML records into structured variant data."""

    variants = []

    for xml_str in xml_records:

        try:
            root = ET.fromstring(xml_str)
        except ET.ParseError:
            continue

        for record in root.findall(".//VariationArchive"):

            cs = record.find(".//Classifications/GermlineClassification/Description")
            clinical_significance = cs.text.strip() if cs is not None else None

            rs = record.find(".//Classifications/GermlineClassification/ReviewStatus")
            review_status = rs.text.strip() if rs is not None else None

            hgvs = {
                "protein_missense": [],
                "protein_frameshift": [],
                "protein_other": [],
                "coding": [],
                "genomic": [],
                "other": [],
            }

            for h in record.findall(".//HGVSlist/HGVS"):

                expr = (
                    h.find("ProteinExpression/Expression")
                    or h.find("NucleotideExpression/Expression")
                    or h.find("Expression")
                )

                if expr is None or not expr.text:
                    continue

                line = expr.text.split(":")[1].strip()

                if "p." in line:
                    if "fs" in line:
                        hgvs["protein_frameshift"].append(line)
                    elif re.match(r"p\.[A-Za-z]{3}\d+[A-Za-z]{3}", line):
                        hgvs["protein_missense"].append(line)
                    else:
                        hgvs["protein_other"].append(line)

                elif "c." in line or "n." in line:
                    hgvs["coding"].append(line)

                elif "g." in line:
                    hgvs["genomic"].append(line)

                else:
                    hgvs["other"].append(line)

            variants.append(
                {
                    "clinical_significance": clinical_significance,
                    "review_status": review_status,
                    "hgvs": hgvs,
                }
            )

    return variants


def get_clinvar_variants(gene, save_json=True):
    """Main pipeline for retrieving ClinVar variants for a gene."""

    uid_list = esearch_clinvar(gene)
    vcv_list = uid_to_vcv(uid_list)
    xml_records = fetch_vcv_xml(vcv_list)

    variants = parse_vcv(xml_records)

    if save_json and variants:

        out_file = CLINVAR_DIR / f"{gene}_clinvar.json"

        with open(out_file, "w") as f:
            json.dump(variants, f, indent=4)

    return variants


if __name__ == "__main__":
    get_clinvar_variants("TP53")