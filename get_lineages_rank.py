from Bio import Entrez
import time
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", required=True, help="Входной CSV-файл")
parser.add_argument("--output", "-o", default="hosts_lineage_rank.csv", help="Выходной CSV-файл")
args = parser.parse_args()

Entrez.email = "adelina.kzn15@gmail.com"
df = pd.read_csv('sequences.csv')
hosts = df["Host"].dropna().unique().tolist()
def get_lineage_with_ranks(organism_name):
    stream = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
    record = Entrez.read(stream)
    stream.close()
    idlist = record.get("IdList") or []
    if not idlist:
        return None
    taxid = idlist[0]
    
    stream = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(stream)
    stream.close()
    
    lineage_list = records[0].get("LineageEx", [])
    lineage_dict = {item["Rank"]: item["ScientificName"] for item in lineage_list if item["Rank"] != "no rank"}
    lineage_dict[records[0]["Rank"]] = records[0]["ScientificName"]
    
    return lineage_dict
res = []
for i, h in enumerate(hosts):
    lineage_dict = get_lineage_with_ranks(h)
    res.append({"host": h, **(lineage_dict or {})})
    print(f"{i+1}/{len(hosts)}: {h} → {lineage_dict}")

df = pd.DataFrame(res)
df.to_csv("hosts_lineage_rank.csv", index=False, encoding='utf-8')