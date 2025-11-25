from Bio import Entrez
import time
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", required=True, help="Входной CSV-файл")
parser.add_argument("--output", "-o", default="hosts_lineage.csv", help="Выходной CSV-файл")
args = parser.parse_args()

Entrez.email = "adelina.kzn15@gmail.com"
df = pd.read_csv('sequences.csv')
hosts = df["Host"].dropna().unique().tolist()
def get_lineage(organism_name):
            stream = Entrez.esearch(db="taxonomy", term=organism_name, retmode="xml")
            record = Entrez.read(stream)
            stream.close()
            stream = None
            idlist = record.get("IdList") or []
            if not idlist:
                return None
            taxid = idlist[0]
            stream = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(stream)
            stream.close()
            lineage = records[0].get("Lineage", "")
            return lineage.split("; ") if lineage else None
res = []
for i, h in enumerate(hosts):
    lineage = get_lineage(h)
    res.append({"host": h, "lineage": lineage})
    print(f"{i+1}/{len(hosts)}: {h} → {lineage}")
max_len = max(len(d.get('lineage') or []) for d in res)
df = pd.DataFrame([
    {'host': d.get('host'),**{f'level_{i+1}': (d.get('lineage') or [None] * max_len)[i] if i < len(d.get('lineage') or []) else None for i in range(max_len)}
    } for d in res])
df.to_csv("hosts_lineage.csv", index=False, encoding='utf-8')