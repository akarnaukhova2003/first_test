import argparse
import os
import pandas as pd
from Bio import SeqIO
from parser_gb import read_csv

def load_extra_genes(extra_file):
    """
    Загружает словарь extra_genes из CSV.
    Формат CSV: id,gene
    """
    df = pd.read_csv(extra_file)
    extra_dict = {}
    for _, row in df.iterrows():
        rec_id = str(row["id"]).strip()
        gene = str(row["gene"]).strip()
        extra_dict[rec_id] = gene
    return extra_dict


def orf_coord(input_file, orf_map, remove_exceptions, extra_file=None):
    """
    Retrieves the coordinates of ORFs
    
    Input:
        input_file - file with nucleotide sequences in GenBank format
        orf_map - csv file annotation of orfs and their codes
        remove_exceptions - skip IDs from exceptions list
        extra_file - csv with additional mapping {id -> gene}
    Output:
        coord_file - file with coordinates
    """

    exceptions_file = "..\\sapovirus\\norovirus_exceptions.csv"
    orf_dict = read_csv(orf_map)

    orf_types_final = ["1A", "1B", "2"]

    out_file_name = os.path.splitext(input_file)[0] + "_orf.txt"
    print("Writing results to:", out_file_name)
    out_file = open(out_file_name, "w")

    out_file.write("id" + "," + ",".join(orf_types_final) + "\n")

    dict_coord = {}

    extra_genes = load_extra_genes(extra_file) if extra_file else {}

    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, "gb"))
        for rec in records:
            is_exception = False
            if remove_exceptions:
                with open(exceptions_file, "r") as exceptions_f:
                    for line in exceptions_f:
                        if line.strip() == rec.name:
                            print("Record", rec.name, "was skipped")
                            is_exception = True
            if is_exception:
                continue

            dict_coord[rec.name] = {}
            pol_count = 0

            for feature in rec.features:
                cod_start = int(feature.qualifiers.get("codon_start", [1])[0]) - 1

                if feature.type == "CDS":
                    product = None
                    if "product" in feature.qualifiers:
                        product = map_feature(feature.qualifiers["product"][0], orf_dict)
                    elif "gene" in feature.qualifiers:
                        product = map_feature(feature.qualifiers["gene"][0], orf_dict)

                    if not product:
                        continue

                    if product in orf_types_final and product not in dict_coord[rec.name]:
                        dict_coord[rec.name][product] = [
                            int(feature.location.start) + cod_start,
                            int(feature.location.end),
                        ]

            if rec.name in extra_genes:
                target_gene = extra_genes[rec.name].upper()
                if target_gene not in dict_coord[rec.name]:
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        found = False
                        if "product" in feature.qualifiers:
                            prod = feature.qualifiers["product"][0].upper()
                            if target_gene in prod:
                                found = True
                        elif "gene" in feature.qualifiers:
                            gene = feature.qualifiers["gene"][0].upper()
                            if target_gene == gene:
                                found = True
                        if found:
                            cod_start = int(feature.qualifiers.get("codon_start", [1])[0]) - 1
                            start = int(feature.location.start) + cod_start
                            end = int(feature.location.end)
                            dict_coord[rec.name][target_gene] = [start, end]
                            break

    # записываем файл
    for rec_id in dict_coord.keys():
        s = rec_id
        for orf in orf_types_final:
            if orf in dict_coord[rec_id]:
                st, e = dict_coord[rec_id][orf]
                s += "," + str(st) + "-" + str(e)
            else:
                s += ",NA-NA"
        s += "\n"
        out_file.write(s)

    out_file.close()


def map_feature(feature, feature_map):
    for k, v in feature_map.items():
        if feature.lower() == k.lower():
            return v
    return feature


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str, required=True,
                        help="Input GenBank file")
    parser.add_argument("-orf_map", "--orf_map_file", type=str, required=True,
                        help="CSV-file with short codes for ORFs")
    parser.add_argument("-extra", "--extra_file", type=str, required=False,
                        help="CSV with extra genes (id,gene)")
    parser.add_argument("-r", "--remove_exceptions", action="store_true",
                        help="Remove exceptions")

    args = parser.parse_args()
    orf_coord(args.input_file, args.orf_map_file, args.remove_exceptions, args.extra_file)
