import argparse
import copy
import os
import pandas as pd
import re
from Bio import SeqIO
from parser_gb import read_csv

def orf_coord(input_file, orf_map, remove_exceptions):
    '''
    Retrieves the coordinates of ORFs

    Input:
        input_file - file with nucleotide sequences in genbank-format
        orf_map - csv file annotation of orfs and their codes
    Output:
        coord_file - file with coordinates
    '''

    exceptions_file = '..\\sapovirus\\norovirus_exceptions.csv'
    orf_dict = read_csv(orf_map)
    orf_types_final = ['1A', '1B', '2']

    out_file_name = os.path.splitext(input_file)[0] + '_orf.txt'
    print(out_file_name)
    out_file = open(out_file_name, 'w')
    out_file.write('id' + ',' + ','.join(orf_types_final) + '\n')
    dict_coord = {}
    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, 'gb'))
        for rec in records:
            is_exception = False
            if remove_exceptions:
                with open(exceptions_file, 'r') as exceptions_f:
                     for line in exceptions_f:
                        if line.strip() == rec.name:
                            print('Record', rec.name, 'was skipped')
                            is_exception = True
            if is_exception == True:
                continue
            dict_coord[rec.name] = {}
            pol_count = 0
            for feature in rec.features:
                print(feature)
                if 'codon_start' in feature.qualifiers.keys():
                    cod_start = int(feature.qualifiers['codon_start'][0]) - 1
                    if cod_start<0:
                        print('here')
                else:
                    cod_start = 0

                if feature.type != 'CDS':
                    continue
                    

                raw_annotations = []
                for el in ('product', 'gene', 'note'):
                    raw_annotations.extend(feature.qualifiers.get(el, []))

                mapped = []
                seen = set()
                for a in raw_annotations:
                    if a is None:
                        continue
                    ma = map_feature(a, orf_dict)
                    if ma not in seen:
                        mapped.append(ma)
                        seen.add(ma)
                if not mapped:
                    continue
                
                for ann in mapped:  
                    if ann in dict_coord[rec.name]:
                        continue
                    def simple_coords(loc):
                        return [int(loc.start) + cod_start, int(loc.end)]
                    success = False
                    if ann in orf_types_final:
                        dict_coord[rec.name][ann] = simple_coords(feature.location)
                        success = True
                    elif ann == '1AB':
                        if pol_count == 1:
                            continue
                        if len(feature.location.parts) > 1:
                            dict_coord[rec.name]['1A'] = simple_coords(feature.location.parts[0])
                            dict_coord[rec.name]['1B'] = simple_coords(feature.location.parts[1])
                        elif len(feature.location.parts) == 1:
                            coords = simple_coords(feature.location.parts[0])
                            dict_coord[rec.name]['1B'] = coords
                        pol_count += 1
                        success = True

                    elif ann == '1AB_ORF':
                        coords = simple_coords(feature.location)
                        if pol_count == 0:
                            dict_coord[rec.name]['1A'] = coords
                            pol_count += 1
                        elif pol_count == 1:
                            dict_coord[rec.name]['1B'] = coords
                            pol_count += 1
                        else:
                            print(f"Too many polymerase ORFs for {rec.name}")
                        success = True

                    else:
                        pass
                    if success:
                        break  

    # Запись файла
    for rec_id, orfs in dict_coord.items():
        s = rec_id
        for orf in orf_types_final:
            if orf in orfs:
                s += f",{orfs[orf][0]}-{orfs[orf][1]}"
            else:
                s += ",NA-NA"
        s += "\n"
        out_file.write(s)
    out_file.close()


def map_feature(feature, feature_map):
    '''
    feature_map - dictionary, e.g. feature_map['Italia']='ITA'
    feature - возможный ключ для feature_map
    '''
    for k, v in feature_map.items():
        try:
            if feature.lower() == k.lower():
                return v
        except Exception:
            continue
    return feature


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str, required=True, help="Input file")
    parser.add_argument("-orf_map", "--orf_map_file", type=str, required=True, help="CSV with ORF codes")
    parser.add_argument("-r", "--remove_exceptions", action="store_true", help="Remove exceptions")
    args = parser.parse_args()

    orf_coord(args.input_file, args.orf_map_file, args.remove_exceptions)
