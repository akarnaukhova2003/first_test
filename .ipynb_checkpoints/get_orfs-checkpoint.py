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
    # dictionary with annotations of ORFs
    orf_dict = read_csv(orf_map)
    # possible ORFs
    #orf_types = ['1A', '1B', '1AB', '1AB_ORF', 'S', 'E', 'M', 'N']

    #orf_types = ['ORF1', 'ORF2', 'ORF3']
    #orf_types = ['1A', '1AB', '1B', '2']
    #orf_types = ['1A', '1B', '1AB', '2', '?', 'X']
    # list of ORFs in final table
    #orf_types_final = ['1A', '1B', 'S', 'E', 'M', 'N']
    #orf_types_final = ['ORF1', 'ORF2', 'ORF3']
    #orf_types_final = ['1A', '1B', '2',]
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
            if remove_exceptions == True:
                with open(exceptions_file, 'r') as exceptions_f:
                    for line in exceptions_f:
                        if line.strip() == rec.name:
                            print('Record', rec.name, 'was skipped')
                            is_exception = True
            if is_exception == True:
                continue
            dict_coord[rec.name] = {}
            # counter for polymerase A and B genes
            # if counter==0, haven't met A gene
            pol_count = 0
            print(rec.name)
            for feature in rec.features:
                if 'codon_start' in feature.qualifiers.keys():
                    cod_start = int(feature.qualifiers['codon_start'][0]) - 1
                    if cod_start<0:
                        print('here')
                else:
                    cod_start = 0
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers.keys():
                        product = map_feature(feature.qualifiers['product'][0], orf_dict)
                        #   print(feature.qualifiers['product'][0], product)
                        
                        if product not in dict_coord[rec.name].keys():
                            print(product)
                            if product in orf_types_final:
                                dict_coord[rec.name][product] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                            elif product == '1AB':
                                if pol_count == 1:
                                    continue
                                else:
                                    print(feature.location.parts)
                                    print([int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)])
                                    dict_coord[rec.name]['1A'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                    if len(feature.location.parts) > 1:
                                        dict_coord[rec.name]['1B'] = [int(feature.location.parts[1].start) + cod_start, int(feature.location.parts[1].end)]
                                    pol_count += 1
                            elif product == '1AB_ORF':
                                if pol_count == 0:
                                    dict_coord[rec.name]['1A'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.name]['1B'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                else:
                                    print('Couldn\'t find annotation \'product\' qualifier for {}'.format(rec.name))
#                            if product in orf_types_final:
#                                dict_coord[rec.name][product] = [int(feature.location.start) + cod_start, int(feature.location.end)]
#                            else:
#                                print('Couldn\'t find annotation \'product\' qualifier for {}'.format(rec.id))

                        
                    elif 'gene' in feature.qualifiers.keys():
                        gene = map_feature(feature.qualifiers['gene'][0], orf_dict)
                        if gene not in dict_coord[rec.name].keys():
                            if gene in orf_types_final:
                                dict_coord[rec.name][gene] = [int(feature.location.start) + cod_start, int(feature.location.end)]

                            elif gene == '1AB':
                                if pol_count == 1:
                                    continue
                                else:
                                    dict_coord[rec.name]['1A'] = [int(feature.location.parts[0].start) + cod_start, int(feature.location.parts[0].end)]
                                    dict_coord[rec.name]['1B'] = [int(feature.location.parts[1].start) + cod_start, int(feature.location.parts[1].end)]
                                    pol_count += 1
                            elif gene == '1AB_ORF':
                                if pol_count == 0:
                                    dict_coord[rec.name]['1A'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.name]['1B'] = [int(feature.location.start) + cod_start, int(feature.location.end)]
                                else:
                                    print('Couldn\'t find annotation in \'gene\' qualifier for {}'.format(rec.name))

#                            else:
#                                print('Couldn\'t find annotation in \'gene\' qualifier for {}'.format(rec.id))

    for id in dict_coord.keys():
        # string to write to out_file
        s = id 
        for orf in orf_types_final:
            if orf in dict_coord[id].keys():
                st = dict_coord[id][orf][0]
                e = dict_coord[id][orf][1]
                s = s + ',' + str(st) + '-' + str(e)
            else:
                s = s + ',NA-NA'
        s += '\n'
        out_file.write(s)
    #print(dict_coord)
    out_file.close()


def map_feature(feature, feature_map):
    '''
    feature_map - dictionary, e.g. feature_map['Italia']='ITA'
    feature - possible key from feature_map
    return value if key is in feature_map
    '''
    for k, v in feature_map.items():
        if feature.lower() == k.lower():
            return v
    return feature

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-orf_map", "--orf_map_file", type=str,
                        help="Csv-file with short codes for ORFs", required=True)
    parser.add_argument("-r", "--remove_exceptions",
                        help="Remove exceptions", action="store_true")
    args = parser.parse_args()

    orf_coord(args.input_file, args.orf_map_file, args.remove_exceptions)