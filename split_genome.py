import argparse
import copy
import os
import pandas as pd
import re
from Bio import SeqIO


def split_genome_corona(input_file, coord_file):
    '''
    Splits nucleotide sequences from input_file into ORFs using their coordinates from 
    coord_file

    Input:
        input_file - file with nucleotide sequences in fasta-format
        coord_file - text file with coordinates of NTRs
        


    '''
    

    coord_df = pd.read_csv(coord_file, sep=",", index_col=0)
    orfs = list(coord_df.columns)


    print(coord_df)
    dict_orfs = {}
    for orf in orfs:
        dict_orfs[orf] = list()

    seqs = SeqIO.parse(open(input_file), format='fasta')

    for seq in seqs:
        # accession number of sequence
        acc = seq.id.split("_")[0]
        
        if acc == 'NC' or acc == 'AC':
            acc = '_'.join([acc,seq.id.split("_")[1]])
        #print(acc)
        for orf in orfs:
            # checks whether all orfs have known coordinates
            
            if ((acc in coord_df.index) and (coord_df.loc[acc] != 'NA-NA').all()):

                # coordinates
                cd = coord_df[orf][acc]
                s, e = [int(x) for x in cd.split('-')]
                seq_orf = copy.deepcopy(seq)
                seq_orf.seq = seq.seq[s:e]
                dict_orfs[orf].append(seq_orf)
            else:
                print('{} not in coord file'.format(acc))
    out_temp = os.path.splitext(input_file)[0]
    for orf in dict_orfs.keys():
        outfile = out_temp + '_' + orf + '.fasta'
        SeqIO.write(dict_orfs[orf], outfile, "fasta")



def split_genome(input_file, coord_file):
    '''
    Splits nucleotide sequences from input_file into 5'NTR,
    one coding region (everything between NTRs), 3'NTR using coordinates of these region from 
    coord_file

    Input:
        input_file - file with nucleotide sequences in fasta-format
        coord_file - text file with coordinates of NTRs
        


    '''
    out_temp = os.path.splitext(input_file)[0]
    
    output_5utr_n = out_temp + '_5utr.fasta'
    output_cds_n = out_temp + '_coding.fasta'
    output_3utr_n = out_temp + '3utr.fasta'
    
    print(output_5utr_n)
    
    list_5utr = list()
    list_cds = list()
    list_3utr = list()

    coord_df = pd.read_csv(coord_file, sep=",", index_col=0)

    seqs = SeqIO.parse(open(input_file), format='fasta')

    for seq in seqs:
        # accession number of sequence
        acc = seq.id.split("_")[0]
        if acc == 'NC' or acc == 'AC':
            acc = '_'.join([acc,seq.id.split("_")[1]])
        if acc in coord_df.index:
            if coord_df['st_5utr'][acc] !=0 and coord_df['e_5utr'][acc]!= 0:
                seq_5utr = copy.deepcopy(seq)
                seq_5utr.seq = seq.seq[:coord_df['e_5utr'][acc]]
                list_5utr.append(seq_5utr)
            seq_cds = copy.deepcopy(seq)
            seq_cds.seq = seq.seq[coord_df['st_cds'][acc]-1:coord_df['e_cds'][acc]]
            list_cds.append(seq_cds)
            
            if coord_df['st_3utr'][acc] !=0 and coord_df['e_3utr'][acc]!= 0:
                seq_3utr = copy.deepcopy(seq)
                seq_3utr.seq = seq.seq[coord_df['st_3utr'][acc]:coord_df['e_3utr'][acc]]
                list_3utr.append(seq_3utr)
        else:
            print('{} not in coord file'.format(acc))
    SeqIO.write(list_5utr, output_5utr_n, "fasta")
    SeqIO.write(list_cds, output_cds_n, "fasta")
    SeqIO.write(list_3utr, output_3utr_n, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file in fasta-format", required=True)
    parser.add_argument("-coord", "--coord_file", type=str,
                        help="Csv-file with coordinates of 5\'UTR, coding region and 3\'UTR", required=True)
    args = parser.parse_args()

    split_genome_corona(args.input_file, args.coord_file)