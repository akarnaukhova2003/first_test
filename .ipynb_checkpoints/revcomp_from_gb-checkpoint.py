from Bio import SeqIO
import argparse

def main(gb_file, fasta_file, output_file):
    compl_ids = []
    with open(gb_file) as handle:
        records = SeqIO.parse(handle, "genbank")
        for rec in records:
            accessions = rec.annotations.get("accessions", [])
            if not accessions:
                continue
            accession = accessions[0]
            for feature in rec.features:
                if feature.type == "CDS" and '(-)' in str(feature.location):
                    compl_ids.append(accession)

    print(f"Найдено {len(compl_ids)} комплементарных ID")

    new_records = []
    for seq in SeqIO.parse(fasta_file, "fasta"):
        acc = seq.id.split("_")[0]
        if acc in compl_ids:
            seq.seq = seq.seq.reverse_complement()
            seq.id = seq.id + "_revcomp"
            seq.description = "reverse complement"
        new_records.append(seq)

    SeqIO.write(new_records, output_file, "fasta")
    print(f"Результат сохранён в {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reverse complement sequences from FASTA based on GenBank CDS strand")
    parser.add_argument("-gb", "--genbank", required=True, help="Input GenBank file")
    parser.add_argument("-fa", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    main(args.genbank, args.fasta, args.output)
