# python alphafold_fasta2json.py -i input.fa -o output.json

import argparse
import json
from Bio import SeqIO

def read_fasta(fasta_file):
    records = []
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), start=1):
        name_id = record.id
        fasta_data = {
            "name": name_id,
            "modelSeeds": [],
            "sequences": [
                {
                    "proteinChain": {
                        "sequence": str(record.seq),
                        "count": 1
                    }
                }
            ]
        }
        records.append(fasta_data)
    return records

def json_output(records, json_file):
    with open(json_file, 'w') as output_file:
        json.dump(records, output_file, indent=4)


def main():
    parser = argparse.ArgumentParser(description='Convert FASTA to JSON suitable for Alphafold.')
    parser.add_argument('-i', '--input', type=str, required=True, help='fasta file')
    parser.add_argument('-o', '--output', type=str, required=True, help='json file')

    args = parser.parse_args()
    records = read_fasta(args.input)
    json_output(records, args.output)

if __name__ == '__main__':
    main()
