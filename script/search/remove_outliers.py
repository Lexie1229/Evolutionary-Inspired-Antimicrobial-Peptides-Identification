# python -m miniprot -i input.tsv -o output.tsv

import argparse

# get options
parser = argparse.ArgumentParser(description="Remove outliers.")
parser.add_argument("-i", "--input", help="Path to the input file", required=True)
parser.add_argument("-o", "--output", help="Path to the output file", required=True)
parser.add_argument("-m", "--mode", help="Data formate, miniprot or blast", required=True)
args = parser.parse_args()

input_file = args.input
output_file = args.output
mode = args.mode

# miniprot
miniprot_thresholds = {
    'A': (0.56, 1.45),
    'C': (0.98, 1.07), 
    'E': (0.99, 1.01),
    'T': (1.00, 1.02),
    'TE': (1.00, 1.01)
}

# blast
blast_thresholds = {
    'A': (0.90, 1.08),
    'C': (0.98, 1.07),
    'E': (0.99, 1.02), 
    'T': (0.99, 1.02),
    'TE': (1.00, 1.01)
}

# get thresholds
def get_thresholds(mode):
    if mode == 'miniprot':
        thresholds = miniprot_thresholds
    elif mode == 'blast':
        thresholds = blast_thresholds
    else:
        raise ValueError("Invalid mode.")
    return thresholds

def main():
    thresholds = get_thresholds(mode)

    with open(input_file, 'r') as inputfile, open(output_file, 'w') as outputfile:
        for line in inputfile:
            # cols: domain_id|domain_len|strand|contig_id|contig_len|hit_start|hit_end
            cols = line.strip().split('\t')
            domain_id = cols[0]
            domain_len = int(cols[1])
            hit_start = int(cols[5])
            hit_end = int(cols[6])

            # domain：A/C/T/E/TE
            domain = domain_id.split('.')[0]

            # nucleotide length
            hit_len = hit_end - hit_start + 1
            domain_len = domain_len * 3

            # ratio of hit length to domain length
            ratio_hit_domain = round(hit_len / domain_len, 4)

            # remove outliers
            if thresholds.get(domain) and thresholds[domain][0] <=  ratio_hit_domain <= thresholds[domain][1]:
                outputfile.write(line)
            else:
                continue

if __name__ == "__main__":
    main()
