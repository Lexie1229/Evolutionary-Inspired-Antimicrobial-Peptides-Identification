# python input.tsv output.tsv

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

# blast
thresholds = {
    'A': (1.45, 2.23),
    'C': (1.54, 2.20),
    'E': (1.29, 2.37),
    'T': (1.76, 2.17),
    'TE': (1.51, 1.85)
}

def main():
    with open(input_file, 'r') as inputfile, open(output_file, 'w') as outputfile:
        for line in inputfile:

            # cols: domain_id|domain_len|strand|contig_id|contig_len|hit_start|hit_end|score
            cols = line.strip().split('\t')
            domain_id = cols[0]
            domain_len = int(cols[1])
            score = float(cols[7])

            # domain：A/C/T/E/TE
            domain = domain_id.split('.')[0]

            # ratio of score to domain length
            ratio_score_domain = round(score / domain_len, 4)

            # quality control
            if thresholds.get(domain) and thresholds[domain][0] <=  ratio_score_domain <= thresholds[domain][1]:
                outputfile.write(line)
            else:
                continue

if __name__ == "__main__":
    main()
