# python -t 0.95 -i input.tsv -o output.tsv

import argparse
from collections import defaultdict

# get options
parser = argparse.ArgumentParser(description="Combine matches at the same location within a specific overlap threshold.")
parser.add_argument("-i", "--input", help="Path to the input file", required=True)
parser.add_argument("-o", "--output", help="Path to the output file", required=True)
parser.add_argument("-t", "--threshold", type=float, default=0.95, help="Overlap threshold (default: 0.95)")
args = parser.parse_args()

input_file = args.input
output_file = args.output
threshold = args.threshold

# define functions
def calculate_overlap_length(hit_start, hit_end, start_curr, end_curr):
    overlap_len = max(0, min(hit_end, end_curr) - max(hit_start, start_curr) + 1)
    return overlap_len

# cluster matches with the same contig and domain and a specific overlap
contig_groups = defaultdict(dict)

with open(input_file, 'r') as inputfile:
    for line in inputfile:

        # cols: pro_id|domain_len|strand|contig_id|contig_len|hit_start|hit_end|score
        cols = line.strip().split('\t')
        domain_id, contig_id, hit_start, hit_end = cols[0], cols[3], int(cols[5]), int(cols[6])

        # domain：A/C/T/E/TE
        domain = domain_id.split('.')[0]

        # nucleotide length
        hit_len = hit_end - hit_start + 1

        # contig
        if contig_id not in contig_groups:
            contig_groups[contig_id] = {}
            count = 0
        
        # initial variable
        if 'domain_curr' not in locals():
           count += 1
        domain_curr = locals().get('domain_curr', domain)
        start_curr = locals().get('start_curr', hit_start)
        end_curr = locals().get('end_curr', hit_end)
        len_curr = locals().get('len_curr', end_curr - start_curr + 1)

        # calculate overlap length
        overlap_len = calculate_overlap_length(hit_start, hit_end, start_curr, end_curr)

        # overlap threshold
        if (overlap_len / len_curr) >= threshold and (overlap_len / hit_len) >= threshold and domain == domain_curr:
                contig_groups[contig_id].setdefault((count, domain_curr), []).append(line)
                start_curr = min(start_curr, hit_start)
                end_curr = max(end_curr, hit_end)
        else:
            domain_curr = domain
            start_curr = hit_start
            end_curr = hit_end
            len_curr = end_curr - start_curr + 1
            count += 1
            contig_groups[contig_id].setdefault((count, domain_curr), []).append(line)

# combine matches
domain_groups = defaultdict(dict)

for contig_id, domains in sorted(contig_groups.items()):
    for (count, domain), lines in sorted(domains.items()):
        domain_name_curr = ''
        domain_len_curr = 0
        strand_curr = ''
        contig_len_curr = 0
        start_min_curr = float('inf')
        end_max_curr = 0
        score_curr = 0
        len_curr = 0

        for line in lines:

            # cols: pro_id|domain_len|strand|contig_id|contig_len|hit_start|hit_end|score
            cols = line.strip().split('\t')
            domain_name = cols[0].rsplit('.',1)[0]
            domain_len, hit_strand, contig_len = int(cols[1]), cols[2], int(cols[4])
            hit_start, hit_end, hit_score = int(cols[5]), int(cols[6]), int(cols[7])

            # score
            score = hit_score / domain_len

            contig_len_curr = contig_len

            if strand_curr == '':
                strand_curr = hit_strand

            if strand_curr == hit_strand:
                start_min_curr = min(start_min_curr, hit_start)
                end_max_curr = max(end_max_curr, hit_end)
                len_curr = end_max_curr - start_min_curr + 1

                if score > score_curr:
                    domain_name_curr = domain_name
                    domain_len_curr = domain_len
                    score_curr = score
            else:
                print(f"Error: Strand mismatch for {input_file}, {count}, {domain}.")

        domain_groups[contig_id].setdefault((count, domain), []).append((domain_name_curr,strand_curr,contig_id,contig_len_curr,start_min_curr,end_max_curr))

with open(output_file, 'w') as outputfile:
    for contig_id, domains in sorted(domain_groups.items()):
        for (count, domain), elements in sorted(domains.items()):
            for element in elements:
                out = '\t'.join(map(str, element)) + '\n'
                outputfile.write(out)