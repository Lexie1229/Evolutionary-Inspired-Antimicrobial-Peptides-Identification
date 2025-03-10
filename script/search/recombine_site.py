# python -t 0.95 -i input.tsv -o output.tsv

import argparse
from collections import defaultdict

# get options
parser = argparse.ArgumentParser(description="Combine matches at the same location from different process within a specific overlap threshold.")
parser.add_argument("-i", "--input", help="Path to the input file", required=True)
parser.add_argument("-o", "--output", help="Path to the output file", required=True)
parser.add_argument("-t", "--threshold", type=float, default=0.95, help="Overlap threshold (default: 0.95)")
args = parser.parse_args()

input_file = args.input
output_file = args.output
threshold = args.threshold

# define functions
def parse_domain_name(domain_name):
    domain = domain_name.split('.')[0]
    attr1 = domain_name.split('.')[1] # function
    attr2 = domain_name.split('.')[2] # site or substrate
    return (domain, attr1, attr2)

def calculate_overlap_length(hit_start, hit_end, start_curr, end_curr):
    overlap_len = max(0, min(hit_end, end_curr) - max(hit_start, start_curr) + 1)
    return overlap_len

# cluster matches with the same contig and domain and a specific overlap
contig_groups = defaultdict(dict)

with open(input_file, 'r') as inputfile:
    for line in inputfile:

        # cols: domain_name|hit_strand|contig_id|contig_len|hit_start|hit_end
        cols = line.strip().split('\t')
        domain_name, contig_id, hit_start, hit_end = cols[0], cols[2], int(cols[4]), int(cols[5])

        # domain：A/C/T/E/TE
        domain = domain_name.split('.')[0]

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

# combine sites
domain_groups = defaultdict(dict)

for contig_id, domains in sorted(contig_groups.items()):
    domain_groups[contig_id] = defaultdict(dict)
    for (count, domain), lines in sorted(domains.items()):
        domain_name_curr = ''
        strand_curr = ''
        contig_len_curr = 0
        start_min_curr = float('inf')
        end_max_curr = 0

        for line in lines:

            # cols: domain_name|strand|contig_id|contig_len|hit_start|hit_end
            cols = line.strip().split('\t')
            domain_name, hit_strand = cols[0], cols[1]
            contig_len, hit_start, hit_end = int(cols[3]), int(cols[4]), int(cols[5])

            contig_len_curr = contig_len

            if strand_curr == '':
                strand_curr = hit_strand

            if strand_curr == hit_strand:
                start_min_curr = min(start_min_curr, hit_start)
                end_max_curr = max(end_max_curr, hit_end)

                if domain_name_curr == '':
                    domain_name_curr = domain_name
                else:
                    if domain_name != domain_name_curr:
                        (domain, attr1, attr2) = parse_domain_name(domain_name)
                        (domain_curr, attr1_curr, attr2_curr) = parse_domain_name(domain_name_curr)

                        if attr1 == attr1_curr:
                            domain_name_curr = domain + '.' + attr1 + '.' + attr2 + ':' + attr2_curr
                        else:
                            if attr2 == attr2_curr:
                                domain_name_curr = domain + '.' + attr1 + ':' + attr1_curr + '.' + attr2
                            else:
                                domain_name_curr = domain + '.' + attr1 + ':' + attr1_curr + '.' + attr2 + ':' + attr2_curr
                    else:
                        domain_name_curr = domain_name
            else:
                print(f"Error: Strand mismatch for {input_file}, {line}")

        domain_groups[contig_id].setdefault((count, domain), []).append((domain_name_curr,strand_curr,contig_id,contig_len_curr,start_min_curr,end_max_curr))

with open(output_file, 'w') as outputfile:
    for contig_id, domains in sorted(domain_groups.items()):
        outputfile.write(f">{contig_id}\n")
        for (count, domain), elements in sorted(domains.items()):
            outputfile.write(f"{count}\t{domain}\t")
            for element in elements:
                out = '\t'.join(map(str, element)) + '\n'
                outputfile.write(out)

# domain_name
# A.L.Dab
# C.DCL.C4
# T.TL.T4_T3
# E.E.E3
# TE.TE.TE10