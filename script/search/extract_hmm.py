import sys

# thresholds
thresholds = {
    "A": [1.0e-225, 750.0],
    "C": [1.0e-150, 590.0],
    "T": [1.0e-40, 135],
    "E": [1.0e-150, 600 ],
    "TE": [1.0e-100, 400],
}

hmmscan_out = sys.argv[1]
output_file = sys.argv[2]

def read_cols(line):
    cols = line.strip().split()
    target_name = cols[0]
    query_name = cols[2]
    seq_evalue = float(cols[4])
    seq_score = float(cols[5])
    domain_evalue = float(cols[7])
    domain_score = float(cols[8])
    domain = target_name.split('.')[0]
    contig = query_name.split(':')[0]
    hit_start = query_name.split(':')[1].split('-')[0]
    hit_end = query_name.split(':')[1].split('-')[1]
    frame_strand = query_name.split(':')[2]

    return(domain, target_name, query_name, contig, hit_start, hit_end, frame_strand, seq_evalue, seq_score, domain_evalue, domain_score)

with open(hmmscan_out, 'r') as scan_out, open(output_file, 'w') as outputfile:
    lines = [line for line in scan_out if not line.startswith('#')]

    for line in lines:
        domain, target_name, query_name, contig, hit_start, hit_end, frame_strand, seq_evalue, seq_score, domain_evalue, domain_score = read_cols(line)
        threshold_evalue, threshold_score = thresholds.get(domain, (None, None))

        if threshold_evalue is not None and threshold_score is not None:
            if seq_evalue <= threshold_evalue and seq_score >= threshold_score:
            # if seq_evalue <= threshold_evalue and seq_score >= threshold_score and domain_evalue <= threshold_evalue and domain_score >= threshold_score:
                outputfile.write(f"{target_name}\t{contig}\t{hit_start}\t{hit_end}\t{frame_strand}\n")