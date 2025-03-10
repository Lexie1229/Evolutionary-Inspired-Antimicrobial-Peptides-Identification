import sys
from collections import defaultdict

# distance thresholds
thresholds = {
        ('A', 'C'): [(450, 560),(1230, 1690)],
        ('A', 'T'): [(310, 340)],
        ('C', 'E'): [(490, 510)],
        ('C', 'T'): [(40, 60)],
        ('C', 'TE'): [(3550, 4040)],
        ('T', 'E'): [(20, 50)],
        ('T', 'TE'): [(70, 80)]
    }

# define functions
def check_distance(prev_domain, domain, distance):
    if (prev_domain, domain) in thresholds:
        for range_min, range_max in thresholds[(prev_domain, domain)]:
            if range_min <= distance <= range_max:
                return True
    elif (domain, prev_domain) in thresholds:
        for range_min, range_max in thresholds[(domain, prev_domain)]:
            if range_min <= distance <= range_max:
                return True
    return False


input_file = sys.argv[1]
output_file = sys.argv[2]

cluster_groups = defaultdict(dict)

with open(input_file, 'r') as inputfile:
    for line in inputfile:
        if line.startswith('>'):
            
            # >contig_id
            contig_id = line.strip()[1:]

            cluster_count = 0
            domain_count = 0
            start_curr = 0
            end_curr = 0
            distance_curr = 0
            domain_curr = '' 
            cluster_groups.setdefault(contig_id, {})

        else:
            # cols: domain_count|domain|domain_name|hit_strand|contig_id|contig_len|hit_start|hit_end
            cols = line.strip().split('\t')
            domain = cols[1]
            hit_start = int(cols[6])
            hit_end = int(cols[7])

            if end_curr == 0:
                end_curr = hit_end
                domain_curr = domain
                cluster_count += 1
                domain_count += 1
                cluster_groups[contig_id].setdefault(cluster_count, {})[domain_count]= line
            else:
                distance_curr = hit_start - end_curr

                if check_distance(domain_curr, domain, distance_curr):
                    end_curr = hit_end
                    domain_curr = domain
                    domain_count += 1
                    cluster_groups[contig_id].setdefault(cluster_count, {})[domain_count]= line
                else:
                    end_curr = hit_end
                    domain_curr = domain
                    cluster_count += 1
                    domain_count = 1
                    cluster_groups[contig_id].setdefault(cluster_count, {})[domain_count]= line

with open(output_file, 'w') as outputfile:
    for contig_id, clusters in sorted(cluster_groups.items()):
        outputfile.write(f">{contig_id}\n")
        for cluster_num, domains in sorted(clusters.items()):
            outputfile.write(f">>cluster.{cluster_num}\n")
            for domain_num, lines in sorted(domains.items()):
                outputfile.write(f"c.{cluster_num}\td.{domain_num}\t")
                for line in lines:
                    outputfile.write(line)
