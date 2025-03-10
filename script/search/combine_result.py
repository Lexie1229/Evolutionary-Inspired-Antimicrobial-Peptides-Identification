import os
import argparse

def load_hmm_result(hmm_file):
    if not os.path.isfile(hmm_file) or os.path.getsize(hmm_file) == 0:
        return {}
    hmm_dict = {}
    with open(hmm_file, 'r') as hmmfile:
        for line in hmmfile:
            if line.startswith(('>', '>>')):
                continue  # Skip lines that start with '>' or '>>'
            parts = line.strip().split('\t')
            key = (parts[1], parts[2], parts[3])
            hmm_dict[key] = parts[0]
    return hmm_dict

def combine_result(search_file, hmm_file, output_file):
    hmm_result = load_hmm_result(hmm_file)

    if not os.path.isfile(search_file) or os.path.getsize(search_file) == 0:
        open(output_file, 'w').close()  # Create an empty output file
        return

    with open(search_file, 'r') as searchfile, open(output_file, 'w') as outputfile:
        for line in searchfile:
            if line.startswith(('>', '>>')):
                outputfile.write(line)  # Obtain lines that start with '>' or '>>'
                continue
            parts = line.strip().split('\t')
            key = (parts[6], parts[8], parts[9])
            if key in hmm_result:
                combine_line = '\t'.join(parts) + f"\t{hmm_result[key]}"
            else:
                combine_line = '\t'.join(parts) + f"\tNA"
            #     combine_line = '\t'.join(parts) + f"\t{parts[3]}"
            outputfile.write(combine_line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge results from search results and hmm results.')
    parser.add_argument('-s', '--search_file', type=str, required=True, help='search results.')
    parser.add_argument('-m', '--hmm_file', type=str, required=True, help='hmm results')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='output files.')
    
    args = parser.parse_args()
    
    combine_result(args.search_file, args.hmm_file, args.output_file)

if __name__ == '__main__':
    main()
