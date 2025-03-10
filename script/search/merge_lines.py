# python merge_lines.py -i input.tsv -o output.tsv

import argparse

def merge_lines(input_file, output_file):
    merged_dict = {}
    
    with open(input_file, 'r') as inputfile:
        lines = inputfile.readlines()
    
    for line in lines:
        parts = line.strip().split('\t')
        key = tuple(parts[1:])
        
        if key in merged_dict:
            merged_dict[key].append(parts[0])
        else:
            merged_dict[key] = [parts[0]]
    
    with open(output_file, 'w', encoding='utf-8') as f:
        for key, values in merged_dict.items():
            merged_line = "{}\t{}".format(':'.join(values), '\t'.join(key))
            f.write(merged_line + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge lines in a file based on column contents.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='The input file path')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='The output file path')
    
    args = parser.parse_args()
    
    merge_lines(args.input_file, args.output_file)

if __name__ == '__main__':
    main()