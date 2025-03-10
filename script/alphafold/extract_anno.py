# python extract_anno.py -i input.tsv -o output.tsv

import argparse
import pandas as pd

def extract_anno(name):
    parts = name.split('_')
    domain = parts[0]
    function = parts[1]
    site = parts[2]
    strain = '_'.join(parts[3:])
    return name, strain, domain, function, site

def read_and_write(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', header=None)
    
    results = []
    for index, row in df.iterrows():
        name = row[0]
        anno = extract_anno(name)
        results.append(anno)
    
    df_result = pd.DataFrame(results, columns=['Name', 'Strain', 'Domain', 'Function', 'Site'])
    df_result.to_csv(output_file, sep='\t', index=False, header=True)

def main():
    parser = argparse.ArgumentParser(description='Extract the anno info.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file.')
    args = parser.parse_args()

    read_and_write(args.input, args.output)

if __name__ == "__main__":
    main()