# python extract_score_mean.py -i input.tsv -o output.tsv

import argparse
import pandas as pd
import re

def rename_column(name):
    name = re.sub(r'\.cif:A$', '', name)
    parts = name.split('_')
    new_parts = []
    for index, value in enumerate(parts):
        if index < 4:
            new_parts.append(value.capitalize())
        elif index == 4:
            new_parts.append(value.lower())
        else:
            new_parts.append(value.upper())
    new_name = '_'.join(new_parts)
    new_name = new_name.replace('_Dcl_', '_DCL_')
    new_name = new_name.replace('_Lcl_', '_LCL_')
    new_name = new_name.replace('_Tl_', '_TL_')
    new_name = new_name.replace('_Td_', '_TD_')
    return new_name

def calculate_mean_of_score(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', header=0)
    df_extracted = df.iloc[:, :4].copy()

    # rename
    df_extracted.iloc[:, 0] = df_extracted.iloc[:, 0].apply(rename_column)
    df_extracted.iloc[:, 1] = df_extracted.iloc[:, 1].apply(rename_column)
    # TM score average
    df_extracted['TM_mean'] = df_extracted.iloc[:, 2:4].mean(axis=1).round(4)
    # distance
    df_extracted['Distance'] = (1 - df_extracted['TM_mean']).round(4)

    df_extracted.to_csv(output_file, sep='\t', index=False, header=True)

def main():
    parser = argparse.ArgumentParser(description='Calculate the mean of the TM score.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file.')
    args = parser.parse_args()

    calculate_mean_of_score(args.input, args.output)

if __name__ == "__main__":
    main()
