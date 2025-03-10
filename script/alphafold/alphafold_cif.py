# python alphafold_cif.py -i input_dir -o output_dir

import argparse
import os
import shutil

def extract_model_0_cif(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith("_model_0.cif"):
                parent_dir = os.path.basename(root)
                new_filename = f"{parent_dir}.cif"
                src_file = os.path.join(root, file)
                dst_file = os.path.join(output_dir, new_filename)
                shutil.copy(src_file, dst_file)

def main():
    parser = argparse.ArgumentParser(description='Extract model_0.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input path')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output path')
    args = parser.parse_args()
    extract_model_0_cif(args.input, args.output)

if __name__ == '__main__':
    main()