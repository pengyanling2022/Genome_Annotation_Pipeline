import re
import argparse

def filter_gff3(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        buffer = []
        features = {'gene': False, 'CDS': False, 'exon': False, 'mRNA': False}

        for line in infile:
            if line.startswith("###"):
                # Only write out the block if all features are present
                if all(features.values()):
                    outfile.write("###\n")
                    outfile.writelines(buffer)
                # Reset for the next block
                buffer = []
                features = {key: False for key in features}
            else:
                buffer.append(line)
                for feature in features:
                    if re.search(rf'\b{feature}\b', line):
                        features[feature] = True
        
        # Write the last block if all features are present
        if all(features.values()):
            outfile.write("###\n")
            outfile.writelines(buffer)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter GFF3 files.')
    parser.add_argument('input_gff3', type=str, help='Path to the input GFF3 file')
    parser.add_argument('output_gff3', type=str, help='Path to the output filtered GFF3 file')
    
    args = parser.parse_args()
    
    filter_gff3(args.input_gff3, args.output_gff3)