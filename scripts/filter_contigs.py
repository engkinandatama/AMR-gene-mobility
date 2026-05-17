#!/usr/bin/env python3
import sys
import os

def filter_fasta(input_file, output_file, min_len):
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist!")
        sys.exit(1)
        
    print(f"Filtering {input_file} for contigs >= {min_len} bp...")
    count_in = 0
    count_out = 0
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        header = None
        seq = []
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    count_in += 1
                    seq_str = ''.join(seq)
                    if len(seq_str) >= min_len:
                        fout.write(header + '\n' + seq_str + '\n')
                        count_out += 1
                header = line
                seq = []
            else:
                seq.append(line)
        # write the last one
        if header:
            count_in += 1
            seq_str = ''.join(seq)
            if len(seq_str) >= min_len:
                fout.write(header + '\n' + seq_str + '\n')
                count_out += 1
                
    print(f"Done. Kept {count_out} / {count_in} contigs.")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: filter_contigs.py <input.fa> <output.fa> <min_len>")
        sys.exit(1)
    filter_fasta(sys.argv[1], sys.argv[2], int(sys.argv[3]))
