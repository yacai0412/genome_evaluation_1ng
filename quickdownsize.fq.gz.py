import sys
import gzip
import random
from Bio import SeqIO

input_file = sys.argv[1]
target_bases = int(sys.argv[2])
output_file = input_file.split("/")[-1].replace("fastq.gz", f"{target_bases}.fastq")

selected = []
total_bases = 0

with gzip.open(input_file, "rt") as handle:
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(handle, "fastq"):
            if random.random() < 0.1:  # random subsampling (adjust as needed)
                seq_len = len(record.seq)
                if total_bases + seq_len > target_bases:
                    break
                SeqIO.write(record, out_handle, "fastq")
                total_bases += seq_len

print(f"Done! Wrote reads with {total_bases:,} bases.")
