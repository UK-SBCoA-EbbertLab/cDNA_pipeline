#!/usr/bin/env python3

import sys
from math import log


def errs_tab(n):

    """Generate list of error rates for qualities less than equal than n."""
    return [10**(q / -10) for q in range(n + 1)]



def calculate_mean_quality(quals, qround=False, tab=errs_tab(128)):
    """Calculate average basecall quality of a read.
    Receive the ascii quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[ord(q) - 33] for q in quals]) / len(quals), 10)

        if qround:
            return round(mq)
        else:
            return mq
    else:
        return 0.0 
    

def filter_fastq(input_file, threshold, output_file):

    num_pass_thresh = 0

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            if calculate_mean_quality(quality) >= threshold:
                outfile.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
                num_pass_thresh = num_pass_thresh + 1
    
    return num_pass_thresh

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_fastq.py <input_file> <quality_threshold> <output_file>")
        sys.exit(1)

    input_fastq = sys.argv[1]
    quality_threshold = float(sys.argv[2])
    output_fastq = sys.argv[3]

    num_pass_thresh = filter_fastq(input_fastq, quality_threshold, output_fastq)
    print(num_pass_thresh)
