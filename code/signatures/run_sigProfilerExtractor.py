import pandas as pd
import random
from SigProfilerExtractor import sigpro as sig
import sys

#### Signature extraction usign SigProfilerExtractor
#
# Run as: python run_sigProfilerExtractor.py "/path/feature_matrix.tsv" "output/path/" 

def run_sigProfiler(input_path, output_path): 

    sig.sigProfilerExtractor('matrix', 
                            output_path, 
                            input_path,
                            opportunity_genome = "GRCh38", 
                            context_type = 'CNV48')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_path> <output_path>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    run_sigProfiler(input_path, output_path)