import pandas as pd
import os
from SigProfilerAssignment import Analyzer as Analyze
import sys

#### Signature attribution usign SigProfilerAssignment
#
# Run as: python run_sigProfilerAssignment.py "/path/feature_matrix.tsv" "path/to/De-Novo_Signatures.txt" "output/path/"

def run_sigAssignment(input_path, signatures_path, output_path):    
    Analyze.denovo_fit(
        samples= input_path,
        output= output_path,
        signatures= signatures_path,
        signature_database=None,
        nnls_add_penalty=0.05,
        nnls_remove_penalty=0.01,
        initial_remove_penalty=0.05,
        genome_build="GRCh38",
        cosmic_version=3.3,
        new_signature_thresh_hold=0.8,
        make_plots=True,
        collapse_to_SBS96=True,
        connected_sigs=True,
        verbose=False,
        devopts=None,
        exome=False,
        input_type="matrix")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_path> <signatures_path> <output_path>")
        sys.exit(1)

    input_path = sys.argv[1]
    signatures_path = sys.argv[2]
    output_path = sys.argv[3]

    run_sigAssignment(input_path, signatures_path, output_path)