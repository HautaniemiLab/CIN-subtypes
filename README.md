# CIN-subtypes
Discovery of ovarian cancer subtypes based on chromosomal instability signatures

This repository contains scripts and certain data files used in the "Decoding the Genomic and Functional Landscape of Emerging Subtypes in Ovarian Cancer " manuscript by Micoli et al. (2025).
To use this repository, copy-number, breakpoint and structural variation data must be obtained using the tools mentioned from the [Hartwig Medical Foundation](https://pages.github.com/) repository.

Required packages and versions can be found in `requirements.txt`

### Features quantification

Features quantification from pre-existing modules can be done by running the R script `extract_features.R`. Input required: 
* sample_info.tsv: data frame with sample information. Specifically, it must have columns "sample", "patient" and "path". The "path" variable consist in the path directing to the output from PURPLE and LINX for each sample.
* models.rds: procided models for discretization of continuous variables (provided in the `data` folder)
* outout path

The resulting file `extraction.df` contains the quantified features for each sample and can be used for *de novo* signature extraction or signature attribution.

### Signatures
This step requires installation of [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) and [SigProfilerAssignment](https://github.com/AlexandrovLab/SigProfilerAssignment). Documentation for the tools in the respective repositories.

Signature extraction can be run with:
```
python run_sigProfilerExtractor.py 'path/to/feature_df.tsv' 'outhput/folder/path'
```

Signature attribution can be run with: 
```
python run_sigProfilerAssignment.py 'path/to/feature_df.tsv' 'path/to/Signatures_definition.txt' 'outhput/folder/path'
```
The input file `Signature_definition-txt` can be obtained by running SigProfilerExtractor or can be found in the `data` folder.

### Clustering
Signature activities are clustered using the R package [ConsensusClusterPlus](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html). Iterative clustering is performed by running:
```
Rscript iterative_clustering.R activities_file.txt sample_info_file.tsv output/folder/path/
```

Input:
* activities_file.txt: signature activities file from SigProfilerExtractor or SigProfilerAssignment
* sample_info_file.tsv: data frame with columns "sample" and "patient"
* output_folder

After that, perform visual evaluation of the output plots from ConsensusClusterPlus to decide at which number of clusters to cut the tree. Then, evaluate samples assignment to clusters in different runs and retrieve the prevailing cluster per sample. Use the script `cluster_analysis_and_selection.R`.
