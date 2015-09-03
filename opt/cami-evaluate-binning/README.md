# Evaluation Framework for Taxonomic Binning by J. Dröge

## Procedure
Both the labels and the predictions provided as NCBI taxon IDs are converted to tab-separated lineages. In these tables, the first column contains the sequence IDs and each following column stands for taxonomic rank and columns are in ascending order, starting ranging a low-level rank like species to a high-level rank like phylum. This conversion requires information from the NCBI taxonomy but all following steps only depend on string comparisons. We refer to this format as the RACOL (rank column) format in the following.

For hierarchical classifications like multiple rank taxa, there are two kinds of evaluation measures. First, there are measures which only consider classifications at a single taxonomic rank (using e.g. only the species column in the two generated input files). The best-known amoung these are simple specificity/sensitivity and precision/recall. Secondly, there are measures which consider classifications at multiple ranks (using more than one input column) simultaneously. The first kind of measures are typically calculated using a confusion matrix for the particular rank, which the provided scripts are doing. Confusion matrices are generated from two RACOL input files and concatenated in a text output. A confusion matrix in text format can be parsed into a Python object which implements the different performance measures for single-rank evaluation. The second kind of measures are directly computed from two tab-separated input files. If the binning predictions are not hierarchical, e.g. not given as taxon IDs, then the two-column prediction file is already in valid RACOL format and can be used to generate a single confusion matrix.

## Files

### classevaltools.py
A python library for evaluation.

### taxonomyncbi.py
A python library for taxonomy access.

### tax2racol
A Python script which takes a tab-separated two-column file where the first columns contains the sequence ID and the second an NCBI taxon ID. The output will be in RACOL format where the first column is the sequence ID and the following columns stand for taxonomic ranks in ascending order and contain the taxon names. In addition to the input (provided as standard input), the script allows to specify for which ranks to generate columns and also requires the user to provide an NCBI taxonomy which must be in SQLite-BioSQL format. These files can be constructed from the raw NCBI taxonomy files (names.dmp, nodes.dmp) by a provided script (available very soon). If this seems a too complicated dependence, this script could easily be replaced by a more lightweight version.

### fasta-seqlen
This is an AWK script to calculate the length of FASTA sequence entries. The FASTA file is streamed via the standard input and the sequence ID and length are printed on the standard output. If piped to a file, this output is a proper weights file for the confusion_matrix script.

### confusionmatrix
This Python script takes two RACOL files as input (representing row and column classes) and prints the confusion matrices on the standard output. Items (classified sequences or other classified objects) can be weighted (usually by sequence length) by providing a two-column tab-separated weights file.

### cmat2*
These Python scripts parse confusion matrices in text form and output statics or plots. They just use the functionality which is implemented in the Python objects.

### count-depth_true_false_unknown
This Python script calculates the amount of false, true and unknown data at each rank (depth) but without mapping lower to higher taxa. More colloquially, this evaluation method doesn't forgive false predictions and is very sensitive to optimistic taxonomic predictions. Because this kind of information spans over multiple ranks, it is calculated directly from the input RACOL files.

## Application

### Generate files in RACOL format
```bash
tax2racol -t taxonomy.sqlite -ranks genus,family,order < labels.tax > labels.racol
tax2racol -t taxonomy.sqlite -ranks genus,family,order < predictions.tax > predictions.racol
```

### Generate weights (sequence length) file
```bash
fasta-seqlen.awk < predictions.fna > predictions.seqlen
```

### Generate confusion matrices
```bash
confusion-matrix --rows labels.racol --columns predictions.racol --weights predictions.seqlen --matrix-form quadratic --allow-missing-columns > predictions.cmat
```

### Display basic statistics
```bash
cmat2stat < preditions.cmat
```

## Measure definitions and plots
For understanding the measure currently implemented in this framework, please refer to http://arxiv.org/abs/1404.1029

## Feedback
Please send feedback and questions to johannes.droege@uni-duesseldorf.de

