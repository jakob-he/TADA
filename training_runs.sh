#!/usr/bin/env bash

#This file annotates TADs and CNVs with different sets of features and trains a classifier based on these features.
#For each setting a new directory is created. Single runs can be started using their ID. If no input is given all runs are executed.

#The tada environment needs to be active
#conda activate tada

#1. Run
#Basic binary features: GnomAD genes overlap, FANTOM5 enhancer overlap, TAD boundary breaking.
#Annotate TADs
mkdir -p classification_run_1
python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes enhancers -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed -o classification_run_1/annotated_TADs.p
python3 annotate_cnvs.py -t classification_run_1/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_variants_audano_gnomad_trimmed_filtered.bed -o classification_run_1/non_pathogenic_cnvs.p
python3 annotate_cnvs.py -t classification_run_1/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_and_number_matched.bed -o classification_run_1/pathogenic_cnvs.p
python3 classification_run.py -c1 classification_run_1/non_pathogenic_cnvs.p -c2 classification_run_1/pathogenic_cnvs.p -f basic_binary -cls lr -kw solver=liblinear -o classification_run_1 > classification_run_1/log.txt
