#!/usr/bin/env bash

#This file annotates TADs and CNVs with different sets of features and trains a classifier based on these features.
#For each setting a new directory is created. Single runs can be started using their ID. If no input is given all runs are executed.

#The tada environment needs to be active
#conda activate tada

#1. Run
#Basic binary features: GnomAD genes overlap, FANTOM5 enhancer overlap, TAD boundary breaking.
mkdir -p Basic_binary
python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes enhancers -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed -o Basic_binary/annotated_TADs.p
python3 annotate_cnvs.py -t Basic_binary/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Basic_binary/non_pathogenic_cnvs.p
python3 annotate_cnvs.py -t Basic_binary/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Basic_binary/pathogenic_cnvs.p
python3 classification_run.py -fs -c1 Basic_binary/non_pathogenic_cnvs.p -c2 Basic_binary/pathogenic_cnvs.p -f basic_binary -cls lr -kw solver="liblinear" -o Basic_binary > Basic_binary/log_lr.txt
python3 classification_run.py -kw n_estimators=10 -c1 Basic_binary/non_pathogenic_cnvs.p -c2 Basic_binary/pathogenic_cnvs.p -f basic_binary -cls rf -o Basic_binary > Basic_binary/log_rf.txt

#2. Run
#Extended binary features: GnomAD genes overlap, DDG2P genes, FANTOM5 enhancer overlap, TAD boundary breaking, CTCF binding site overlap, TADs with high pLI genes (=1), TADs with highly conserved enhancer (>=0.9)
mkdir -p Extended_binary
python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes DDG2P enhancers CTCF -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/GENES/DDG2P_genes.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed /project/GENTLE_SCRATCH/DATA/CTCF/ctcf_peaks_fibroblasts.bed -o Extended_binary/annotated_TADs.p
python3 annotate_cnvs.py -t Extended_binary/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Extended_binary/non_pathogenic_cnvs.p
python3 annotate_cnvs.py -t Extended_binary/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Extended_binary/pathogenic_cnvs.p
python3 classification_run.py -fs -c1 Extended_binary/non_pathogenic_cnvs.p -c2 Extended_binary/pathogenic_cnvs.p -f extended_binary -cls lr -kw solver="liblinear" -o Extended_binary > Extended_binary/log_lr.txt
python3 classification_run.py -kw n_estimators=10 -c1 Extended_binary/non_pathogenic_cnvs.p -c2 Extended_binary/pathogenic_cnvs.p -f extended_binary -cls rf -o Extended_binary > Extended_binary/log_rf.txt

#3. Run
#Basic continuous features: GnomAD genes distance, FANTOM5 enhancer distance, TAD boundariy distance
mkdir -p Basic_continuous
python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes enhancers -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed -o Basic_continuous/annotated_TADs.p
python3 annotate_cnvs.py -t Basic_continuous/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Basic_continuous/non_pathogenic_cnvs.p
python3 annotate_cnvs.py -t Basic_continuous/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Basic_continuous/pathogenic_cnvs.p
python3 classification_run.py -fs -c1 Basic_continuous/non_pathogenic_cnvs.p -c2 Basic_continuous/pathogenic_cnvs.p -f basic_continuous -cls lr -kw solver="liblinear" -o Basic_continuous > Basic_continuous/log_lr.txt
python3 classification_run.py -kw n_estimators=10 -c1 Basic_continuous/non_pathogenic_cnvs.p -c2 Basic_continuous/pathogenic_cnvs.p -f basic_continuous -cls rf -o Basic_continuous > Basic_continuous/log_rf.txt

#4. Run
#Extended continuous features: GnomAD genes distance, FANTOM5 enhancer distance, TAD boundariy distance
mkdir -p Extended_continuous
python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes DDG2P enhancers CTCF -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf_HI.bed /project/GENTLE_SCRATCH/DATA/GENES/DDG2P_genes.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed /project/GENTLE_SCRATCH/DATA/CTCF/ctcf_peaks_fibroblasts.bed -o Extended_continuous/annotated_TADs.p
python3 annotate_cnvs.py -t Extended_continuous/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Extended_continuous/non_pathogenic_cnvs.p
python3 annotate_cnvs.py -t Extended_continuous/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Extended_continuous/pathogenic_cnvs.p
python3 classification_run.py -fs -c1 Extended_continuous/non_pathogenic_cnvs.p -c2 Extended_continuous/pathogenic_cnvs.p -f extended_continuous -cls lr -kw solver="liblinear" -o Extended_continuous > Extended_continuous/log_lr.txt
python3 classification_run.py -kw n_estimators=10 -c1 Extended_continuous/non_pathogenic_cnvs.p -c2 Extended_continuous/pathogenic_cnvs.p -f extended_continuous -cls rf -o Extended_continuous > Extended_continuous/log_rf.txt
python3 classification_run.py -kw n_estimators=10 -c1 Extended_continuous/non_pathogenic_cnvs.p -c2 Extended_continuous/pathogenic_cnvs.p -f extended_continuous -cls brf -o Extended_continuous > Extended_continuous/log_brf.txt


# #5. Run
# #This Run includes the UK biobank variants
# #Basic binary features: GnomAD genes overlap, FANTOM5 enhancer overlap, TAD boundary breaking.
# mkdir -p Basic_binary_biobank
# python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes enhancers -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed -o Basic_binary_biobank/annotated_TADs.p
# python3 annotate_cnvs.py -t Basic_binary_biobank/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Basic_binary_biobank/non_pathogenic_cnvs.p
# python3 annotate_cnvs.py -t Basic_binary_biobank/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Basic_binary_biobank/pathogenic_cnvs.p
# python3 classification_run.py -c1 Basic_binary_biobank/non_pathogenic_cnvs.p -c2 Basic_binary_biobank/pathogenic_cnvs.p -f basic_binary -cls lr -kw solver="liblinear" -o Basic_binary_biobank > Basic_binary_biobank/log_lr.txt
# python3 classification_run.py -c1 Basic_binary_biobank/non_pathogenic_cnvs.p -c2 Basic_binary_biobank/pathogenic_cnvs.p -f basic_binary -cls rf -o Basic_binary_biobank > Basic_binary_biobank/log_rf.txt
#
# #6. Run
# #This Run includes the UK biobank variants
# #Extended binary features: GnomAD genes overlap, DDG2P genes, FANTOM5 enhancer overlap, TAD boundary breaking, CTCF binding site overlap, TADs with high pLI genes (=1), TADs with highly conserved enhancer (>=0.9)
# mkdir -p Extended_binary_biobank
# python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes DDG2P enhancers CTCF -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/GENES/DDG2P_genes.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed /project/GENTLE_SCRATCH/DATA/CTCF/ctcf_peaks_fibroblasts.bed -o Extended_binary_biobank/annotated_TADs.p
# python3 annotate_cnvs.py -t Extended_binary_biobank/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Extended_binary_biobank/non_pathogenic_cnvs.p
# python3 annotate_cnvs.py -t Extended_binary_biobank/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Extended_binary_biobank/pathogenic_cnvs.p
# python3 classification_run.py -c1 Extended_binary_biobank/non_pathogenic_cnvs.p -c2 Extended_binary_biobank/pathogenic_cnvs.p -f extended_binary -cls lr -kw solver="liblinear" -o Extended_binary_biobank > Extended_binary_biobank/log_lr.txt
# python3 classification_run.py -c1 Extended_binary_biobank/non_pathogenic_cnvs.p -c2 Extended_binary_biobank/pathogenic_cnvs.p -f extended_binary -cls rf -o Extended_binary_biobank > Extended_binary_biobank/log_rf.txt

# #7. Run
# #This Run includes the UK biobank variants
# #Basic continuous features: GnomAD genes distance, FANTOM5 enhancer distance, TAD boundariy distance
# mkdir -p Basic_continuous_biobank
# python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes enhancers -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed -o Basic_continuous_biobank/annotated_TADs.p
# python3 annotate_cnvs.py -t Basic_continuous_biobank/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Basic_continuous_biobank/non_pathogenic_cnvs.p
# python3 annotate_cnvs.py -t Basic_continuous_biobank/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Basic_continuous_biobank/pathogenic_cnvs.p
# python3 classification_run.py -c1 Basic_continuous_biobank/non_pathogenic_cnvs.p -c2 Basic_continuous_biobank/pathogenic_cnvs.p -f basic_continuous -cls lr -kw solver="liblinear" -o Basic_continuous_biobank > Basic_continuous_biobank/log_lr.txt
# python3 classification_run.py -c1 Basic_continuous_biobank/non_pathogenic_cnvs.p -c2 Basic_continuous_biobank/pathogenic_cnvs.p -f basic_continuous -cls rf -o Basic_continuous_biobank > Basic_continuous_biobank/log_rf.txt
#
# #7. Run
# #This Run includes the UK biobank variants
# #Extended continuous features: GnomAD genes distance, DDG2P genes distance, FANTOM5 enhancer distance, TAD boundary distance, CTCF binding site distance, pLi of the closest gene, enhancer conservation of the closest gene, haploinsufficiency of the closest gene
# mkdir -p Extended_continuous_biobank
# python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes DDG2P enhancers CTCF -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf_HI.bed /project/GENTLE_SCRATCH/DATA/GENES/DDG2P_genes.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed /project/GENTLE_SCRATCH/DATA/CTCF/ctcf_peaks_fibroblasts.bed -o Extended_continuous_biobank/annotated_TADs.p
# python3 annotate_cnvs.py -t Extended_continuous_biobank/annotated_TADs.p -c /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed -o Extended_continuous_biobank/non_pathogenic_cnvs.p
# python3 annotate_cnvs.py -t Extended_continuous_biobank/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o Extended_continuous_biobank/pathogenic_cnvs.p
# python3 classification_run.py -c1 Extended_continuous_biobank/non_pathogenic_cnvs.p -c2 Extended_continuous_biobank/pathogenic_cnvs.p -f extended_continuous -cls lr -kw solver="liblinear" -o Extended_continuous_biobank > Extended_continuous_biobank/log_lr.txt
# python3 classification_run.py -c1 Extended_continuous_biobank/non_pathogenic_cnvs.p -c2 Extended_continuous_biobank/pathogenic_cnvs.p -f extended_continuous -cls rf -o Extended_continuous_biobank > Extended_continuous_biobank/log_rf.txt

#8.Run
#This run tests multiple allele frequency cutoffs for the UK Biobank data and produces a plot of the average precision means
# mkdir -p biobank_AF_cutoff_comparision
# python3 annotate_tads.py -t /project/GENTLE_SCRATCH/DATA/TAD_BOUNDARIES/Dali/Dixon/TADcalls_TopDom_Dixon_500M_50000.bed -an genes DDG2P enhancers CTCF -af /project/GENTLE_SCRATCH/DATA/GENES/gnomad_genes_pli_loeuf_HI.bed /project/GENTLE_SCRATCH/DATA/GENES/DDG2P_genes.bed /project/GENTLE_SCRATCH/DATA/ENHANCERS/OTHER/FANTOM5/fantom5_enhancer_phastcon_average.bed /project/GENTLE_SCRATCH/DATA/CTCF/ctcf_peaks_fibroblasts.bed -o biobank_AF_cutoff_comparision/annotated_TADs.p
# python3 annotate_cnvs.py -t biobank_AF_cutoff_comparision/annotated_TADs.p -c /project/GENTLE_BKP/DECIPHER_CNV/mundlos_and_decipher_size_matched_biobank_max_500kb.bed -o biobank_AF_cutoff_comparision/pathogenic_cnvs.p
# for I in 1 2 3 4 5 6 7 8 9 10
# do
# awk -F "\t" "{if((\$7>=$I)&&(\$3-\$2>=50)){printf\"%s\t%s\t%s\n\",\$1,\$2,\$3}}" /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/BIOBANK/biobank_DEL.bed > "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/BIOBANK/biobank_DEL_AF_$I.bed"
# cat "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/BIOBANK/biobank_DEL_AF_$I.bed" /project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_gnomad_audano_biobank_size_matched_max_500kb.bed > "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_AF_$I.bed"
# python3 annotate_cnvs.py -t biobank_AF_cutoff_comparision/annotated_TADs.p -c "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_AF_$I.bed" -o biobank_AF_cutoff_comparision/non_pathogenic_cnvs.p
# python3 classification_run.py -c1 biobank_AF_cutoff_comparision/non_pathogenic_cnvs.p -c2 biobank_AF_cutoff_comparision/pathogenic_cnvs.p -f extended_continuous -cls rf -o biobank_AF_cutoff_comparision > "biobank_AF_cutoff_comparision/log_rf_AF_$I.txt"
# rm "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/BIOBANK/biobank_DEL_AF_$I.bed" "/project/GENTLE_SCRATCH/DATA/REFERENCE_STRUCTURAL_VARIANTS/shared_AF_$I.bed" biobank_AF_cutoff_comparision/non_pathogenic_cnvs.p
# done
