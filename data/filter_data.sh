#!usr/bin/bash
#This file filters contains the filter processes for the data sources described in data_sources.txt

## CNVs

#G nomad CNVs

# filter for PASS variants and exlcude header
zcat gnomad_v2_sv.sites.bed.gz | awk -F "\t" '{if(($1 ~ /^#.*/)||($7=="PASS")) {print}}' > gnomad_v2_sv.sites.PASS.bed

# filter for CNVs with a size greater than 50bp and a allele frequency greater than 0.1
awk -F "\t" '{if(($5=="DUP" || $5=="DEL") && ($31>0.1) && ($6>50)) {print}}' gnomad_v2_sv.sites.PASS.bed > gnomad_v2_sv.sites.PASS.frequent.bed

# filter for CNVs with a size greater than 50bp and a allele frequency less than 0.1
awk -F "\t" '{if(($5=="DUP" || $5=="DEL") && ($31>0.1) && ($6>50)) {print}}' gnomad_v2_sv.sites.PASS.bed > gnomad_v2_sv.sites.PASS.rare.bed

# Fix positions in scientific format to decimal numbers
awk -F "\t" '{if($3~/^[0-9]*.*[eE][+-][0-9]*$/) {$3=sprintf("%.0f", $3)};};1' gnomad_v2_sv.sites.PASS.frequent.bed > gnomad_v2_sv.sites.PASS.frequent.fixedlines.bed
awk -F "\t" '{if($3~/^[0-9]*.*[eE][+-][0-9]*$/) {$3=sprintf("%.0f", $3)};};1' gnomad_v2_sv.sites.PASS.rare.bed > gnomad_v2_sv.sites.PASS.rare.fixedlines.bed

# Replace whitespaces with tabs in malformated lines
awk -v OFS="\t" '$1=$1' gnomad_v2_sv.sites.PASS.frequent.bed > gnomad_v2_sv.sites.PASS.frequent.fixedlines.formatted.bed
awk -v OFS="\t" '$1=$1' gnomad_v2_sv.sites.PASS.rare.bed > gnomad_v2_sv.sites.PASS.rare.fixedlines.formatted.bed

#Eicher et al. Variants

#filter for MAJOR and SHARED CNVs larger than 50bp
zcat EEE_SV-Pop_1.ALL.sites.20181204.bed.gz | awk -F "\t" '{if(($19=="MAJOR" || $19=="SHARED") && ($5=="DEL") && ($3-$2>=50)) {print}}' > EEE_SV-Pop_1.ALL.sites.20181204.shared.major.bed

# LiftOver variants to hg19 using the USCS LiftOver tool available at: https://genome-store.ucsc.edu/ (Free for non-comercial users)
# Chain files can be obainted from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# The path to the liftOver has to be adapted to the user-specific path
./liftOver -bedPlus=3 EEE_SV-Pop_1.ALL.sites.20181204.shared.major.bed hg38ToHg19.over.chain.gz EEE_SV-Pop_1.ALL.sites.20181204.shared.major.hg19.bed EEE_SV-Pop_1.ALL.sites.20181204.shared.major.liftoverfailed.bed

# Combinedfrequenct variants
cat EEE_SV-Pop_1.ALL.sites.20181204.shared.major.hg19.bed gnomad_v2_sv.sites.PASS.frequent.fixedlines.formatted.bed > combined.frequenct.Eichler.GnomAD.bed

## Genes

# Gnomad Genes

# change suffix to .gz
mv gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz gnomad.v2.1.1.lof_metrics.by_gene.txt.gz

# filter for gene position, gene name, pLI and LOEUF
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.gz | awk -F "\t" '{if(NR==1){printf "#CHROM\tSTART\tEND\tGENE\tpLI\tLOEUF\n"}else{printf "chr%s\t%s\t%s\t%s\t%s\t%s\n",$(NF-2),$(NF-1),$NF,$1,$21,$30'}} > GnomAD.genes.LOEUF.pLI.bed

# DDG2P Genes

# filter for gene names
zcat DDG2P.csv.gz | awk -F "," '{printf "%s\n",$1}' > DDG2P.Genes.bed

# Add positions to genes
python3 add_positions_to_DDG2P_genes.py

# Add haploinsuffiency scores to GnomAD genes
gunzip HI_Predictions_Version3.bed.gz
python3 add_HI_to_gnomad_genes.py

## Enhancer

# FANTOM5 Enhancer

# trim enhancer to chromosome, start and end position
zcat human_permissive_enhancers_phase_1_and_2.bed.gz | awk -F "\t" '{printf "%s\t%s\t%s\n",$1,$2,$3}' > FANTOM5_enhancer.bed

# Add PhastCon and PhyloP scores (requires the 100way bigwig files)
# python3 examples/conservation_annotation.py -b enhancers.bed -c 100way.bw -o annotated_enhancers.bed

## TAD boundaries

# decrompress the supplement
unzip gkx145_Supp.zip

# cd into the Supplement folder
cd Supplementary

# convert the Dixon TAD boundary file into Unix format
dos2unix -c mac TADcallsByTool_Dixon.bed

# extract TopDom calls with 500M coverage and 50kb resolution and move the bed file to the data diirectory
awk -F "\t" '{if(($4=="TopDom") && ($5=="500M") && ($6==50000)) {print}}' TADcallsByTool_Dixon.bed > ../TAD.Dixon.500M.50kb.bed

## Recombination Hotspots
# This requires a fai file for the hg19 reference genome

# Seperating maternal and paternal crossovers and  discard all random chromsomes
# zcat aau1043_DataS4.gz | awk -F "\t" ('{if(($5=='M')&&($1!~/chr.._.*/)){printf "%s\t%s\t%s\n",$1,$2,$3}}') > crossovers_locations_maternal.bed
# zcat aau1043_DataS4.gz | awk -F "\t" ('{if(($5=='P')&&($1!~/chr.._.*/)){printf "%s\t%s\t%s\n",$1,$2,$3}}') > crossovers_locations_paternal.bed
#
# # liftOver crossover locations to hg19
# ./liftOver crossovers_locations_maternal.bed hg38ToHg19.over.chain.gz crossovers_locations_maternal_hg19.bed
# ./liftOver crossovers_locations_paternal.bed hg38ToHg19.over.chain.gz crossovers_locations_paternal_hg19.bed
#
# # create windowed hg19 reference genome
# # this requires a fai file for the hg19 genome
# bedtools makewindows -g hg19.fasta.fai -w 10000 > hg19.10kb.windows.bed
#
# # compute coverage
# bedtools coverage -a hg19.10kb.windows.bed -b crossover_locations_maternal_hg19.bed > crossover_maternal_coverage_hg19.bed
# bedtools coverage -a hg19.10kb.windows.bed -b crossover_locations_paternal_hg19.bed > crossover_paternal_coverage_hg19.bed
