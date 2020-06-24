##Deletions


#Final: Extended non-pathogenic
#Binary isochore Enrichment analysis of the non-pathogenic variants including Biobank data for FANTOM5 enhancer and Genes.
gat-run.py --verbose=5 --segment-file=non_pathogenic.DEL.no_overlap.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_75.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_90.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_pli_loeuf_HI.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_loeuf_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_loeuf_0_1.bed  --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_HI_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_HI_0_1.bed --annotation-file=../GAT_ANNOTATIONS/GENES/DDG2p_genes.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Dixon_2015_stability_formatted_GAT.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/IMR90_CTCF_peak_idr_optimal.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/extended_telomeres_5MB_ungapped_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=del.extended.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-samples=10000 --num-threads=4 --log=del_extended_non_pathogenic.log > del_extended_non_pathogenic.tsv

#Final: Extended pathogenic
#Binary isochore Enrichment analysis of the non-pathogenic variants including Biobank data for FANTOM5 enhancer and Genes.
gat-run.py --verbose=5 --segment-file=pathogenic.DEL.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_75.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_90.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_pli_loeuf_HI.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_loeuf_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_loeuf_0_1.bed  --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_HI_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_HI_0_1.bed --annotation-file=../GAT_ANNOTATIONS/GENES/DDG2p_genes.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Dixon_2015_stability_formatted_GAT.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/IMR90_CTCF_peak_idr_optimal.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/extended_telomeres_5MB_ungapped_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed  --random-seed=7  --counter=segment-overlap --nbuckets=10000 --bucket-size=960  --output-counts-pattern=del.extended.pathogenic.%s.counts.tsv.gz --num-samples=10000 --num-threads=4  --log=del_extended_pathogenic.log > del_extended_pathogenic.tsv

gat-compare.py del.extended.non-pathogenic.segment-overlap.counts.tsv.gz del.extended.pathogenic.segment-overlap.counts.tsv.gz > del.extended.compared_calls.txt

python3 plot_main_figures.py -a del_extended_pathogenic.tsv -b del_extended_non_pathogenic.tsv -c del.extended.compared_calls.txt -o Deletions_extended_01 -q 0.01


##Duplications

#Final: Extended non-pathogenic
#Binary isochore Enrichment analysis of the non-pathogenic variants including Biobank data for FANTOM5 enhancer and Genes.
gat-run.py --verbose=5 --segment-file=non_pathogenic.DUP.no_overlap.size_matched.bed --annotation-file=/=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_75.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_90.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_pli_loeuf_HI.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_loeuf_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_loeuf_0_1.bed  --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_HI_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_HI_0_1.bed --annotation-file=../GAT_ANNOTATIONS/GENES/DDG2p_genes.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Dixon_2015_stability_formatted_GAT.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/IMR90_CTCF_peak_idr_optimal.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/extended_telomeres_5MB_ungapped_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed  --random-seed=7 --counter=segment-overlap --output-counts-pattern=dup.extended.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-samples=10000 --num-threads=4  --log=dup_extended_non_pathogenic.log > dup_extended_non_pathogenic.tsv

#Final: Extended pathogenic
#Binary isochore Enrichment analysis of the non-pathogenic variants including Biobank data for FANTOM5 enhancer and Genes.
gat-run.py --verbose=5 --segment-file=pathogenic.DUP.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_75.bed --annotation-file=../GAT_ANNOTATIONS/ENHANCER/fantom5_enhancer_phastcon_average_filtered_90.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_pli_loeuf_HI.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_loeuf_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_loeuf_0_1.bed  --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_high_HI_0_9.bed --annotation-file=../GAT_ANNOTATIONS/GENES/gnomad_genes_low_HI_0_1.bed --annotation-file=../GAT_ANNOTATIONS/GENES/DDG2p_genes.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Dixon_2015_stability_formatted_GAT.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/IMR90_CTCF_peak_idr_optimal.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/extended_telomeres_5MB_ungapped_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed  --random-seed=7  --counter=segment-overlap --nbuckets=10000 --bucket-size=960  --output-counts-pattern=dup.extended.pathogenic.%s.counts.tsv.gz --num-samples=10000 --num-threads=4  --log=dup_extended_pathogenic.log > dup_extended_pathogenic.tsv

gat-compare.py dup.extended.non-pathogenic.segment-overlap.counts.tsv.gz dup.extended.pathogenic.segment-overlap.counts.tsv.gz > dup.extended.compared_calls.txt

python3 plot_main_figures.py -a dup_extended_pathogenic.tsv -b dup_extended_non_pathogenic.tsv -c dup.extended.compared_calls.txt -o Duplications_extended_01 -q 0.01