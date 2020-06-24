# CHROM HMM
#DEL Pathogenic
gat-run.py --verbose=5 --segment-file=non_pathogenic.DEL.no_overlap.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_transcribed.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Poised_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Strong_Enhancer.bed   --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Enhancer.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Segmental_duplications.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/HARs_hg19.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Insulator.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_transition.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_elongation.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Repetitive.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Heterochrom.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Polycomb.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Active_Enhancer.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=del.chromhmm.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=del_chromhmm_non_pathogenic.log > del_chromhmm_non_pathogenic.tsv

# DEL Non-Pathogenic
gat-run.py --verbose=5 --segment-file=pathogenic.no_overlap.DEL.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_transcribed.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Poised_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Strong_Enhancer.bed   --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Enhancer.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Segmental_duplications.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/HARs_hg19.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Insulator.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_transition.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_elongation.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Repetitive.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Heterochrom.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Polycomb.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Active_Enhancer.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=del.chromhmm.pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=del_chromhmm_pathogenic.log > del_chromhmm_pathogenic.tsv


gat-compare.py del.chromhmm.non-pathogenic.segment-overlap.counts.tsv.gz del.chromhmm.pathogenic.segment-overlap.counts.tsv.gz > del.chromhmm.compared_calls.txt

python3 plot_supp_figures.py -a del_chromhmm_pathogenic.tsv -b del_chromhmm_non_pathogenic.tsv -c del.chromhmm.compared_calls.txt -o Deletions_ChromHMM_01 -q 0.01


# DUP Npn-Pathogenic
gat-run.py --verbose=5 --segment-file=non_pathogenic.DUP.no_overlap.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_transcribed.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Poised_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Strong_Enhancer.bed   --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Enhancer.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Segmental_duplications.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/HARs_hg19.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Insulator.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_transition.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_elongation.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Repetitive.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Heterochrom.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Polycomb.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Active_Enhancer.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=dup.chromhmm.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=dup_chromhmm_non_pathogenic.log > dup_chromhmm_non_pathogenic.tsv

# DUP Pathogenic
gat-run.py --verbose=5 --segment-file=pathogenic.no_overlap.DUP.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_transcribed.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Poised_Promotor.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Strong_Enhancer.bed   --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Weak_Enhancer.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/Segmental_duplications.bed --annotation-file=../GAT_ANNOTATIONS/OTHER/HARs_hg19.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Insulator.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_transition.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Transciptional_elongation.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Repetitive.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Heterochrom.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Polycomb.bed --annotation-file=../GAT_ANNOTATIONS/ChromHMM/BroadHMM_hESC_Active_Enhancer.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=dup.chromhmm.pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=dup_chromhmm_pathogenic.log > dup_chromhmm_pathogenic.tsv

gat-compare.py dup.chromhmm.non-pathogenic.segment-overlap.counts.tsv.gz dup.chromhmm.pathogenic.segment-overlap.counts.tsv.gz > dup.chromhmm.compared_calls.txt

python3 plot_supp_figures.py -a dup_chromhmm_pathogenic.tsv -b dup_chromhmm_non_pathogenic.tsv -c dup.chromhmm.compared_calls.txt -o Duplications_ChromHMM_01 -q 0.01

# TADs

# # DELs Non-Pathogenic
gat-run.py --verbose=5 --segment-file=merged_gnomad_eichler_biobank_DGV.DEL.no_overlap.size_matched.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_without_functional_elements.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_1_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.9_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.5_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.9_genes.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.5_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.1_enhancer.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=del.tads.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=del_tads_non_pathogenic.log > del_tads_non_pathogenic.tsv

# DELs Pathogenic
gat-run.py --verbose=5 --segment-file=pathogenic.no_overlap.DEL.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_without_functional_elements.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_1_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.9_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.5_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.9_genes.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.5_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.1_enhancer.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=del.tads.pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=del_tads_pathogenic.log > del_tads_pathogenic.tsv

gat-compare.py del.tads.non-pathogenic.segment-overlap.counts.tsv.gz del.tads.pathogenic.segment-overlap.counts.tsv.gz > del.tads.compared_calls.txt

python3 plot_supp_figures.py -a del_tads_pathogenic.tsv -b del_tads_non_pathogenic.tsv -c del.tads.compared_calls.txt -o Deletions_TADs_01 -q 0.01

# DUPs Non-Pathogenic
gat-run.py --verbose=5 --segment-file=merged_gnomad_DGV.DUP.no_overlap.size_matched.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_without_functional_elements.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_1_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.9_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.5_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.9_genes.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.5_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.1_enhancer.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=dup.tads.non-pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=dup_tads_non_pathogenic.log > dup_tads_non_pathogenic.tsv

# DUPs Pathogenic
gat-run.py --verbose=5 --segment-file=pathogenic.no_overlap.DUP.size_matched.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_without_functional_elements.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_1_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.9_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.5_enhancer.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.9_genes.bed  --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.5_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_pli_0.1_genes.bed --annotation-file=../GAT_ANNOTATIONS/TADs/tads_with_phastcon_0.1_enhancer.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --workspace-file=../gat/tutorial/TutorialIntervalOverlap/contigs_ungapped_autosomes.bed --isochore-file=../GAT_ANNOTATIONS/Constantini_isochores_autosomes.bed --random-seed=7 --counter=segment-overlap --output-counts-pattern=dup.tads.pathogenic.%s.counts.tsv.gz --nbuckets=10000 --bucket-size=960 --num-threads=4 --num-samples=10000  --log=dup_tads_pathogenic.log > dup_tads_pathogenic.tsv

gat-compare.py dup.tads.non-pathogenic.segment-overlap.counts.tsv.gz dup.tads.pathogenic.segment-overlap.counts.tsv.gz > dup.tads.compared_calls.txt

python3 plot_supp_figures.py -a dup_tads_pathogenic.tsv -b dup_tads_non_pathogenic.tsv -c dup.tads.compared_calls.txt -o Duplications_TADs_01 -q 0.01