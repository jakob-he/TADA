.. _pathogenicityprediction:

Pathogenicity Prediction
========================

In this page I describe the data collection, feature selection and model selection that resulted in the default classification model. 

Data Collection
---------------

There are two sets of CNVs that I used for training of the default classification model: Pathogenic variants and non-pathogenic variants. The non-pathogenic variants are a combined set from three different data sources: A publication by Audano et al. :cite:`c-audano2019characterizing`, the recently published collection of SVs provided by the GnomAD consortium :cite:`c-collins2019open` and a set of CNVs called from the UK Biobank data set :cite:`c-aguirre2019phenome`.

The Audano set originally contained 99604 insertions, deletions and inversions. The variants were called from high coverage long-read sequencing experiments of 15 individuals. Since almost all of the annotations and the other CNVs are called based on the GRCh37 reference genome, I used the UCSC tool *LiftOver* :cite:`c-kent2002human` to map the Audano variants from GRCh38 to GRCh37. In this process 3748 variants were lost. In the original publication, the variants were partitioned by their population frequency into four groups: *Shared* - Present in all individuals, *Major* - Present in half or more of all samples, *Polymorphic* - present in more than one but less than 50 \% of all samples and *Singleton* - present in only one sample. Even though the individuals are supposedly healthy, I only included the Shared and Major variant groups in further analysis, to avoid rare and potentially pathogenic variants. Additionally, I discarded all variant types other than CNVs. This resulted in a set of 5605 deletions.

The GnomAD SVs were called from deep whole genome sequencing using short reads of 14,981 individuals across multiple global populations. The total set of 498,257 SVs is comprised of 13 different mutational classes. In the first step of filtering, I discarded all variants which are not sequence resolved or do not pass any of the other filters described in the original publication. I reduced this set of 347,390 variants by filtering for variants with an allele frequency greater or equal to 0.1, to avoid the inclusion of rare and potentially pathogenic variants. The last filter I applied removed all variants smaller than 50 bp and those which were labeled as mutational classes other than duplications and deletions i.e. CNVs. The final filtered set of GnomAD variants includes a total 6933 CNVs.

The set of variants called from the UK biobank data is not yet publicly available. However, the authors of the paper generously provided me with a curated population-wide call set comprised of 275180 CNVs (171825 deletions, 103355 insertions). The calls are based on SNP-chip experiments and included the number of individuals in the cohort supporting each CNV. As a recent publication :cite:`c-weedon2019very` shows, the application of SNP-chips for the identification of rare variants is extremely unreliable. This study compared SNP calls from both SNP-chips and next-generation sequencing. The authors observed a very low positive predictive value for heterozygous genotypes when investigating variants with a frequency less or equal to 0.001 in the UK biobank cohort. Even though these result can not be directly applied to CNV calls based on SNP-chips data, they do provide an indication that low frequency calls should be treated with caution in further analysis. I therefore discarded all variants of the UK biobank CNV call set, which are supported by less than 3 individuals. The final set of UK biobank variant includes 30,532 variants.

The combined set of 63870 non-pathogenic variants was further filtered by discarding all variants larger than 500kb and matching the size distribution of pathogenic variants to the non-pathogenic variants. The size distributions were matched using an empirical cumulative distribution function (ECDF). I sorted an array containing all unique sizes of the non-pathogenic set and created forty evenly spaced bins in the range of the unique sizes. For each of these bins, I counted the number of non-pathogenic CNVs that are of less or equal size. If the number of pathogenic variants was less than the number of non-pathogenic variants in the individual bin, I included all pathogenic variants in the size matched set. If more pathogenic variants satisfied the size requirements of the bin, I randomly sampled the same number of pathogenic variants as non-pathogenic variants in the bin. The total number of size matched and filtered non-pathogenic variants is 6826.

The pathogenic variants are a combined set of 20990 DECIPHER CNVs :cite:`c-firth2009decipher` and 723 variants provided by the Mundlos AG :cite:`c-Mundlos`. The variants from DECIPHER are categorized by their pathogenic effect. I included all variants that are categorized as pathogenic, likely pathogenic or unknown. Any variants with a length smaller than 50 bp were discarded during filtering. In a second filtering step the variants were binned based on the size distribution of non-pathogenic variants and then an equal number of pathogenic and non-pathogenic variants were drawn from the respective data sources. The total number of pathogenic variants is therefore equal to the number of non-pathogenic variants i.e. 6826. 

Feature Selection
-----------------

The feature set I used to annotate the pathogenic and non-pathogenic CNVs consist of the following annotations:

	* GnomAD gene distance
	* FANTOM5 enhancer distance
	* CTCF binding site distance
	* TAD boundary distance
	* DDG2P gene distance
	* Gene LOEUF score
	* Enhancer PhastCon score
	* Gene Haploinsufficiency
	* Haploinsufficiency Log-Odds score

The distance always refers to the closest functional element located in the same TAD/TADs as the corresponding CNV. Similarly, the intolerance scores were taken from the closest Gene and Enhancer. I provide a more detailed description of the annotations and intolerance metrics in the :ref:`background` section.

Model Selection
---------------

I tested the performance of multiple machine learning approaches i.e. logistic regression, Random Forest and Balanced Random Forest. Each of the classifiers is still available when using TADA. However, the balanced Random Forest provided the best mean average precision of a 10-fold cross-validation. The default model is therefore a Balanced Random Forest trained on the previously described feature set. 


.. bibliography:: references.bib
   :style: plain
   :labelprefix: C
   :keyprefix: c-
