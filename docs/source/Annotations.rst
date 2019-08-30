:: _Annotations:

Annotations
===========

In this page I eleobarte on the annotations and intolerance metrics. Each of the annotations can be downloaded using the ``data_sources.sh`` script. 

Enhancers
---------

One of the most commonly used sets of predicted enhancers are the FANTOM5 enhancers cite:`d-andersson2014atlas`. The FANTOM5 expression atlas aimed at creating a comprehensive overview of functional annotations with focus on promoters cite:`d-forrest2014promoter`. Using so called *cap analysis gene expression* (CAGE) the authors identified transcription start sites (TSS) and created a cell-type specific functionally annotated mammalian transcriptome. In this process they discovered that experimentally determined active enhancers and predicted enhancers with increased reporter activity show increased CAGE activity. Based on the patterns of CAGE activity, bidirectoral capped RNA seemed to be a signature feature of active enhancers. Given this feature, the authors predicted in total 43,011 cell-type specific and ubiquitous activate enhancer candidates.

Genes
-----
The majority of genes that I used as annotations are genes annotated with loss-of-function metrics i.e. LOEUF and pLI scores, originally provided by GnomAD cite:`d-collins2019open. The DDG2P genes are reportedly associated with developmental disorders and have been curated by DECIPHER cite:`d-firth2009decipher`. 

CTCF Binding Sites
------------------
I obtained the CTCF binding sites from the USCS genome browser cite:`d-kent2002human` (Cell Line: AG09319 Replicate 2).


pLI
---
The *probability of being LoF intolerant* or pLI is originally based on the ExAC study cite:`d-lek2016analysis`. The study analyzed high-quality exomes of 60,706 individuals. The variants called from the exomes of these individuals were used to investigate the frequency with which specific genes are hit by mutations of different classes. The frequency of variants hitting genes was compared to the expectation under a model of random mutation. The deviation from the random mutation model was quantified using a Z-score, which was observed to be increased for missense and protein-truncating mutations compared to synonymous variants, indicating a higher constraint of genes for LoF mutations. To account for sequence length the authors employed an expectation maximization (EM) algorithm, separating the genes into three categories: null (observed number of LoF mutations is approximately equal to the expected number), recessive (the observed number of mutations is less or equal to 50\% of the expected number) and haploinsufficient (the observed number is less or equal to 10\% of the expected number). The pLI is equal to the probability of the gene belonging to the haploinsufficient class. Scores less or equal than 0.1 indicate a LoF tolerant gene and scores greater or equal to 0.9 indicate LoF intolerant gene. The pLI scores were later updated using the data from the GnomAD database cite:`d-karczewski2019variation`.  


LOEUF
-----
The *LoF observed/expected upper bound fraction* or LOEUF is a recently developed measurement of LoF intolerance. The measurement is a product of the GnomAD database cite:`d-karczewski2019variation` and as suggested by the authors an improved version of the pLI score. The pLI score is ideally used as a dichotomous measurement where genes are either intolerant (pLI greater or equal to 0.9) or tolerant  (pLI less or equal to 0.1) to LoF mutations. The increase of potential LoF variants due to the 125,748 exomes and 15,708 genomes sequenced in the GnomAD project, in combination with a refined model, enabled the authors to construct a continuous rather than dichotomous measurement. They computed a gene-specific expected/observed ratio of LoF mutations and a confidence interval around the ratio. For the final annotation they used the upper bound of this confidence interval since it represents the most conservative measurement. So even though pLI and LOEUF scores can be used in a very similar fashion the actual values are almost opposites. While lower LOEUF scores indicate a lower tolerance to LoF mutation, low pLi scores indicate a high tolerance to LoF mutations.


HI
--
Haploinsufficiency is the inability of a gene to maintain normal function, when one of the copies is missing. It can be considered as a major contributing factor to genetic disease. When investigating the effect of CNVs spanning one or more genes, haploinsufficiency can be a viable predictive factor of the functional impact of a CNV cite:`d-huang2010characterising`. Assuming a gene is of high importance to the normal function of an organism, the impact of a CNV such as a deletion can dramatically change depending on the ability of the gene to function with only a single copy. Huang et al. cite:`d-huang2010characterising` curated a list of 1,079 haplosufficient (HS) genes based on CNVs detected in 8,548 supposedly healthy individuals. If a gene was unambiguously and repeatedly hit by CNVs, it was considered to be haploinsufficent, since the individuals did not show any apparent effect of the missing or duplicated gene copy. The HS genes were compared to a set of known haploinsufficent (HI) genes described in literature cite:`d-dang2008identification,seidman2002transcription`. The comparison was based on a selection of functional, genomic and evolutionary features which included: The length of the gene, the spliced transcript, the 3' UTR and and coding sequence, the number of exons, the median Genomic Evolutionary Rate Profiling (GERP) score cite:`d-cooper2005distribution` of both promotor and gene sequence, tissue-specific expression patterns based on the GNF atlas cite:`d-su2004gene` and  protein-protein as well as gene interaction networks. In a first analysis the difference in distributions of the features were investigated using Mann-Whitney-U and Fisher-exact tests. The authors observed significant differences between HS and HI genes. The results supported the possibility of constructing a classifier to distinguish between the two groups. A linear discriminant analysis (LDA) was used as a supervised machine learning algorithm and showed reasonably accurate results (AUC=0.81). This model was applied to all genes in order to calculate a gene specific probability of being haploinsufficient (p(HI)). Similar to the pLI and LOEUF values, I used this measurement in the prioritization of CNVs to further partition the set of genes into more specific and potentially more informative groups. Furthermore, the authors proposed to compute a so called log-odds (LOD) score to assess the pathogenicity of deletions. The score compares the probability that none of the genes affected by the deletion is haploinsufficent with the probability that at least one of genes is haploinsufficent.


PhastCons
---------
The PhastCons score is a metric that is designed to identify the conservation of individual loci in the genome cite:`d-siepel2005evolutionarily`. The computation is based on a so called phylogenetic hidden Markov model or Phylo-HMM cite:`d-felsenstein1996hidden` with a conserved state and a nonconserved state. Each of the states is associated with a phylogentic model i.e. a probability distribution over alignment columns. The branches of the model associated with the conserved state are scaled by a parameter $p$, which refers to the rate of substitutions in conserved regions in proportion to the rate of substitutions in nonconserved regions. The computation of substitution probabilities of the phylogenetic model depends on the length of the branches, since they represent the expected number of changes per site. Different branch length therefore result in different substitution rates. The free parameters i.e. :math: `\mu` and :math:`v` were estimated using an EM algorithm. In the expectation step posterior expected counts are computed based on the number of times each distinct alignment column is produced by each state of the HMM and the number of times a transition between the two states occurs. In the maximization step the free parameters are optimized with regards to a measurement combining the previously computed posterior expected counts. The resulting optimized parameters can then be used to compute posterior probabilities that each site was generated by the conserved state using the the forward/backward algorithm cite:`d-durbin1998biological`. 
I used PhastCons conservation scores based on a Phylo-HMM trained on the multiple alignment of ninety-nine vertebrate genomes. The model is also referred to as the PhastCons-100way model. I averaged the single-nucleotide conservation cores over the sequence of enhancers to quantify their conservation. 

.. bibliography:: references.bib
   :style: plain
   :labelprefix: D
   :keyprefix: d-

