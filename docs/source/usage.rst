Usage Instructions
==================

This page contains basic and more advanced instructions for the use of TADA.

Basic Instructions
------------------

The basic usage of TADA is described in the Introduction. I included annotated TADs, CNVs and a pretrained model in addition to the command line scripts. These files are used by default for the annotation, prediction and classification using TADA. They include TAD boundaries as defined by Dixon et al. called with TopDom which are annotated with the extended continuous feature set. I also added annotated frequent and rare CNVs. The frequent CNVs are a combined set of GnomAD (AF > 0.01) and Audano variants (Major and Shared). The rare CNVs are GnomAD variants with a frequency less or equal to 0.01. The pretrained default model is the result on pathogenic and non-pathogenic data-sets as described in :ref:`pathogenicityprediction`. The user can annotate a set of CNVs using the annotated TADs and predict pathogenicity with the default model. 

Advanced Usage
--------------

TADA also allows the user to annotate TADs and CNVs with BED-files that are not included in any of the default annotations. The script ``annotate_TADs`` requires a set of TAD boundaries which can be defined by the user and a collection of annotations. The TAD boundaries need to be BED-files, including the chromosome, start and end position, delimited by tabs. Annotations also require these three features but can include additional information. For instance, I used genes annotated with pLI scores in the annotation of the default TADs. In the current version, TADA only support intolerance metrics e.g. pLI and PhastCons that have been used for the default TADs and CNVs. However, it allows the user to add an unlimited amount of annotations. For each additional annotation a name has to be given to the command line script. The general usage and an example on how to add annotations is shown below::

    general usage: annotate_tads [-h] -t TADS [-an [ANNOTATION_NAMES [ANNOTATION_NAMES ...]]] [-af [ANNOTATION_FILES [ANNOTATION_FILES ...]]] [-o OUTPUT]
    
    #DDG2P example:
    annotate_tads -t Dixon_TADs.bed -an DDG2P -af DDG2P.bed -o DDG2P_annotated_TADs.p

This annotated set of TADs can then be used to annotated CNVs as described in the Introduction. However, including the new annotations in the classification is currently not supported. Only annotations of the four predefined feature set can be used for predicting pathogenicity or comparable measurements. In the **Tutorials** section the download and preprocessing of annotation is explained in more detail.


