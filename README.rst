TADA
====

Introduction
------------

The TAD annotation (TADA) tools is designed to annotate CNVs based on functional annotation with respect to their regulatory context i.e. TADs. TADA allows to determine the functional impact of CNVs either by annotation and manual filtering or automated classification. The default Random Forest models are trained on curated sets of pathogenic and non-pathogenic variants (see our preprint for details (https://www.biorxiv.org/content/10.1101/2020.06.30.180711v1.abstract)). New CNVs can therefore be annotated with the probability of being pathogenic i.e. a *pathogenicity score*. A simplified workflow of TADA is shown in the figure below.

.. figure:: TADA_Workflow.jpg
   :scale: 30 %
   :alt: TADA Workflow


Installation
------------

- Install TADA via GitHub::

     git clone https://github.com/jakob-he/TADA
     python setup.py install

- To ensure a successfull installtion a test protocol can be executed::

      python setup.py test

Usage
-----

There are three use cases for TADA that can be executed using the command line:

	* Annotation of CNVs for manual filtering
	* Prediction of pathogenicity
	* Training a classification model

The first use case i.e. annotation of CNVs for manual filtering requires a BED- or VCF-file with CNVs. The path to the CNV needs to specified in a *config* file. TADA provides two default config files for deletions and duplications, respectively. The default output of the annotation is a pickled CNV file containing all the annotated CNV objects (*Annotated_CNVs.p*) and a CSV-file where each line refers to an annotated CNV (*Annotated_CNVs.csv*). The annotation is executed as follows::

    general usage: annotate_cnvs [-h] [-p] [-c CONFIG] [-o OUTPUT]

    default: annotate_cnvs -c config_file -o output_directory

The process of pathogenicity prediction is similar to the annotation. In a first step, the CNVs passed to the script are annotated. Then a previously defined classification model defined in the config file is used to compute a *pathogenicity score*. This output refers to the probability that the CNV is pathogenic. The score is added to the resulting CSV-file. The default deletion and duplication models are Random Forests trained on curated sets of pathogenic and non-pathogenic variants. However, the user can also specify a different model. It is recommend to train the model using TADA as described below, since it needs to satisfy certain requirements. The recommended way of using the prediction script is::

    general usage: predict_variants [-h] [-c CONFIG] [-o OUTPUT]

    default: predict_variants -c config_file -o output_directory

The training of a classification model is an advanced use of TADA. In addition to two sets of CNVs e.g. pathogenic and non-pathogenic, a set of annotations needs to be specified. It is also possible to use the default annotation set. Based on these annotation a feature set is computed. For user defined annotations the features are distances of a CNV to the individual genomic elements. The model can the be trained on the CNVs set with the given feature set as follows::

   general usage: classification_run  [-h] [-c CONFIG] [-fs] [-o OUTPUT] [-gridcv] [-rs RANDOM_SEED]

   default usage: classification_run -c config_file -o output_directory

The *fs* option allows to produce multiple visualizations for feature selection. The visualizations include the permutation based feature importance and a partial correlation based node graph.
If required the classification run can be executed with the *gridcv* option to find an optimal parameter set. For reproducability the *rs* option can be set to a specific integer.
