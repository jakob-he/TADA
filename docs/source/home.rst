Getting Started
===============

Introduction
------------

The TAD annotation (TADA) tools is designed to annotate CNVs based on functional annotation with respect to their regulatory context i.e. TADs. TADA allows to determine the functional impact of CNVs either by annotation and manual filtering or automated classification. The default model provided in the ``data`` sub-directory is trained on pathogenic and non-pathogenic variants. New CNVs can therefore be annotated with the probability of being pathogenicity. The origin and methodology of the default model is further described in :ref:`pathogenicityprediction`. For a detailed description of CNVs, their pathogenicity and TADs, see :ref:`background`.

Installation
------------

- Install TADA via GitHub::

     git clone https://github.com/jakob-he/TADA
     python setup.py install

Usage
-----

There are three use cases for TADA that can be executed using the command line:

	* Annotation of CNVs for manual filtering
	* Prediction of pathogenicity
	* Training a classification model

The first use case i.e. annotation of CNVs for manual filtering only requires a BED- or VCF-file with CNVs. By default the CNVs are annotated using the feature set also used for the pathogenicity prediction. The individual annotations are listed in :ref:`pathogenicityprediction`. The default output of the annotation are a pickled CNV file containing all the annotated CNV objects (*Annotated_CNVs.p*) and a CSV-file where each line refers to an annotated CNV (*Annotated_CNVs.csv*). It is important to mention, that the location of the default TADs needs to changed if the script is not executed in the TADA directory. The output location is by default the current directory but can also be defined by the user. The annotation is executed as follows::

    general usage: annotate_cnvs [-h] [-t TADS] -c CNVS [-vcf] [-p] [-csv] [-f FEATURES] [-o OUTPUT]

    default: annotate_cnvs -c path_to_CNV_file -o example_annotated_CNVs.p
	
The prediction of pathogenicity is similar to the annotation. In a first step, the CNVs passed to the script are annotated. Then a previously defined classification model is used to compute the output of a sigmoid function. This output refers, for instance, to the probability that the CNV is pathogenic. This probability is added to the resulting CSV-file. The default model is a Balanced Random Forest trained on pathogenic and non-pathogenic variants. However, the user can also specify a different model. It is recommend to train the model using TADA, since it needs to satisfy certain requirements and can currently only be based on one of the four available feature sets. The recommended way of using the prediction script is::

    general usage: predict_variants [-h] -c CNVs [-t TADS] [-m MODEL] [-o OUTPUT] [-l] [-f FEATURE] [-csv]

    default: predict_variants -c path_to_CNV-file -o example_annotated_CNVs.p

The training of a classification model is a more advanced use of TADA. In addition to two sets of annotated CNVs e.g. pathogenic and non-pathogenic, a classifier and feature set have to be specified. The user can choose between three different classifiers: **lr**-Logistic Regression, **rf**-Random Forest and **brf**-Balanced Random Forest. The predefined feature sets are: **basic_binary**, **extended_binary**, **basic_continuous** and **extended_continuous**. Basic feature sets only use enhancers and genes while extended feature sets include CTCF binding sites and DDG2P genes. The most comprehensive feature set i.e. the extended continuous feature set also included intolerance metrics. If necessary additional parameters for the classifier e.g. number of trees can be passed to the **kwargs** option. The script can also be used to generate multiple plots to explore the feature set such as partial correlation node graphs (the **fs** option). Unfortunately the pathogenic variant are not yet publicly available. As an alternative a model can be trained on frequent and rare CNVs (results in a very inaccurate model). The pickled classifier is saved as *classifier-model.p* to the user-defined output location .e.g. *brf-model.p*. The general and default usage is::

   general usage: classification_run [-h] -c1 CNVS1 -c2 CNVS2 -f FEATURE [-fs] -cls CLS [-kw KWARGS] [-o OUTPUT]

   default usage: classification_run -c1 frequenct_CNVs.p -c2 rare_CNVs.p -f extended_continuous -cls brf -o path_to_outputdir


Tests
-----

The package includes several test to ensure a sucsessfull installation::

    git clone https://github.com/jakob-he/TADA
    python setup.py test


