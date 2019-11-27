"""Helper functions to parse BED files."""
# standart libraries
import pathlib
import json
import gzip

# own libraries
from .bed import Bed
from .bed_class import BedClass

# third party libararies
import numpy as np
import pandas as pd
from scipy import stats, linalg


def objects_from_file(path, cls_string, column_names=[], **kwargs):
    """Load a BED file and return a list of Bed objects""
    Args:
        path: Path to the BED file.
        cls_string: A string that matches one of the bed classes (e.g. Gene).
    """
    path = validate_file(path)
    bed_class = BedClass.from_str(cls_string).get_class()


    if path.suffix == '.gz':
        open_path = gzip.open(path.absolute())
    else:
        open_path = path.open()

    bed_objects = []
    vcf = False
    for line in open_path:
        if not type(line) == str:
            line = line.decode('utf-8')
        #handle headers of bed CNV files
        if line.startswith('##'):
            vcf = True
            continue
        if line.startswith('#'):
            column_names = line.strip().split('\t')[3:]
            continue
        bed_objects.append(bed_class(line, column_names, **dict(kwargs,vcf=vcf))
    open_path.close()
    return bed_objects


def validate_file(path):
    """Check if the path is a valid BED file and return a pathlib.Path object."""
    path = pathlib.Path(path)

    # check if path is a valid file
    if not path.is_file():
        raise Exception(f'{path} is not a valid path')

    # check if path is a bed or txt file
    # TODO this is just prelimenary check
    if path.suffix not in ['.bed','.txt','.gz','.vcf']:
        raise Exception(f'{path} is not a bed,txt or gz file')

    return path


def read_result_file(path):
    """Load the data from a result file as a pandas dataframe.
    Currently only a vector containing the 10-fold CV results is returned."""
    path = pathlib.Path(path)

    # check if path is a valid file
    if not path.is_file():
        raise Exception(f'{path} is not a valid path')

    with path.open() as results:
        for line in results:
            if line.startswith('10'):
                cv_avg = float(line.strip().split(':')[1])
            elif line.startswith('non-pathogenic'):
                support = int(line.strip().split('     ')[-1])
    return cv_avg, support


def create_chr_dictionary_from_beds(beds: [Bed]):
    """Create a dictionary based on bed objects with the chromsomes as keys."""
    bed_dict = {}
    for bed in beds:
        if not bed.chr in bed_dict:
            bed_dict[bed.chr] = [bed]
            continue
        bed_dict[bed.chr].append(bed)
    bed_dict = {key: sorted(bed_dict[key]) for key in bed_dict}
    return bed_dict


def is_in(bed_list, reference_bed):
    """Returns True if the first element of the list of bed objects is in the reference_bed"""
    # return False if the list of bed objects contains no element
    if not bed_list:
        return False

    # check if the first element is in the reference bed object
    if bed_list[0].start < reference_bed.end:
        return True
    else:
        return False

def to_bed(bed_elements,output,label=''):
    """Saves the input as bed file.
    The input can either be a dict with chromsomes as keys and list of bed elements as items or a list of bed elements.
    """
    if type(bed_elements) == list:
        bedlist_to_bed(bed_elements,output,label)
    elif type(bed_elements) == dict:
        chrom_dict_to_bed(bed_elements,output,label)
    else:
        print('The input has to be a dictionary or list with Bed elements.')

def bedlist_to_bed(bedlist,output,label):
    with open(output,'w') as output:
        for bed in bedlist:
            output.write(f'{bed.chr}\t{bed.start}\t{bed.end}')
            if label:
                output.write(f'\t{label}')
            output.write('\n')

def chrom_dict_to_bed(chrom_dict,output,label):
    with open(output,'w') as output:
        for chrom in chrom_dict:
            for bed in chrom_dict[chrom]:
                output.write(f'{bed.chr}\t{bed.start}\t{bed.end}')
                if label:
                    output.write(f'\t{label}')
                output.write('\n')



def reduce_dict(dictionary, keys):
    """Returns a dictionary containing only the input keys"""
    return {key: (dictionary[key] if key in dictionary else []) for key in keys}


def create_annotated_bed_dict(bed_dict, annotation_dicts, annotate=False,filer_exons=False,feature_type='binary'):
    """Annotates every BED flement with the overlapping annotation.
    For each BED element in a chromosome the function iterates through the sorted annotations as long as the
    start position of any of the first annotations is less than the end position of the BED element.
    If one of the elements satisfies this condition there are three options:
        1. The annotation ends before the BED element -> discard it, since there is no overlapping BED element in the data set.
        2. The annotation starts in or before the BED element and ends in the BED element -> Add to the BED element's annotations and remove annotation from the list.
        2. The annotation starts in or before the BED element and does not end in the BED element -> Add to BED element's annotations but keep it for other BED elements.

    Args:
        bed_dict: A dictionary with chromsomes as keys and the corresponding BED elements as values.
        annotation_dict: A list of dictionaries with chromosomes as keys and the corresponding annotation elements as values.
    """
    # reduce genes and enhancers to chromsomes were tads are available
    annotation_dicts = {key:reduce_dict(
        dictionary, bed_dict.keys()) for key, dictionary in annotation_dicts.items()}

    # iterate through chromsomes
    for chrom in bed_dict:
        for bed_element in bed_dict[chrom]:
            annotation_queue = {}

            for annotation_name, annotation_dict in annotation_dicts.items():
                bed_element.annotations[annotation_name] = []
                annotation_queue[annotation_name] = []

            while any(is_in(annotation_dict[chrom],bed_element) for annotation_dict in annotation_dicts.values()):
                for annotation_name, annotation_dict in annotation_dicts.items():
                    if is_in(annotation_dict[chrom],bed_element):
                        if annotation_dict[chrom][0].end < bed_element.start:
                            annotation_dict[chrom].pop(0)
                        elif annotation_dict[chrom][0].end <= bed_element.end:
                            bed_element.annotations[annotation_name].append(annotation_dict[chrom].pop(0))
                        elif annotation_dict[chrom][0].end > bed_element.end:
                            bed_element.annotations[annotation_name].append(annotation_dict[chrom][0])
                            annotation_queue[annotation_name].append(annotation_dict[chrom].pop(0))

            annotation_dict[chrom] = annotation_queue[annotation_name] + annotation_dict[chrom]
            if annotate:
                bed_element.annotate(feature_type)
            if filer_exons:
                bed_element.filter_exons()

    return bed_dict


def annotate_cnvs(tad_dict, cnv_dict):
    """Finds all TADs overlapping with the CNV, then constructs a new chrom dictionary with the annotated CNVs.
    The function iterates through the TADs one chromsome at a time. For each TAD it checks every CNVs were
    the start position is either less or equal to the TADs start position. If that is the case there are four possiblities:
        1. The CNV ends before the TAD -> this CNV is either not in any of the available TADs or ended in between TADs.
        2. The CNV ends in the TAD but starts before it -> append the TAD to the CNV and move the CNV to the list of annotated CNVs and change boundary_spanning to True.
        3. The CNV ends and starts in the TAD -> append the TAD to the CNV and move the CNV to the list of annotated CNVs.
        4. Tne CNVs starts before the TAD but ends after it -> append the TAD to the CNV, keep it in the CNV dict and change the boundary_spanning attribute to True.
    """
    # reduce the cnvs to chromsomes were tads are available
    cnv_dict = reduce_dict(cnv_dict,tad_dict.keys())

    # create empty list
    annotated_cnvs = []

    #iterate through CNVs one chromsome at a time
    for chrom in cnv_dict:
        for tad in tad_dict[chrom]:
            cnv_queue = []
            while is_in(cnv_dict[chrom],tad):
                if cnv_dict[chrom][0].end <= tad.end and cnv_dict[chrom][0].end >= tad.start:
                    if cnv_dict[chrom][0].start < tad.start:
                        cnv_dict[chrom][0].boundary_spanning = True
                    cnv_dict[chrom][0].tads.append(tad)
                    annotated_cnvs.append(cnv_dict[chrom].pop(0))
                elif cnv_dict[chrom][0].end > tad.end:
                    cnv_dict[chrom][0].boundary_spanning = True
                    cnv_dict[chrom][0].tads.append(tad)
                    cnv_queue.append(cnv_dict[chrom].pop(0))
                else:
                    annotated_cnvs.append(cnv_dict[chrom].pop(0))
            cnv_dict[chrom] = cnv_queue + cnv_dict[chrom]
        annotated_cnvs.extend(cnv_dict[chrom])

    return create_chr_dictionary_from_beds(annotated_cnvs)


def getOverlap(interval_a,interval_b):
    """Returns the overlap of two intervals"""
    return max(0,min(interval_a[1],interval_b[1])-max(interval_a[0],interval_b[0]))

def getDistance(interval_a,interval_b):
    """Returns the distance between two intervals"""
    return max(interval_a[0]-interval_b[1],interval_b[0]-interval_a[1])

def phi_coeff(array_1,array_2):
    """Implementation of the phi coefficient computation."""
    cont_tab = pd.crosstab(array_1,array_2)
    print(cont_tab)
    return (cont_tab[1][1]*cont_tab[0][0] - cont_tab[0][1]*cont_tab[1][0])/np.sqrt(sum(cont_tab[0])*sum(cont_tab[1])*sum(cont_tab[:][1])*sum(cont_tab[:][0]))

def partial_corr(feature_df,feature_type):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in the feature dataframe, controlling
    for the remaining variables in the dataframe. The implementation is based on https://gist.github.com/fabianp/9396204419c7b638d38f.
    Args:
    feature_df = pandas DataFrame object, shape (n, p)
        Array with the different variables. Each column of C is taken as a variable
    Output:
    P_corr = pandas DataFrame object, shape (p, p)
        P_corr[i, j] contains the partial correlation of feature_df[:, i] and feature_df[:, j] controlling
        for the remaining variables in the feature dataframe.
    """
    C = feature_df.values
    p = C.shape[1]
    P_corr = np.zeros((p, p), dtype=np.float)

    for i in range(p):
        P_corr[i, i] = 1
        for j in range(i+1, p):
            idx = np.ones(p, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j = linalg.lstsq(C[:, idx], C[:, i])[0]

            res_j = C[:, j] - C[:, idx].dot(beta_i)
            res_i = C[:, i] - C[:, idx].dot(beta_j)
            corr = stats.pearsonr(res_i, res_j)[0]
            P_corr[i, j] = corr
            P_corr[j, i] = corr
    P_corr = pd.DataFrame(P_corr,columns=feature_df.columns)
    P_corr.rename({idx:column for idx,column in enumerate(feature_df.columns)},inplace=True)
    return P_corr
