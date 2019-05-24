"""Helper functions to parse BED files."""
# standart libraries
import pathlib
import json
import gzip

# own libraries
from .bed import Bed
from .bed_class import BedClass


def objects_from_file(path, cls_string, column_names=[]):
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
    for line in open_path:
        if not type(line) == str:
            line = line.decode('utf-8')
        #handle headers of bed CNV files
        if line.startswith('#'):
            column_names = line.strip().split('\t')[3:]
            continue
        bed_objects.append(bed_class(line, column_names))

    open_path.close()
    return bed_objects


def validate_file(path):
    "Check if the path is a valid BED file and return a pathlib.Path object."
    path = pathlib.Path(path)

    # check if path is a valid file
    if not path.is_file():
        raise Exception(f'{path} is not a valid path')

    # check if path is a bed or txt file
    # TODO this is just prelimenary check
    if path.suffix not in ['.bed','.txt','.gz']:
        raise Exception(f'{path} is not a bed,txt or gz file')

    return path


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


def create_annotated_tad_dict(tad_dict, gene_dict, enhancer_dict):
    """Annotates every TAD with the overlapping genes and enhancers.
    For each TAD in a chromosome the function iterates through the sorted gene and enhancer lists as long as the
    start position of either the first gene or first enhancer is less than the end position of the TAD.
    If one of the elements satisfies this condition there are three options:
        1. The element ends before the TAD -> discard it, since there is no overlapping TAD in the data set.
        2. The element starts in or before the TAD and ends in the TAD -> Add to TADs elements and remove element from the list.
        2. The element starts in or before the TAD and does not in the TAD -> Add to TAD elements but keep it for other TADs.

    Args:
        tad_dict: A dictionary with chromsomes as keys and the corresponding Tad elements as values.
        gene_dict: The same as tad_dict with Gene objects.
        enhancer_dict: The same as tad_dict with Enhancer objects.
    """
    # reduce genes and enhancers to chromsomes were tads are available
    gene_dict, enhancer_dict = [reduce_dict(
        dictionary, tad_dict.keys()) for dictionary in [gene_dict, enhancer_dict]]

    # iterate through chromsomes
    for chrom in tad_dict:
        for tad in tad_dict[chrom]:
            gene_queue = []
            enhancer_queue = []
            while is_in(gene_dict[chrom], tad) or is_in(enhancer_dict[chrom], tad):
                if is_in(gene_dict[chrom], tad):
                    if gene_dict[chrom][0].end < tad.start:
                        gene_dict[chrom].pop(0)
                    elif gene_dict[chrom][0].end <= tad.end:
                        tad.genes.append(gene_dict[chrom].pop(0))
                    elif gene_dict[chrom][0].end > tad.end:
                        tad.genes.append(gene_dict[chrom][0])
                        gene_queue.append(gene_dict[chrom].pop(0))

                if is_in(enhancer_dict[chrom], tad):
                    if enhancer_dict[chrom][0].end < tad.start:
                        enhancer_dict[chrom].pop(0)
                    elif enhancer_dict[chrom][0].end <= tad.end:
                        tad.enhancers.append(enhancer_dict[chrom].pop(0))
                    elif enhancer_dict[chrom][0].end > tad.end:
                        tad.enhancers.append(enhancer_dict[chrom][0])
                        enhancer_queue.append(enhancer_dict[chrom].pop(0))

            enhancer_dict[chrom] = enhancer_queue + enhancer_dict[chrom]
            gene_dict[chrom] = gene_queue + gene_dict[chrom]

    return tad_dict


def annotate_cnvs(tad_dict, cnv_dict):
    """Finds all TADs overlapping with the CNV, then constructs a new dictionary with the annotated CNVs.
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

    for cnv in annotated_cnvs:
        cnv.calculate_overlap_and_distances()

    return create_chr_dictionary_from_beds(annotated_cnvs)


def getOverlap(interval_a,interval_b):
    """Returns the overlap of two intervals"""
    return max(0,min(interval_a[1],interval_b[1])-max(interval_a[0],interval_b[0]))

def getDistance(interval_a,interval_b):
    """Returns the distance between two intervals"""
    return max(interval_a[0]-interval_b[1],interval_b[0]-interval_a[1])
