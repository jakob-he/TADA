"""This script can be used to annotate Enhancers or any other bed element with PhastCon conservation scores."""
import pathlib
import argparse

from TADA.lib import utils

import pyBigWig

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bedfile',help='Input bedfile which is going to be annotated.')
    parser.add_argument('-c','--conservation_file', help='BigWig file containing the conservation scores.')
    parser.add_argument('-o','--output',help='Output name and location for the annotated bed file.')
    args = parser.parse_args()
    return args



def main():
    #read input arguments
    args = argparser()

    #read bed elements from bedfile
    beds = utils.objects_from_file(args.bedfile,'enhancer')

    #open BigWig file
    bw = pyBigWig.open(args.conservation_file)

    #iterate through bedfile and find annotation
    conservation_scores = [bw.stats(bed.chr,bed.start,bed.end) for bed in beds]

    #save annotated bed elements
    with open(args.output,'w') as output:
        for idx,bed in enumerate(beds):
            additional_columns = "\t".join(bed.data.values())
            output.write(f'{bed.chr}\t{bed.start}\t{bed.end}{additional_columns}\t{",".join([str(score) for score in conservation_scores[idx]])}\n')


















if __name__ == '__main__':
    main()
