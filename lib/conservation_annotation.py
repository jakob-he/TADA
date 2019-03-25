"""This script can be used to annotate Enhancers or any other bed element with PhastCon conservation scores."""
import pathlib
import argparse
import .utils as utils

import pyBigWig

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bedfile',help='Input bedfile which is going to be annotated.')
    parser.add_argument('-c','--conservation_file',default='../data/ANNOTATION/hg19.100way.phastCons.bw',help='BigWig file containing the conservation scores.')
    parser.add_argument('-o','--output',help='Output name and location for the annotated bed file.')
    args = parser.parse_args()
    return args





def main():
    #read input arguments
    args = argparser()

    #read bed elements from bedfile
    beds = utils.objects_from_file(args.bedfile)

    #open BigWig file
    bw = pyBigWig.open(args.conservation_file)

    #iterate through bedfile and find annotation
    conservation_scores = [bw.stats(bed.chr,bed.start,bed.end) for bed in beds]

    #save annotated bed elements
    with open(args.output,'w') as output:
        for idx,bed in enumerate(beds):
            additional_columns = "\t".join(self.data.values())
            output.write(f'{bed.chr}\t{bed.start}\t{bed.end}\t{additional_columns}\t{conservation_scores[idx]}')


















if __name__ == '__main__':
    main()
