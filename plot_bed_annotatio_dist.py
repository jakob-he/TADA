"""Plots the distribtion of a given column in a bed file"""
import pathlib
import argparse
from lib import plotting, utils


def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bedfile',help='Input bedfile which is going to be annotated.')
    parser.add_argument('-c','--columns',help='Comma seperated column names for the bed file except chr, start and end')
    parser.add_argument('-a','--annotation',help='The name of the column that should be plotted')
    parser.add_argument('-o','--output',help='Output name and location for the annotated bed file.')
    parser.add_argument('-p','--plottype',help='Type of the plot. "h"-histogram   "b"-boxplot   "v"-violinplot')
    args = parser.parse_args()
    return args





def main():
    #read input arguments
    args = argparser()

    #get column names
    column_names = args.columns.split(",")

    #read bed elements from bedfile
    beds = utils.objects_from_file(args.bedfile,'bed',column_names)

    #plotting
    plotting.plot_annotation_dist(beds, args.annotation, args.output, args.plottype, save=True)





if __name__ == '__main__':
    main()
