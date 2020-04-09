"""
The Gene class is an extension of the Bed class.
"""
from .bed import Bed
from . import utils
import numpy as np

class Gene(Bed):
    def __init__(self, line, column_names, vcf=False):
        """A Class containing the information about a Gene.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. sequence)
        Returns:
            A new Gene object.
        """
        super().__init__(line,column_names)
        self.annotations = {}

    def __str__(self):
        """CNV Objects are representated by all the columns entered in the bed file"""
        additional_columns = "\t".join(self.data.values())
        exons = "\n\t".join([str(exon) for exon in self.annotations['EXONS']])
        interactions = "\n\t".join([str(fragment) for fragment in self.annotations['fragments']])
        return f'{self.chr}\t{self.start}\t{self.end}\t{additional_columns}\n\tExons:\n\t{exons}\n\tInteracting Fragments:\n\t{interactions}' 

    def filter_exons(self):
        """Discard all exons corresponding to a different gene name"""
        matching_exons = []
        for exon in self.annotations['EXONS']:
            if self.data['gene_name'] in exon.data['gene_name'].split(','):
                matching_exons.append(exon)
        self.annotations['EXONS'] = matching_exons

    def filter_interactions(self):
        """Check for overlapping promotors and append interacting fragments"""
        self.annotations['fragments'] = []
        # Set interacting fragments to zero if more than one promotor overlaps
        if len(self.annotations['POINT']) == 1:
            for fragment in self.annotations['POINT'][0].data['fragments'].split(','):
                chr, location = fragment.split(':')
                start, end = location.split('-')
                if chr == self.chr:
                    self.annotations['fragments'].append((chr,int(start),int(end)))


    def get_exon_overlap(self,cnv):
        """Return proportion of exons overlapping with the given cnv"""
        if self.annotations['EXONS']:
            exon_prop = sum([1 if utils.getOverlap((cnv.start,cnv.end),(exon.start,exon.end)) > 0 else 0 for exon in self.annotations['EXONS']])/len(self.annotations['EXONS'])
            return exon_prop
        else:
            return 0

    def get_interaction_overlap(self,cnv):
        """Return proportion of interacting fragments overlapping with the given cnv normalized by the LOEUF score"""
        if self.annotations['fragments']:
            fragment_prop = sum([1 if utils.getOverlap((cnv.start,cnv.end),(fragment[1],fragment[2])) > 0 else 0 for fragment in self.annotations['fragments']])/len(self.annotations['fragments'])
            return float(fragment_prop)
        else:
            return 0
