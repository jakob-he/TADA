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
        exons = "\n".join([str(exon) for exon in self.annotations['exons']])
        return f'{self.chr}\t{self.start}\t{self.end}\t{additional_columns}\nExons:\n{exons}'

    def filter_exons(self):
        """Discard all exons corresponding to a different gene name"""
        matching_exons = []
        for exon in self.annotations['exons']:
            if self.data['gene_name'] in exon.data['gene_name'].split(','):
                matching_exons.append(exon)
        self.annotations['exons'] = matching_exons

    def get_exon_overlap(self,cnv):
        """Return proportion of exons overlapping with the given cnv"""
        if self.annotations['exons']:
            exon_prop = sum([1 if utils.getOverlap((cnv.start,cnv.end),(exon.start,exon.end)) > 0 else 0 for exon in self.annotations['exons']])/len(self.annotations['exons'])
            return exon_prop
        else:
            return 0
