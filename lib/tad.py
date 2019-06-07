"""
The Tad class is an extension of the bed class.
"""
from .bed import Bed


class Tad(Bed):
    def __init__(self, line, column_names, **kwargs):
        """Basic class containing the information from a BED file element.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. enhancer name)
        Returns:
            A new Tad object.
        """
        super().__init__(line, column_names, **kwargs)
        self.annotations = {}

    def __str__(self):
        """An object of the Tad class is representeted by the chromsome, start and end position
        as well as the recursive representation of its genes and enhancer"""
        representation = f'{self.chr}\t{self.start}\t{self.end}\n'
        for name, annotations in self.annotations.items():
            annotation_repr = ','.join([str(annotation) for annotation in annotations])
            representation += f'{name}:{annotation_repr}'
        return representation

    def contains_high_pLI_gene(self):
        """Returns True if the TAD contains a gene with pLi greater or equal to 0.9.
        This requires genes with pLI values in the range (0,1)."""
        if self.annotations['genes']:
            pLIs = [gene.data['pLI'] for gene in self.annotations['genes']]
            return any(float(pLI) == 1 for pLI in pLIs if pLI != 'NA')

    def contains_highly_conserved_enhancer(self):
        """Returns True if the TAD contains an enhancers with Phastcon value greater or equal to 0.6.
        This is the largest threshold showing significant enrichment of pathogenic variants (GAT)."""
        if self.annotations['enhancers']:
            phastcons = [enhancer.data['Phastcon'] for enhancer in self.annotations['enhancers']]
            return any(float(phastcon) == 1 for phastcon in phastcons if phastcon != 'None')

    def annotate(self, feature_type):
        """Annotates the TAD with a set of features."""
        if feature_type == 'binary':
            self.high_pLI = self.contains_high_pLI_gene()
            self.high_Phastcon = self.contains_highly_conserved_enhancer()
