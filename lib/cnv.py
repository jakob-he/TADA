"""
The CNV class is an extension of the Bed class.
"""
from .bed import Bed
from . import utils
import numpy as np

class CNV(Bed):
    def __init__(self, line, column_names, vcf=False):
        """A Class containing the information about a copy number variant.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. sequence)
        Returns:
            A new Gene object.
        """
        if vcf:
            self.data = line.strip().split('\t')
            location = self.data[1].split(':')
            self.chr = location[0]
            position = location[1].split('-')
            self.start = int(position[0])
            self.end = int(position[1])
        else:
            super().__init__(line,column_names)
        self.tads = []
        self.boundary_spanning = False
        self.annotations = {}
        self.annotation_distances = {}

    def __str__(self):
        """CNV Objects are representated by all the columns entered in the bed file"""
        tads = "\n".join([str(tad) for tad in self.tads])
        return f'{self.chr}\t{self.start}\t{self.end}\nTADS\n{tads}'


    def calculate_overlap_and_distances(self, feature_type):
        """Calculates the distance and overlap for each annotation in the same TAD as the CNV. Overlap is always binary."""
        binary = feature_type == 'binary'
        self.annotation_distances['TAD_boundaries']=[]
        if self.tads:
            for tad in self.tads:

                if self.boundary_spanning:
                    self.annotation_distances['TAD_boundaries'].append(0)
                else:
                    self.annotation_distances['TAD_boundaries'].append(min(self.start-tad.start,tad.end-self.end))

                for annotation_name, annotations in tad.annotations.items():
                    self.annotations[annotation_name] = annotations
                    self.annotation_distances[annotation_name] = []
                    for annotation in annotations:
                        overlap = utils.getOverlap([self.start,self.end],[annotation.start,annotation.end])
                        if overlap > 0:
                            overlap = 1
                            distance = 0
                        else:
                            distance = utils.getDistance([self.start,self.end],[annotation.start,annotation.end])

                        if binary:
                            self.annotation_distances[annotation_name].append(overlap)
                        else:
                            self.annotation_distances[annotation_name].append(distance)

    def annotate(self, feature_type):
        """Return a vector containing the feature vector for this annotated cnv"""
        if 'binary' in feature_type:
            self.calculate_overlap_and_distances('binary')
        else:
            self.calculate_overlap_and_distances('continuous')

        if feature_type=='basic_binary':
            features = self.get_basic_binary_features()
        if feature_type=='extended_binary':
            features = self.get_extended_binary_features()
        if feature_type=='basic_continuous':
            features = self.get_basic_continuous_features()
        if feature_type=='extended_continuous':
            features = self.get_extended_continuous_features()

        self.features = features
        return features

    def get_basic_binary_features(self):
        """Returns basic binary features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [self.boundary_spanning,any(overlap for overlap in self.annotation_distances['genes']),any(overlap for overlap in self.annotation_distances['enhancers'])]
        return np.array(features)

    def get_extended_binary_features(self):
        """Returns extended binary features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [self.boundary_spanning,any(overlap for overlap in self.annotation_distances['genes']),any(overlap for overlap in self.annotation_distances['enhancers']),any(overlap for overlap in self.annotation_distances['DDG2P']),any(overlap for overlap in self.annotation_distances['CTCF'])]
        return np.array(features)

    def get_basic_continuous_features(self):
        """Returns basic continuous features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [min(self.annotation_distances['TAD_boundaries']),min(self.annotation_distances['genes'],default=np.nan),min(self.annotation_distances['enhancers'],default=np.nan)]
        return np.array(features,dtype=float)

    def get_extended_continuous_features(self):
        """Returns extended continuous features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        LOEUF = np.nan
        HI = np.nan
        phastcon = np.nan
        Log_odd_HI = np.nan
        try:
            overlapping_genes = np.array(self.annotations['genes'])[np.where(np.array(self.annotation_distances['genes']) == 0)]
            HIs  = [float(gene.data['HI']) for gene in overlapping_genes]
            HSs = [1-hi for hi in HIs]
            HI_division = np.divide(1-np.prod(HSs),np.prod(HSs))
            if HI_division != 0:
                Log_odd_HI = np.log(HI_division)
            LOEUF = float(self.annotations['genes'][np.argmin(self.annotation_distances['genes'])].data['LOEUF'])
            HI = float(self.annotations['genes'][np.argmin(self.annotation_distances['genes'])].data['HI'])
            phastcon = float(self.annotations['enhancers'][np.argmin(self.annotation_distances['enhancers'])].data['Phastcon'])
        except (ValueError,KeyError,IndexError) as e:
            err_message = e

        features = [min(self.annotation_distances['TAD_boundaries']),min(self.annotation_distances['genes'],default=np.nan),min(self.annotation_distances['enhancers'],default=np.nan),min(self.annotation_distances['DDG2P'],default=np.nan),LOEUF,phastcon,HI,min(self.annotation_distances['CTCF'],default=np.nan),Log_odd_HI]
        return np.array(features,dtype=float)
