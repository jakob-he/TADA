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
            A new CNV object.
        """
        if vcf:
            self.data = line.strip().split('\t')

            #correct for chromosome notation
            if self.data[0][:3]=='chr':
                self.chr = self.data[0]
            else:
                self.chr = 'chr' + self.data[0]

            self.start = int(self.data[1])
            self.data = {column_name: self.data[3 + idx]
                         for idx, column_name in enumerate(column_names)}
            self.info = {split(element,'=')[0]:split(element,'=')[0] for element in split(self.data['INFO'],'')}
            self.end = self.start + np.abs(self.info['SVLEN'])
        else:
            super().__init__(line, column_names)
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
        self.annotation_distances['TAD_boundaries'] = []
        self.annotations['TAD_contact_pvalue'] = []
        if self.tads:
            for tad in self.tads:
                tad_distances = [abs(self.start - tad.start),
                                 abs(tad.end - self.end)]

                if self.boundary_spanning:
                    self.annotation_distances['TAD_boundaries'].append(0)
                else:
                    self.annotation_distances['TAD_boundaries'].append(
                        min(tad_distances))

                if tad.data['pvalue_1']:
                    self.annotations['TAD_contact_pvalue'].append(
                        tad.data[f'pvalue_{np.argmin(tad_distances)+1}'])

                for annotation_name, annotations in tad.annotations.items():
                    self.annotations[annotation_name] = annotations
                    self.annotation_distances[annotation_name] = []
                    for annotation in annotations:
                        overlap = utils.getOverlap([self.start, self.end], [
                                                   annotation.start, annotation.end])
                        if overlap > 0:
                            overlap = 1
                            distance = 0
                        else:
                            distance = utils.getDistance([self.start, self.end], [
                                                         annotation.start, annotation.end])

                        if binary:
                            self.annotation_distances[annotation_name].append(
                                overlap)
                        else:
                            self.annotation_distances[annotation_name].append(
                                distance)

    def annotate(self, feature_type):
        """Return a vector containing the feature vector for this annotated cnv"""
        if 'binary' in feature_type:
            self.calculate_overlap_and_distances('binary')
        else:
            self.calculate_overlap_and_distances('continuous')

        if feature_type == 'basic_binary':
            features = self.get_basic_binary_features()
        if feature_type == 'extended_binary':
            features = self.get_extended_binary_features()
        if feature_type == 'basic_continuous':
            features = self.get_basic_continuous_features()
        if feature_type == 'extended_continuous':
            features = self.get_extended_continuous_features()

        self.features = features
        return features

    def get_basic_binary_features(self):
        """Returns basic binary features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [self.boundary_spanning, any(overlap for overlap in self.annotation_distances['genes']), any(
            overlap for overlap in self.annotation_distances['enhancers'])]
        return np.array(features)

    def get_extended_binary_features(self):
        """Returns extended binary features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [self.boundary_spanning, any(overlap for overlap in self.annotation_distances['genes']), any(overlap for overlap in self.annotation_distances['enhancers']), any(
            overlap for overlap in self.annotation_distances['DDG2P']), any(overlap for overlap in self.annotation_distances['CTCF'])]
        return np.array(features)

    def get_basic_continuous_features(self):
        """Returns basic continuous features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy boolean feature vector."""
        features = [min(self.annotation_distances['TAD_boundaries']), min(
            self.annotation_distances['genes'], default=np.nan), min(self.annotation_distances['enhancers'], default=np.nan)]
        return np.array(features, dtype=float)

    def get_extended_continuous_features(self):
        """Returns extended continuous features which are either directly derived from the TADs or         print(self.data)based on the CNV itself.
        The output is a numpy boolean feature vector."""
        LOEUF = np.nan
        HI = np.nan
        phastcon = np.nan
        Log_odd_HI = np.nan
        exon_overlap = 0
        enhancer_overlap = 0
        try:
            overlapping_genes = np.array(self.annotations['genes'])[np.where(
                np.array(self.annotation_distances['genes']) == 0)]
            HIs = [float(gene.data['HI']) for gene in overlapping_genes]
            HSs = [1 - hi for hi in HIs]
            HI_division = np.divide(1 - np.prod(HSs), np.prod(HSs))
            if HI_division != 0:
                Log_odd_HI = np.log(HI_division)
            LOEUF = float(self.annotations['genes'][np.argmin(
                self.annotation_distances['genes'])].data['LOEUF'])
            HI = float(self.annotations['genes'][np.argmin(
                self.annotation_distances['genes'])].data['HI'])
            phastcon = float(self.annotations['enhancers'][np.argmin(
                self.annotation_distances['enhancers'])].data['Phastcon'])
            exon_overlap = self.annotations['genes'][np.argmin(
                self.annotation_distances['genes'])].get_exon_overlap(self)
            closest_enhancer = self.annotations['enhancers'][np.argmin(
                self.annotation_distances['enhancers'])]
            enhancer_overlap = utils.getOverlap((self.start, self.end), (
                closest_enhancer.start, closest_enhancer.end)) / (closest_enhancer.end - closest_enhancer.start)
        except (ValueError, KeyError, IndexError) as e:
            err_message = e

        features = [len(np.array(self.annotations['genes'])[np.where(np.array(self.annotation_distances['genes']) == 0)]), len(np.array(self.annotations['enhancers'])[np.where(np.array(self.annotation_distances['enhancers']) == 0)]), min(self.annotation_distances['TAD_boundaries']), self.annotations['TAD_contact_pvalue'][np.argmin(self.annotation_distances['TAD_boundaries'])], min(self.annotation_distances['genes'], default=np.nan), min(self.annotation_distances['enhancers'], default=np.nan), min(self.annotation_distances['DDG2P'], default=np.nan), LOEUF, phastcon, HI, min(self.annotation_distances['CTCF'], default=np.nan), Log_odd_HI, exon_overlap]
        return np.array(features, dtype=float)
