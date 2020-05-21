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

            # correct for chromosome notation
            if self.data[0][:3] == 'chr':
                self.chr = self.data[0]
            else:
                self.chr = 'chr' + self.data[0]

            self.start = int(self.data[1])
            self.data = {column_name: self.data[3 + idx]
                         for idx, column_name in enumerate(column_names)}
            self.info = {element.split('=')[0]: element.split(
                '=')[1] for element in self.data['INFO'].split(';')}
            self.end = self.start + np.abs(int(self.info['SVLEN']))
        else:
            super().__init__(line, column_names)
        self.tads = []
        self.annotations = {}
        self.annotation_distances = {}

    def __str__(self):
        """CNV Objects are represented by overlapping TADs which in turn are represented by their annotations."""
        tads = "\n".join([str(tad) for tad in self.tads])
        return f'{self.chr}\t{self.start}\t{self.end}\nTADS\n{tads}'

    def calculate_overlap_and_distances(self, feature_type):
        """Calculates the distance and overlap for each annotation in the same TAD as the CNV."""
        self.annotation_distances['TAD_boundaries'] = []
        self.annotations['TAD_stability'] = []
        if self.tads:
            for tad in self.tads:
                boundary_overlaps = [self.start <
                                     tad.start, self.end > tad.end]
                boundary_distances = [self.start -
                                      tad.start, tad.end - self.end]

                self.annotation_distances['TAD_boundaries'].extend(
                    [0 if overlap else boundary_distances[i] for i, overlap in enumerate(boundary_overlaps)])

                if feature_type == 'extended':
                    self.annotations['TAD_stability'].append(
                        float(tad.data['stability_1']))
                    self.annotations['TAD_stability'].append(
                        float(tad.data['stability_2']))

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
                        self.annotation_distances[annotation_name].append(
                            distance)

            # get the TAD stability for the closest TAD boundaries (accounts for CNVs overlapping multiple boundaries)
            if feature_type == 'extended':
                minimal_distances = [i for i, distance in enumerate(self.annotation_distances['TAD_boundaries']) if distance == min(
                    self.annotation_distances['TAD_boundaries'])]
                self.annotations['TAD_stability'] = [
                    self.annotations['TAD_stability'][idx] for idx in minimal_distances]

    def annotate(self, feature_type):
        """Return a vector containing the feature vector for this annotated cnv"""
        self.calculate_overlap_and_distances(feature_type)

        if feature_type == 'distance':
            features = self.get_distance_features()
        if feature_type == 'extended':
            features = self.get_extended_features()

        self.features = features
        return features

    def get_distance_features(self):
        """Returns distance features which are either directly derived from the TADs.
        The output is a numpy feature vector."""
        features = [min(self.annotation_distances[annotation])
                    for annotation in self.annotations]
        return np.array(features, dtype=float)

    def get_extended_features(self):
        """Returns extended features which are either directly derived from the TADs or based on the CNV itself.
        The output is a numpy feature vector."""
        LOEUF = np.nan
        HI = np.nan
        phastcon = np.nan
        Log_odd_HI = np.nan
        exon_overlap = 0
        enhancer_overlap = 0
        overlapping_genes = []
        overlapping_enhancers = []
        minimal_normalized_interaction_overlap = np.nan

        try:
            # get the maximal proportional overlap with interacting fragments of the genes in the same TAD normalized by the gene's LOEUF
            minimal_normalized_interaction_overlap = max([gene.get_interaction_overlap(self) / float(gene.data['LOEUF']) for gene in self.annotations['GENES'] if float(gene.data['LOEUF']) != 0], default=np.nan)
        except (ValueError, KeyError, IndexError) as e:
            pass

        try:
            #  get overlapping genes
            overlapping_genes = np.array(self.annotations['GENES'])[np.where(
                np.array(self.annotation_distances['GENES']) == 0)]
        except (ValueError, KeyError, IndexError) as e:
            pass

        if len(overlapping_genes) > 0:
            try:
                # get the minimal LOEUF score of the overlapping genes
                LOEUF = np.nanmin([float(gene.data['LOEUF'])
                                   for gene in overlapping_genes])
            except (ValueError, KeyError, IndexError) as e:
                pass

            try:
                HIs = [float(gene.data['HI']) for gene in overlapping_genes]
                HSs = [1 - hi for hi in HIs]
                HI_division = np.divide(1 - np.prod(HSs), np.prod(HSs))
                # compute the HI log Odds score
                if HI_division != 0:
                    Log_odd_HI = np.log(HI_division)
                # get the maximum probability of being haploinsufficient for the overlapping gnes
                HI = np.nanmax(HIs)
            except (ValueError, KeyError, IndexError) as e:
                pass

            try:
                # get the highest proportional exon overlap of all overlapping genes
                exon_overlap = np.max([gene.get_exon_overlap(self)
                                       for gene in overlapping_genes])
            except (ValueError, KeyError, IndexError) as e:
                pass
        else:
            try:
                # get the LOEUF score of the closest gene
                LOEUF = float(self.annotations['GENES'][np.argmin(
                    self.annotation_distances['GENES'])].data['LOEUF'])
            except (ValueError, KeyError, IndexError) as e:
                pass
            try:
                # get the HI score of the closest gene
                HI = float(self.annotations['GENES'][np.argmin(
                    self.annotation_distances['GENES'])].data['HI'])
            except (ValueError, KeyError, IndexError) as e:
                pass

        try:
            # get overlapping enhancers
            overlapping_enhancers = np.array(self.annotations['ENHANCERS'])[np.where(
                np.array(self.annotation_distances['ENHANCERS']) == 0)]
        except (ValueError, KeyError, IndexError) as e:
            pass

        if len(overlapping_enhancers) > 0:
            try:
                # get the maximum enhancer overlap of all overlapping enhancers
                enhancer_overlap = np.max([utils.getOverlap((self.start, self.end), (enhancer.start, enhancer.end)) / (
                    enhancer.end - enhancer.start) for enhancer in overlapping_enhancers])
            except (ValueError, KeyError, IndexError) as e:
                pass

            try:
                # get the maximal Phastcon value of all overlapping enhancer
                phastcon = np.nanmax([float(enhancer.data['Phastcon'])
                                      for enhancer in overlapping_enhancers])
            except (ValueError, KeyError, IndexError) as e:
                pass
        else:
            try:
                phastcon = float(self.annotations['ENHANCERS'][np.argmin(
                    self.annotation_distances['ENHANCERS'])].data['Phastcon'])
            except (ValueError, KeyError, IndexError) as e:
                pass

        features = [len(overlapping_genes), len(overlapping_enhancers), min(self.annotation_distances['TAD_boundaries']), max(self.annotations['TAD_stability']), min(self.annotation_distances['GENES'], default=np.nan), min(
            self.annotation_distances['ENHANCERS'], default=np.nan), min(self.annotation_distances['DDG2P'], default=np.nan), LOEUF, phastcon, HI, min(self.annotation_distances['CTCF'], default=np.nan), Log_odd_HI, exon_overlap, minimal_normalized_interaction_overlap]
        return np.array(features, dtype=float)
