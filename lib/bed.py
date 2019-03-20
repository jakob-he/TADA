"""
A basic class to handle elements in BED files. This class serves as a parent class for TADs, Enhancer and genes.
This class assumes that the line passed to the constructor contains at least the following elements in this exact format:
chr \t  start \t  end
"""

class Bed:
    def __init__(self, line, column_names):
        """Basic class containing the information from a BED file element.
        Args:
            line: The line of bed file (This is just a string).
            column_names: Additional columns except chr,start and end in the bed file sorted by appearence. (e.g. gene name)
        Returns:
            A new Bed object.
        """
        self.data = line.strip().split('\t')
        self.chr = self.data[0]
        self.start = int(self.data[1])
        self.end = int(self.data[2])
        self.data = {column_name: self.data[3 + idx]
                     for idx, column_name in enumerate(column_names)}

    def __eq__(self, other):
        """Objects of this class are equal if the start and end position are equal."""
        return self.end == other.end and self.start == other.end

    def __lt__(self, other):
        """Bed objects are sorted based on the end position."""
        return self.start < other.start

    def __str__(self):
        """Bed objects are repsented by their chromsome end and start position"""
        return f'{self.chr}\t{self.start}\t{self.end}'
