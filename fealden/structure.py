import os
import sys
import configparser
import math
import itertools

dir_path = os.path.dirname(os.path.abspath(__file__))
config = configparser.ConfigParser()
config.read(os.path.join(dir_path, "config.ini"))
RNA_PATH = config["RNAstructure"]["Path"]

try:
    sys.path.append(os.path.join(RNA_PATH, "exe"))
    os.environ["DATAPATH"] = os.path.join(RNA_PATH, "data_tables")
    import RNAstructure
except Exception as error:
    print("RNAstructure could not be found; check config.ini")
    print(error)


class RNAfolder:
    """
    RNAfolder composes an instance of the RNA class in RNAstructure
    and generates useful attribute data
    reflecting the secondary structure(s) of the input sequence.

    NOTE: currently hard coded to interpret input as DNA

    Parameters:

        seq --> The sequence of nucleotides that you wish to
            be queried for structure prediction

    """

    def __init__(self, seq):
        self.seq = seq.upper()
        self.RNAobj = RNAstructure.RNA.fromString(f"{self.seq}", backbone="dna")
        self.RNAobj.FoldSingleStrand(percent=15, window=0)
        self.number_folds = self.RNAobj.GetStructureNumber()
        self.point_list = [
            self.get_coordinate_list(structure_num=i + 1)
            for i in range(self.number_folds)
        ]
        self.structure_dict = [dict() for _ in range(self.number_folds)]
        self.make_fold_dict()

    def make_fold_dict(self):
        """
        make_fold_dict populates the list of dictionaries
        (self.structure_dict) defined in __init__
        """
        seq_len = len(self.RNAobj)
        for i, each_dict in enumerate(self.structure_dict):
            each_dict["deltaG"] = self.RNAobj.GetFreeEnergy(i + 1)
            pairs_list = []
            for base in range(seq_len):
                pairs_list.append(
                    [base + 1, self.RNAobj.GetPair(base + 1, structurenumber=i + 1)]
                )
            each_dict["bps"] = pairs_list
            # dict of deltaG's with corresponding folding list for each fold

    def get_coordinate_list(self, h=10, w=10, structure_num=1) -> list:
        self.RNAobj.DetermineDrawingCoordinates(h, w, structure_num)
        comp = [
            [
                self.RNAobj.GetNucleotideXCoordinate(i),
                self.RNAobj.GetNucleotideYCoordinate(i),
            ]
            for i in range(1, len(self.RNAobj) + 1)
        ]
        return comp

    def dist_from_index(
        self, index1: int, index2: int, structure_num=1
    ) -> float:  # INDEX STARTS @ 1
        """
        Returns the distance between two nucleotides, index1 and index2,
        where index2 > index1 and the index values of draw_list start at 1 not 0
        (i.e. base pair 12 corresponds to index 12


        PARAMETERS:
            draw_list --> a list of nucleotide coordinates [x, y] returned from
                find_coordinates for a given structure where each [x, y] pair
                corresponds to a base pair of the structure

            index1 --> the nucleotide index of the first base pair used for the
                distance calculation (index starts at 1)
            index2 --> the nucleotide index of the second base pair used for the
                distance calculation (index starts at 1).

        """
        points = (
            self.point_list[structure_num - 1][index1 - 1]
            + self.point_list[structure_num - 1][index2 - 1]
        )
        return RNAfolder.dist(*points)

    def find_best_tag(self, segment=3) -> list:
        tags_for_structure = []
        tail = [
            i
            for i, j in zip(self.RNAobj.iterIndices(), self.RNAobj.iterNucs())
            if j == "T"
        ]
        # tail = range(len(self)-segment, len(self))
        tag_list = [dict() for _ in range(self.number_folds)]
        prod = [[o, p] for o, p in itertools.product([1], tail[-segment:])]
        for i, a_dict in enumerate(tag_list):
            for [o, p] in prod:
                tags_for_structure.append((o, p, self.dist_from_index(o, p, i + 1)))

            # sort list of tuples (bp1, bp2, magnitude) based off largest magnitude
            # and set highest as the dictionary value for tag_list with the
            # structure number as the key
            a_dict[f"{i+1}"] = sorted(
                tags_for_structure, reverse=True, key=lambda points: points[2]
            )[0]

        return tag_list

    @staticmethod
    def dist(x1, y1, x2, y2, extra=False):
        dX = x2 - x1
        dY = y2 - y1
        mag = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if extra:
            return dX, dY, mag
        else:
            return mag

    def __len__(self):
        return self.RNAobj.GetSequenceLength()

    def __str__(self):
        return self.seq

    def __repr__(self):
        return f"RNAfolder instance\n sequence input: {self.seq}\n \
            number of structures: {self.number_folds}"
