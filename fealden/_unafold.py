import itertools
import math
import os
import re
import subprocess
import tempfile
from collections.abc import Iterator


class RNAfolder:

    """
    RNAfolder composes an instance of the RNA class using unafold/mfold routines
    and generates useful attribute data
    reflecting the secondary structure(s) of the input sequence.

    NOTE: currently hard coded to interpret input as DNA

    Parameters:

        seq --> The sequence of nucleotides that you wish to
            be queried for structure prediction

    """

    def __init__(self, seq: str) -> None:
        """Initialize RNAfolder object."""
        self.seq = seq.upper()
        self.ct_output = self.collect_unafold_ct(self.seq)
        self.number_folds = self.ct_output.count("dG")
        self.point_list = [
            self.get_coordinate_list(structure_num=i + 1)
            for i in range(self.number_folds)
        ]
        self.structure_dict: list[dict[str, float | list[list[int]]]] = [
            {} for _ in range(self.number_folds)
        ]
        self.make_fold_dict()

    @staticmethod
    def collect_unafold_ct(sequence: str) -> str:
        """Python wrapper to call UNAfold hybrid-ss-min without reading or writing files

        Uses modified hybrid-ss-min that prints .ct file contents to stdout."""

        with tempfile.TemporaryDirectory() as tmpdirname:
            hybrid_ss_min_location = os.getenv("HYBRID_SS_MIN")

            command = [
                f"{hybrid_ss_min_location}",
                "--mfold=15",
                "--sodium=0.15",
                "--magnesium=0.0005",
                "--NA=DNA",
                "/dev/stdin",
            ]
            hybrid_process = subprocess.Popen(
                command,
                text=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=tmpdirname,
            )

            _ = hybrid_process.communicate(input=f"{sequence}")[0]

            with open(f"{tmpdirname}/stdin.ct") as f:
                hybrid_output = f.read()

        return hybrid_output

    @staticmethod
    def parse_ct_to_folds(sequence: str) -> tuple[list[str], list[list[str]]]:
        def group_by_heading(
            some_source: list[str], heading: str = "dG"
        ) -> Iterator[list[str]]:
            buffer: list[str] = []
            for line in some_source:
                if heading in line:
                    if buffer:
                        yield buffer
                    buffer = [line]
                else:
                    buffer.append(line)
            yield buffer

        headings = []
        lines_list = []
        for heading_and_lines in group_by_heading(sequence.split("\n")):
            heading = heading_and_lines[0]
            headings.append(heading)
            lines = heading_and_lines[1:]
            lines_list.append(lines)

        return (headings[0:], lines_list[0:])

    @staticmethod
    def return_basepair(position: int, ct_lines: list[str]) -> int:
        my_line = ct_lines[position - 1]

        return int(my_line.split()[4])

    @staticmethod
    def return_free_energy(header_line: str) -> float:
        return float(header_line.split()[3])

    def make_fold_dict(self) -> None:
        """
        make_fold_dict populates the list of dictionaries
        (self.structure_dict) defined in __init__
        """
        seq_len = len(self.seq)
        headings, list_lines = self.parse_ct_to_folds(self.ct_output)
        for i, each_dict in enumerate(self.structure_dict):
            each_dict["deltaG"] = self.return_free_energy(headings[i])
            pairs_list: list[list[int]] = []
            for base in range(seq_len):
                pairs_list.append(
                    [base + 1, self.return_basepair(base + 1, list_lines[i])]
                )
            each_dict["bps"] = pairs_list
            # dict of deltaG's with corresponding folding list for each fold

    def get_coordinate_list(self, structure_num: int = 1) -> list[list[int]]:
        headings, list_lines = self.parse_ct_to_folds(self.ct_output)

        sir_graph_input = "\n".join(
            [headings[structure_num - 1], "\n".join(list_lines[structure_num - 1])]
        )

        def call_sir_graph(tmpdirname: str, ct_output: str) -> None:
            sir_graph_location = os.getenv("SIR_GRAPH")
            command = [
                f"{sir_graph_location}",
                "-p",
                "-o",
                f"{tmpdirname}/new.ps",
                "/dev/stdin",
            ]

            sir_graph_process = subprocess.Popen(
                command,
                text=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=tmpdirname,
            )

            _ = sir_graph_process.communicate(input=f"{ct_output}")

        with tempfile.TemporaryDirectory() as tmpdirname:
            call_sir_graph(tmpdirname, sir_graph_input)
            with open(f"{tmpdirname}/new.ps", encoding="ISO-8859-1") as f:
                contents = f.read()

        raw_pattern = r"(\d*\.\d*) (\d*\.\d*) m\n\(\w\)"

        pattern = re.compile(raw_pattern, re.MULTILINE)

        m = pattern.findall(contents)

        return [[int(float(i[0]) / 10), int(float(i[1]) / 10)] for i in m]

    def dist_from_index(
        self, index1: int, index2: int, structure_num: int = 1
    ) -> float | tuple[int, int, float]:  # INDEX STARTS @ 1
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
        return RNAfolder.dist(*points)  # type: ignore[arg-type]

    def find_best_tag(
        self, segment: int = 3
    ) -> list[dict[str, tuple[int, int, float]]]:
        tags_for_structure = []
        tail = [i for i, j in enumerate(self.seq) if j == "T"]
        # tail = range(len(self)-segment, len(self))
        tag_list: list[dict[str, tuple[int, int, float]]] = [
            {} for _ in range(self.number_folds)
        ]
        prod = [[o, p] for o, p in itertools.product([1], tail[-segment:])]
        for i, a_dict in enumerate(tag_list):
            for [o, p] in prod:
                tags_for_structure.append((o, p, self.dist_from_index(o, p, i + 1)))

            # sort list of tuples (bp1, bp2, magnitude) based off largest magnitude
            # and set highest as the dictionary value for tag_list with the
            # structure number as the key
            a_dict[f"{i+1}"] = sorted(  # type: ignore
                tags_for_structure, reverse=True, key=lambda points: points[2]
            )[0]

        return tag_list

    @staticmethod
    def dist(
        x1: int, y1: int, x2: int, y2: int, extra: bool = False
    ) -> float | tuple[int, int, float]:
        dX = x2 - x1
        dY = y2 - y1
        mag = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        if extra:
            return dX, dY, mag
        return mag

    def __len__(self) -> int:
        """Return sequence length."""
        return len(self.seq)

    def __str__(self) -> str:
        """Return string representation."""
        return self.seq

    def __repr__(self) -> str:
        """Return representation."""
        return f"RNAfolder instance\n sequence input: {self.seq}\n \
            number of structures: {self.number_folds}"
