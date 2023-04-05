from unittest.mock import MagicMock, patch

from fealden._unafold import RNAfolder

SAMPLE_CT = """13	dG = 0.892	stdin
1	C	0	2	0	1	0	0
2	A	1	3	0	2	0	0
3	T	2	4	0	3	0	0
4	G	3	5	0	4	0	5
5	C	4	6	12	5	4	6
6	T	5	7	11	6	5	0
7	A	6	8	0	7	0	0
8	G	7	9	0	8	0	0
9	C	8	10	0	9	0	0
10	T	9	11	0	10	0	0
11	A	10	12	6	11	0	12
12	G	11	13	5	12	11	13
13	T	12	0	0	13	12	0
13	dG = 1.407	stdin
1	C	0	2	0	1	0	0
2	A	1	3	0	2	0	0
3	T	2	4	0	3	0	4
4	G	3	5	13	4	3	5
5	C	4	6	12	5	4	6
6	T	5	7	11	6	5	0
7	A	6	8	0	7	0	0
8	G	7	9	0	8	0	0
9	C	8	10	0	9	0	0
10	T	9	11	0	10	0	0
11	A	10	12	6	11	0	12
12	G	11	13	5	12	11	13
13	T	12	0	4	13	12	0"""

EXPECTED_CT_LINES = [
    [
        "1	C	0	2	0	1	0	0",
        "2	A	1	3	0	2	0	0",
        "3	T	2	4	0	3	0	0",
        "4	G	3	5	0	4	0	5",
        "5	C	4	6	12	5	4	6",
        "6	T	5	7	11	6	5	0",
        "7	A	6	8	0	7	0	0",
        "8	G	7	9	0	8	0	0",
        "9	C	8	10	0	9	0	0",
        "10	T	9	11	0	10	0	0",
        "11	A	10	12	6	11	0	12",
        "12	G	11	13	5	12	11	13",
        "13	T	12	0	0	13	12	0",
    ],
    [
        "1	C	0	2	0	1	0	0",
        "2	A	1	3	0	2	0	0",
        "3	T	2	4	0	3	0	4",
        "4	G	3	5	13	4	3	5",
        "5	C	4	6	12	5	4	6",
        "6	T	5	7	11	6	5	0",
        "7	A	6	8	0	7	0	0",
        "8	G	7	9	0	8	0	0",
        "9	C	8	10	0	9	0	0",
        "10	T	9	11	0	10	0	0",
        "11	A	10	12	6	11	0	12",
        "12	G	11	13	5	12	11	13",
        "13	T	12	0	4	13	12	0",
    ],
]


def test_parse_ct_to_folds() -> None:
    EXPECTED_HEADINGS = ["13\tdG = 0.892\tstdin", "13\tdG = 1.407\tstdin"]
    actual_headings, actual_lines = RNAfolder.parse_ct_to_folds(SAMPLE_CT)

    assert actual_headings == EXPECTED_HEADINGS
    assert actual_lines == EXPECTED_CT_LINES


@patch("fealden._unafold.open")
@patch("fealden._unafold.subprocess.Popen")
def test_collect_unafold_ct(mock_run: MagicMock, mock_open: MagicMock) -> None:
    _ = RNAfolder.collect_unafold_ct("CATGCTAGCTAGT")
    mock_run.assert_called_once()
    mock_open.assert_called_once()


def test_return_basepair() -> None:
    EXPECTED_PAIR_NUMBER = 13
    actual = RNAfolder.return_basepair(4, EXPECTED_CT_LINES[1])
    assert actual == EXPECTED_PAIR_NUMBER
