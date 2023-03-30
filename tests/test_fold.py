from fealden.fold import Fold


def test_Fold() -> None:
    INPUT_BPS = [
        [1, 0],
        [2, 0],
        [3, 0],
        [4, 0],
        [5, 0],
        [6, 0],
        [7, 0],
        [8, 0],
        [9, 25],
        [10, 24],
        [11, 0],
        [12, 22],
        [13, 21],
        [14, 20],
        [15, 0],
        [16, 0],
        [17, 0],
        [18, 0],
        [19, 0],
        [20, 14],
        [21, 13],
        [22, 12],
        [23, 0],
        [24, 10],
        [25, 9],
        [26, 0],
        [27, 42],
        [28, 41],
        [29, 40],
        [30, 39],
        [31, 38],
        [32, 37],
        [33, 0],
        [34, 0],
        [35, 0],
        [36, 0],
        [37, 32],
        [38, 31],
        [39, 30],
        [40, 29],
        [41, 28],
        [42, 27],
        [43, 0],
        [44, 0],
        [45, 0],
        [46, 0],
        [47, 0],
        [48, 0],
        [49, 0],
        [50, 0],
    ]
    INPUT_DELTAG = -4.155

    EXPECTED = "Fold: head=SSNode: length=8,sequence=, deltaG=-4.155,\
 conc=1114.7313270339453, foldData=[[1, 0], [2, 0], [3, 0], [4, 0], [5, 0], [6, 0],\
 [7, 0], [8, 0], [9, 25], [10, 24], [11, 0], [12, 22], [13, 21], [14, 20], [15, 0],\
 [16, 0], [17, 0], [18, 0], [19, 0], [20, 14], [21, 13], [22, 12], [23, 0], [24, 10],\
 [25, 9], [26, 0], [27, 42], [28, 41], [29, 40], [30, 39], [31, 38], [32, 37],\
 [33, 0], [34, 0], [35, 0], [36, 0], [37, 32], [38, 31], [39, 30], [40, 29], [41, 28],\
 [42, 27], [43, 0], [44, 0], [45, 0], [46, 0], [47, 0], [48, 0], [49, 0], [50, 0]],\
 ptrList=[SSNode: length=8,sequence=, SSNode: length=8,sequence=, SSNode:\
 length=8,sequence=, SSNode: length=8,sequence=, SSNode: length=8,sequence=,\
 SSNode: length=8,sequence=, SSNode: length=8,sequence=, SSNode: length=8,sequence=,\
 DSNode: length=2,sequence=, DSNode: length=2,sequence=, SSNode: length=1,sequence=,\
 DSNode: length=3,sequence=, DSNode: length=3,sequence=, DSNode: length=3,sequence=,\
 SSNode: length=5,sequence=, SSNode: length=5,sequence=, SSNode: length=5,sequence=,\
 SSNode: length=5,sequence=, SSNode: length=5,sequence=, DSNode: length=3,sequence=,\
 DSNode: length=3,sequence=, DSNode: length=3,sequence=, SSNode: length=1,sequence=,\
 DSNode: length=2,sequence=, DSNode: length=2,sequence=, SSNode: length=1,sequence=,\
 DSNode: length=6,sequence=, DSNode: length=6,sequence=, DSNode: length=6,sequence=,\
 DSNode: length=6,sequence=, DSNode: length=6,sequence=, DSNode: length=6,sequence=,\
 SSNode: length=4,sequence=, SSNode: length=4,sequence=, SSNode: length=4,sequence=,\
 SSNode: length=4,sequence=, DSNode: length=6,sequence=, DSNode: length=6,sequence=,\
 DSNode: length=6,sequence=, DSNode: length=6,sequence=, DSNode: length=6,sequence=,\
 DSNode: length=6,sequence=, SSNode: length=8,sequence=, SSNode: length=8,sequence=,\
 SSNode: length=8,sequence=, SSNode: length=8,sequence=, SSNode: length=8,sequence=,\
 SSNode: length=8,sequence=, SSNode: length=8,sequence=, SSNode: length=8,sequence=],\
 recSeq={'start': 21, 'end': 27}, recSeqState=0"

    actual = Fold(INPUT_BPS, INPUT_DELTAG, {"start": 21, "end": 27})

    assert repr(actual) == EXPECTED
