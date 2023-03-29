from fealden.seed import Seed


def test_Seed() -> None:
    actual = Seed(
        [
            "2 1 3 11 0",
            "3 2 4",
            "4 3 5 5 7",
            "5 4 4",
            "7 4 6",
            "6 7 9 9 11",
            "9 6 6",
            "11 6 2",
        ],
        "7",
        "CACGTG",
        1,
        "Graph 2",
        50,
    )

    EXPECTED = "Seed: name=Graph 2, head=SSNode: length=0,sequence=,\
        nodes={'1': 'SSNode: length=0,sequence=', '2': 'DSNode: length=-1,sequence=',\
 '3': 'SSNode: length=-1,sequence=', '11': 'SSNode: length=-1,sequence=', '0': 'None',\
 '4': 'DSNode: length=-1,sequence=', '5': 'SSNode: length=-1,sequence=', '7': 'SSNode:\
 length=-1,sequence=', '6': 'DSNode: length=-1,sequence=', '9': 'SSNode: length=-1,\
sequence='}, recNodeName=7,        recSeq=CACGTG, bindingState=1,        max_size=50"

    assert repr(actual) == EXPECTED
