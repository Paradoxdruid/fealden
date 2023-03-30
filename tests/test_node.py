from fealden.node import DSNode, SSNode


def test_DSNode_init() -> None:
    actual = DSNode(None, 3)

    assert repr(actual) == "DSNode: length=3,sequence="


def test_SSNode_init() -> None:
    actual = SSNode(None, 3)

    assert repr(actual) == "SSNode: length=3,sequence="
