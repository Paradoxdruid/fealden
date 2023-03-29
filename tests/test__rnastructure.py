import os
import sys
from unittest import mock

RNA_PATH = os.getenv("RNASTRUCTURE", "/home/andrew/RNAstructure")
sys.path.append(os.path.join(RNA_PATH, "exe"))


def test__init__() -> None:
    with mock.patch.dict(
        "sys.modules",
        {
            "RNAstructure": mock.MagicMock(),
            "_RNAstructure_wrap": mock.MagicMock(),
            "RNAstructure_wrap": mock.MagicMock(),
            "Error_Handling": mock.MagicMock(),
        },
    ):
        from fealden._rnastructure import RNAfolder, RNAstructure  # type: ignore

        _ = RNAfolder("CACGTGGTGCAC")

        RNAstructure.RNA.fromString.assert_called_once_with(
            "CACGTGGTGCAC", backbone="dna"
        )
