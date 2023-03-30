import os
import sys
from unittest import mock

from dotenv import load_dotenv

load_dotenv()

RNA_PATH = os.getenv("RNASTRUCTURE", "/home")
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
