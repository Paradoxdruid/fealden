import random
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from fealden.fealden import Fealden


def test_Fealden() -> None:
    with mock.patch("fealden.fealden.seed.random") as mock_random:
        mock_random.randint = random.Random(0).randint
        mock_random.choice = random.Random(0).choice

        with TemporaryDirectory() as tmpdirname:
            my_file = Path(tmpdirname).joinpath("test-results.csv")
            Fealden("cacgtg", 1, 50, 500, False, str(my_file.resolve()))

            assert "aaacgcagtatgtggccatacacgtgatgtgagagacatccgttagttt" in str(
                my_file.read_text()
            )
