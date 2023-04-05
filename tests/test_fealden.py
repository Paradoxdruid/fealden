import argparse
import random
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from fealden.fealden import Fealden, generate_sensor, main
from fealden.seed import Seed


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


def test_generate_sensor() -> None:
    with mock.patch("fealden.fealden.seed.random") as mock_random:
        mock_random.randint = random.Random(0).randint
        mock_random.choice = random.Random(0).choice

        EXPECTED_SENSOR = (
            "acttcgggacttgcttgaagcacgtgctattggtaccaatagtgagaagt,"
            "1.130352975773343,Graph 2,16,673.9453257075509,1114.7313270339453,"
            "1.6540382201829638,0.0,0.0,0.0,10.0,50,2,CACGTG"
        )

        actual = generate_sensor(
            Seed(
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
            ),
            "CACGTG",
            500,
            1,
        )
        assert repr(actual[0]) == EXPECTED_SENSOR


@mock.patch("fealden.fealden.Fealden")
@mock.patch("argparse.ArgumentParser.parse_args")
def test__main__(mock_arg: mock.Mock, mock_fealden: mock.Mock) -> None:
    mock_arg.return_value = argparse.Namespace(
        recSeq="CACGTG",
        bindingState=1,
        ms=50,
        sps=500,
        interactive=None,
        out="test.csv",
        v=None,
    )
    main()
    mock_fealden.assert_called_once_with("cacgtg", 1, 50, 500, None, "test.csv")
