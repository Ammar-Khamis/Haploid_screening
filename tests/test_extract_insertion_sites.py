import pytest
from Extract_Insertion_Sites import parse_cigar, calculate_5prime_position


@pytest.mark.parametrize(
    "cigar, expected",
    [
        ("10M5I5M", [(10, "M"), (5, "I"), (5, "M")]),
        ("5S10M2D", [(5, "S"), (10, "M"), (2, "D")]),
    ],
)
def test_parse_cigar(cigar, expected):
    assert parse_cigar(cigar) == expected


def test_calculate_5prime_position_forward_insertion():
    cigar = parse_cigar("10M5I5M")
    assert calculate_5prime_position(100, cigar, False) == 99


def test_calculate_5prime_position_reverse_insertion():
    cigar = parse_cigar("10M5I5M")
    assert calculate_5prime_position(100, cigar, True) == 113


def test_calculate_5prime_position_forward_softclip_del():
    cigar = parse_cigar("5S10M2D")
    assert calculate_5prime_position(50, cigar, False) == 49


def test_calculate_5prime_position_reverse_softclip_del():
    cigar = parse_cigar("5S10M2D")
    assert calculate_5prime_position(50, cigar, True) == 60

