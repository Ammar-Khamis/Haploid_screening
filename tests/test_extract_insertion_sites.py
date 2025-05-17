import sys
import importlib
import importlib.util
from pathlib import Path
import pytest

# Get the path to the Haploid_screening_ONT directory (parent of parent)
MODULE_DIR = Path(__file__).resolve().parent.parent.parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))
importlib.invalidate_caches()

# Try to import directly first
try:
    from Extract_Insertion_Sites import parse_cigar, calculate_5prime_position
except ImportError:
    # If direct import fails, try loading from file
    module_path = str(MODULE_DIR / "Extract_Insertion_Sites.py")
    spec = importlib.util.spec_from_file_location("Extract_Insertion_Sites", module_path)
    if spec is not None and hasattr(spec, 'loader') and spec.loader is not None:
        Extract_Insertion_Sites_module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(Extract_Insertion_Sites_module)
        parse_cigar = Extract_Insertion_Sites_module.parse_cigar
        calculate_5prime_position = Extract_Insertion_Sites_module.calculate_5prime_position
    else:
        raise ImportError(f"Could not import Extract_Insertion_Sites module from {module_path}")

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
