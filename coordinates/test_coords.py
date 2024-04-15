import pytest
from coordinates.Coordinates import RandomCoordinates
from coordinates.Coordinates import HandMadeCoordinates
import sys

sys.path.append("/Users/yulkorolenko/Desktop/Diploma/dna-loader")

chrom_info = [(100, 200), (300, 400), (700, 800), (1100, 100)]


def test_hand_made_coord():
    coords = HandMadeCoordinates(
        chrom_info, 100, [
            i for i in range(4)], "helper/handMade.txt")
    assert coords.get_next_coord() == (1, 23, 43)
    assert coords.get_next_coord() == (0, 100, 120)
    assert coords.get_next_coord() == (2, 10, 30)
    assert coords.get_next_coord() == (1, 2, 23)
    assert coords.get_next_coord() == (1, 23, 43)


def test_wrong_hand_made():
    with pytest.raises(Exception):
        coords = HandMadeCoordinates(
            chrom_info, 100, [
                i for i in range(4)], "helper/wrongHandMade.txt")


def test_random_coord():
    coords = RandomCoordinates(chrom_info, 100, [i for i in range(4)], 3)
    first_c = coords.get_next_coord()
    second_c = coords.get_next_coord()
    assert first_c != second_c
