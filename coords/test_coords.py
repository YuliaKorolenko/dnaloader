import sys

sys.path.append("/Users/yulkorolenko/Desktop/Diploma/dna-loader")
from coords.Coordinates import HandMadeCoordinates
from coords.Coordinates import RandomCoordinates
from common import Chromosome_Info
import pytest

chrom_info = [Chromosome_Info(1, 100, 200),
                    Chromosome_Info(2, 300, 400),
                    Chromosome_Info(3, 700, 800),
                    Chromosome_Info(4, 1100, 100)]  

def test_hand_made_coord():
    coords = HandMadeCoordinates(chrom_info, 100, [i for i in range(4)], "helper/handMade.txt")
    assert coords.get_next_coord() == (1, 23, 43)
    assert coords.get_next_coord() == (0, 100, 120)
    assert coords.get_next_coord() == (2, 10, 30)
    assert coords.get_next_coord() == (1, 2, 23)
    assert coords.get_next_coord() == (1, 23, 43)

def test_wrong_hand_made():
    with pytest.raises(Exception):
        coords = HandMadeCoordinates(chrom_info, 100, [i for i in range(4)], "helper/wrongHandMade.txt")
    
def test_random_coord():  
    coords = RandomCoordinates(chrom_info, 100, [i for i in range(4)], 3)
    first_c = coords.get_next_coord()
    second_c = coords.get_next_coord()
    assert first_c != second_c