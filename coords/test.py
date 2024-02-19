from Coordinates import HandMadeCoordinates
from Coordinates import RandomCoordinates
from common import Chromosome_Info
import unittest

class TestReadLinesWithTwoNumbers(unittest.TestCase):

    chrom_info = [Chromosome_Info(1, 100, 200),
                        Chromosome_Info(2, 300, 400),
                        Chromosome_Info(3, 700, 800),
                        Chromosome_Info(4, 1100, 100)]  

    def test_hand_made_coord(self):
        coords = HandMadeCoordinates(self.chrom_info, 100, [i for i in range(4)], "coords/test/handMade.txt")
        self.assertEqual(coords.get_next_coord(), (1, 23, 43))
        self.assertEqual(coords.get_next_coord(), (0, 100, 120))
        self.assertEqual(coords.get_next_coord(), (2, 10, 30))
        self.assertEqual(coords.get_next_coord(), (1, 2, 23))
        self.assertEqual(coords.get_next_coord(), (1, 23, 43))

    def test_wrong_hand_made(self):
        with self.assertRaises(Exception):
            coords = HandMadeCoordinates(self.chrom_info, 100, [i for i in range(4)], "coords/test/wrongHandMade.txt")
        
    def test_random_coord(self):  
        coords = RandomCoordinates(self.chrom_info, 100, [i for i in range(4)], 3)
        first_c = coords.get_next_coord()
        second_c = coords.get_next_coord()
        self.assertNotEqual(first_c, second_c)

if __name__ == '__main__':
    unittest.main()