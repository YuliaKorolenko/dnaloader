import unittest
from Sequences import Sequences, DNASequence

class TestSequences(unittest.TestCase):

    def test_chr_bounds(self):
        seq = Sequences("sequences/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai")
        print(seq.chr_bound[0])

    def test_dna_sequence(self):
        dna_seq = DNASequence("sequences/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                           "sequences/test/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
        
        print(dna_seq.get_lines(0, 0, 10, 10))
        



if __name__ == '__main__':
    unittest.main()