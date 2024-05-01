from sequences.Sequences import Sequences, DNASequence
import sys

sys.path.append("/Users/yulkorolenko/Desktop/Diploma/dna-loader")

# def test_chr_bounds():
#     seq = Sequences("helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai")
#     print(seq.chr_bound[0])


def test_dna_sequence():
    dna_seq = DNASequence("helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
                          "helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")

    print(dna_seq.get_lines(0, 0, 10, 10))
