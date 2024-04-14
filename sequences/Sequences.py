from common import Chromosome_Info, Chromosome_Info_For_Simple
import math
import torch
import pyfaidx
import numpy as np
from sequences.toonehot import toonehot, read_n_values_from_dna

class Sequences():

    def __init__(self, fai_path : str):
        file_fai = open(fai_path)
        self.start_chr_positions = [None] * 24
        self.chr_bound = [None] * 24
        for line in file_fai:
            splitted = line.split()
            chr_num = splitted[0][3:]
            if chr_num == 'X':
                chr_num = 23
            elif chr_num == 'Y':
                chr_num = 24

            chr_num = int(chr_num)
            lenght = int(splitted[1])
            row_length = int(splitted[3])
            with_sep_row_length = int(splitted[4])
            wit_sep_lenght = lenght + math.ceil(lenght / row_length) * (with_sep_row_length - row_length)

            self.start_chr_positions[chr_num-1]=Chromosome_Info_For_Simple(chr_num, int(splitted[2]), lenght, wit_sep_lenght, row_length, with_sep_row_length)
            self.chr_bound[chr_num-1]=Chromosome_Info(chr_num, int(splitted[2]), lenght)

    def get_location_in_fasta(self, chr_info : Chromosome_Info, start, WINDOW_SIZE : int):
        assert start >= 0, 'Start position must be non-negative'

        diff = chr_info.with_sep_row_length - chr_info.row_length
        start_pos = chr_info.start_pos + start +  start // chr_info.row_length * diff
        lengh_window =  WINDOW_SIZE + (WINDOW_SIZE + start % chr_info.row_length) // chr_info.row_length * diff
        end_pos = start_pos + lengh_window

        assert end_pos <= chr_info.start_pos + chr_info.wit_sep_lenght , 'Make the starting position smaller or choose a different chromosome'

        return start_pos, lengh_window

    def get_name(self):
        return "dnaseq"

class DNASequence(Sequences):

    def __init__(self, fai_path : str, fa_path : str):
        super().__init__(
            fai_path,
        )
        self.fa_path = fa_path

    def to_one_hot(self, lines : str, WINDOW_SIZE : int):
        G = torch.zeros((WINDOW_SIZE,))
        T = torch.zeros((WINDOW_SIZE,))
        A = torch.zeros((WINDOW_SIZE,))
        C = torch.zeros((WINDOW_SIZE,))
        
        j = 0
        for i in lines:
            if i == 'G':
                G[j] = 1
            elif i == 'T':
                T[j] = 1
            elif i == 'A':
                A[j] = 1
            elif i == 'C':
                C[j] = 1
            elif i != 'N':
                continue
            j += 1
        
        return torch.stack([G, T, A, C])

    def get_lines(self, chr: int, start: int, end : int, WINDOW_SIZE: int):
            file_fa = open(self.fa_path)

            start_pos, lengh_window = self.get_location_in_fasta(self.start_chr_positions[chr], start, WINDOW_SIZE)

            file_fa.seek(start_pos, 0)
            lines = file_fa.read(lengh_window).replace("\n", "")
            assert len(lines) == WINDOW_SIZE, 'MISTAKE'
                

            return self.to_one_hot(lines, WINDOW_SIZE)

class BlankSequence(Sequences):

    def __init__(self, fai_path : str):
        super().__init__(
            fai_path,
        )

    def get_lines(self, chr: int, start: int, end : int, WINDOW_SIZE: int):
        return

    def get_name(self):
        return "blankseq"

class DNASequenceWithFasta(Sequences):

    def __init__(self, fai_path : str, fa_path : str):
        super().__init__(
            fai_path,
        )
        self.fa_path = fa_path
        self.genome = pyfaidx.Fasta(self.fa_path)

    def get_lines(self, chr: int, start: int, end : int, WINDOW_SIZE: int):
            cur_chr = "chr%d" % (chr + 1)
            if (chr == 22):
                cur_chr = "chrX"
            if (chr == 23):
                cur_chr = "chrY" 
            seq = self.genome[cur_chr][start:start+WINDOW_SIZE].seq
            return torch.tensor(toonehot(seq, {"A" : 0, "C" : 1, "G" : 2, "T": 3}))

class DNASequenceBase(Sequences):
    def __init__(self, fai_path : str, fa_path : str):
        super().__init__(
            fai_path,
        )
        self.fa_path = fa_path

    def get_lines(self, chr: int, start: int, end : int, WINDOW_SIZE: int):
            cur_chr = "chr%d" % (chr + 1)
            if (chr == 22):
                cur_chr = "chrX"
            if (chr == 23):
                cur_chr = "chrY" 
            start_pos, lengh_window = self.get_location_in_fasta(self.start_chr_positions[chr], start, WINDOW_SIZE)
            return torch.tensor(read_n_values_from_dna(self.fa_path, lengh_window, start_pos, np.array([0, 1, 2, 3])))

if __name__ == '__main__':
    dna_seq = DNASequence("helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                    "helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dna_seq.get_name()