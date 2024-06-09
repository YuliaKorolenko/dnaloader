import sys
from .characteristics import Characteristic
from .sequences import Sequences
from torch.utils.data import Dataset
from .coordinates import RandomCoordinates
import numpy as np
from typing import List

class DnaDataset(Dataset):
    def __init__(self, window_size: int, dna_seq: Sequences,
                 char_list: List[Characteristic], seed=3):
        super(DnaDataset, self).__init__()
        self.window_size = window_size
        self.dna_seq = dna_seq
        self.char_list = char_list
        self.coord = RandomCoordinates(
            self.dna_seq.get_chr_bounds(), self.window_size, seed)

class DnaDataset1d(DnaDataset):
    def __init__(self, window_size: int, dna_seq: Sequences,
                 char_list: List[Characteristic], seed=3):
        super().__init__(
            window_size = window_size,
            dna_seq = dna_seq,
            char_list = char_list,
            seed = seed
        )

    def __getitem__(self, index):
        chr_num, start_ps, end_ps = self.coord.get_next_coord()
        seq = self.dna_seq.get_lines(
            chr_num, start_ps, end_ps, self.window_size)
        answer = {}
        for char_class in self.char_list:
            answer[char_class.get_name()] = char_class.get_lines(
                chr_num, start_ps, self.window_size)
        if self.dna_seq.get_name() != "blankseq":
            answer[self.dna_seq.get_name()] = seq
        return answer

    def __len__(self):
        return sys.maxsize
    

class DnaDataset2d(DnaDataset):
    def __init__(self, window_size: int, dna_seq: Sequences,
                 char_list: List[Characteristic], seed=3):
        super().__init__(
            window_size = window_size,
            dna_seq = dna_seq,
            char_list = char_list,
            seed = seed
        )

    def __getitem__(self, index):
        chr_num_1, start_ps_1, end_ps_1 = self.coord.get_next_coord()
        chr_num_2, start_ps_2, end_ps_2 = self.coord.get_next_coord()

        seq_1 = self.dna_seq.get_lines(
            chr_num_1, start_ps_1, end_ps_1, self.window_size)
        seq_2 = self.dna_seq.get_lines(
            chr_num_2, start_ps_2, end_ps_2, self.window_size)
        answer = {}
        for char_class in self.char_list:
            answer[char_class.get_name()] = char_class.get_lines(
                chr_num_1, start_ps_1, start_ps_2, self.window_size)
        if self.dna_seq.get_name() != "blankseq":
            answer[self.dna_seq.get_name()] = np.array([seq_1, seq_2])
        return answer

    def __len__(self):
        return sys.maxsize
