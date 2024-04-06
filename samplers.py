import sys
from characteristics import СharacteristicBigWig, Characteristics, СharacteristicBigWigCSR
from sequences import Sequences, DNASequence, BlankSequence
from torch.utils.data import Dataset, DataLoader
from coordinates import RandomCoordinates
from common import Chromosome_Info
import numpy as np

class _DnaDataset(Dataset):
    def __init__(self, window_size : int, dna_seq : Sequences, big_wigs : Characteristics, seed = 3):
        super(_DnaDataset, self).__init__()
        self.window_size = window_size
        self.dna_seq = dna_seq
        self.big_wigs = big_wigs
        self.coord = RandomCoordinates(self.big_wigs.get_chr_info(), self.window_size, seed)

    def __getitem__(self, index):
        chr_num, start_ps, end_ps = self.coord.get_next_coord()
        # print(chr_num, start_ps, self.window_size)
        seq = self.dna_seq.get_lines(chr_num, start_ps, end_ps, self.window_size)
        bw = self.big_wigs.get_lines(chr_num, start_ps, self.window_size)
        if self.dna_seq.get_name() == "blankseq":
            return {self.big_wigs.get_name() : bw}
        return {self.big_wigs.get_name() : bw, self.dna_seq.get_name() : seq}

    def __len__(self):
        return sys.maxsize
        


if __name__ == '__main__':
    chr_big_wig = СharacteristicBigWigCSR("/home/ojpochemy/dnaloader/resfold", 1)
    dna_seq = DNASequence("/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                        "/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dna_dataset = _DnaDataset(3, dna_seq, chr_big_wig, 4)
    datal = DataLoader(dna_dataset, num_workers = 2, batch_size = 1)
    for k in datal:
        print(k)