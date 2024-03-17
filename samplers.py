import sys
from characteristics import СharacteristicBigWig, Characteristics
from sequences import Sequences, DNASequence
from torch.utils.data import Dataset, DataLoader
from coordinates import RandomCoordinates
from common import Chromosome_Info
import numpy as np

class _DnaDataset(Dataset):
    def __init__(self, window_size : int, dna_seq : Sequences, big_wigs : Characteristics, seed = 3):
        super(_DnaDataset, self).__init__()
        chrom_info = [Chromosome_Info(0, 181390, 182390)] 
        self.window_size = window_size
        self.dna_seq = dna_seq
        self.big_wigs = big_wigs
        self.coord = RandomCoordinates(chrom_info, self.window_size, [i for i in range(1)], seed)

    def __getitem__(self, index):
        chr_num, start_ps, end_ps = self.coord.get_next_coord()
        seq = self.dna_seq.get_lines(chr_num, start_ps, end_ps, self.window_size)
        bw = self.big_wigs.get_lines(chr_num, start_ps, self.window_size)
        return {self.big_wigs.get_name() : bw, self.dna_seq.get_name() : seq}

    def __len__(self):
        return sys.maxsize
        


if __name__ == '__main__':
    chr_big_wig = СharacteristicBigWig("/home/ojpochemy/dna-loader/resfold", 1)
    dna_seq = DNASequence("/home/ojpochemy/dna-loader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                        "/home/ojpochemy/dna-loader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dna_dataset = _DnaDataset(3, dna_seq, chr_big_wig, 4)
    datal = DataLoader(dna_dataset, num_workers = 2, batch_size = 3)

    # print(dna_seq.get_lines(0, 1900000, 1901000, 1000))


    j = 0
    for i in datal:
        print(i)
        # print(np.shape(i))
        # print(len(i))
        # print("bw", len(i["bw"]))
        # print("seq", len(i["dnaseq"]))
        # print("seq", i["dnaseq"][0])
        break