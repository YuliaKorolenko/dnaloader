import time
from torch.utils.data import DataLoader
from dnaloader.samplers import DnaDataset1d
from dnaloader.characteristics import Characteristic1d, Characteristic
from dnaloader.sequences import DNASequenceWithFasta

def get_my_time(track: Characteristic, chr: int, start: int, window: int):
    start_time = time.time()
    track.get_indixies(WINDOW_SIZE=window, start=start, chr=chr)
    return time.time() - start_time

if __name__ == '__main__':
    window_len = 1000000
    batch_size = 2 
    num_of_iterations = 1
    num_worker = 4
    # You can preprocess bigwigs as follows:
    # chr_big_wig = Characteristic1d(folder_res="/Users/yuliya.karalenka/Desktop/dnaloader/resfold", size_of_bw=2, path="/Users/yuliya.karalenka/Desktop/dnaloader/helper/filepaths.txt", type_of_loader="medium",)
    # This is an example for bigWig, which has already been pre-processed.
    chr_big_wig = Characteristic1d(folder_res="/Users/yuliya.karalenka/Desktop/dnaloader/resfold", size_of_bw=2, type_of_loader="medium")
    print(get_my_time(chr_big_wig, 1, 2, 1000000))
