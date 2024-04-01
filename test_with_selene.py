from characteristics import СharacteristicBigWig, Characteristics, СharacteristicBigWigCSR
from sequences import DNASequence, BlankSequence, DNASequenceWithFasta
from samplers import _DnaDataset
from torch.utils.data import DataLoader
import time
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

import os

import tabix
import pyfaidx
import pyBigWig
import pkg_resources

from functools import wraps

from selene_sdk.targets import Target
from selene_sdk.sequences import Genome
from selene_sdk.samplers import RandomPositionsSampler
from selene_sdk.samplers import SamplerDataLoader
import threading
import matplotlib.pyplot as plt

class GenomicSignalFeatures(Target):
    """
    #Accept a list of cooler files as input.
    """

    def __init__(self, input_paths, features, shape, blacklists=None, blacklists_indices=None,
                 replacement_indices=None, replacement_scaling_factors=None):
        """
        Constructs a new `GenomicFeatures` object.
        """
        self.input_paths = input_paths
        self.initialized = False
        self.blacklists = blacklists
        self.blacklists_indices = blacklists_indices
        self.replacement_indices = replacement_indices
        self.replacement_scaling_factors = replacement_scaling_factors

        self.n_features = len(features)
        self.feature_index_dict = dict(
            [(feat, index) for index, feat in enumerate(features)])
        self.shape = (len(input_paths), *shape)

    def get_feature_data(self, chrom, start, end, nan_as_zero=True, feature_indices=None):
        if not self.initialized:
            self.data = [pyBigWig.open(path) for path in self.input_paths]
            if self.blacklists is not None:
                self.blacklists = [tabix.open(blacklist) for blacklist in self.blacklists]
            self.initialized = True

        if feature_indices is None:
            feature_indices = np.arange(len(self.data))

        wigmat = np.zeros((len(feature_indices), end - start), dtype=np.float32)
        for i in feature_indices:
            try:
                wigmat[i, :] = self.data[i].values(chrom, start, end, numpy=True)
            except:
                print(chrom, start, end, self.input_paths[i], flush=True)
                raise

        if self.blacklists is not None:
            if self.replacement_indices is None:
                if self.blacklists_indices is not None:
                    for blacklist, blacklist_indices in zip(self.blacklists, self.blacklists_indices):
                        for _, s, e in blacklist.query(chrom, start, end):
                            wigmat[blacklist_indices, np.fmax(int(s) - start, 0): int(e) - start] = 0
                else:
                    for blacklist in self.blacklists:
                        for _, s, e in blacklist.query(chrom, start, end):
                            wigmat[:, np.fmax(int(s) - start, 0): int(e) - start] = 0
            else:
                for blacklist, blacklist_indices, replacement_indices, replacement_scaling_factor in zip(
                        self.blacklists, self.blacklists_indices, self.replacement_indices,
                        self.replacement_scaling_factors):
                    for _, s, e in blacklist.query(chrom, start, end):
                        wigmat[blacklist_indices, np.fmax(int(s) - start, 0): int(e) - start] = wigmat[
                                                                                                replacement_indices,
                                                                                                np.fmax(int(s) - start,
                                                                                                                                                                                                       0): int(
                                                                                                    e) - start] * replacement_scaling_factor

        if nan_as_zero:
            wigmat[np.isnan(wigmat)] = 0
        return wigmat

def get_my_time(window_len, chr_big_wig, dna_seq,  batch_size, num_of_iterations, num_worker):
    dna_dataset = _DnaDataset(window_len, dna_seq, chr_big_wig, 59)
    datal = DataLoader(dna_dataset, num_workers = num_worker, batch_size = batch_size)

    j = 0
    start_time = time.time() 
    for cur_res in datal:
        j += 1
        if (j == num_of_iterations):
            break
    end_time = time.time() 
    print("end time my dna-loader ", end_time - start_time)
    del datal
    return end_time - start_time

def get_selene_time(num_of_iterations):
    j = 0
    start_time = time.time() 
    for seq, x in loader:
        j += 1
        if (j == num_of_iterations):
            break
    end_time = time.time() 
    print("end time selene ", end_time - start_time)
    return end_time - start_time


if __name__ == '__main__':

    chr_big_wig = СharacteristicBigWigCSR("/home/ojpochemy/dnaloader/resfold", 1)
    dna_seq = DNASequenceWithFasta("/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                    "/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")

    bigwig_file_1 = "/home/ojpochemy/SamplerBigWig/foldbigwig0/interval.all.obs_1.bw" 
    bigwig_file_2 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_2.bw" 
    bigwig_file_3 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_3.bw" 
    bigwig_file_4 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_4.bw" 
    bigwig_file_5 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_5.bw" 
    bigwig_file_6 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_6.bw" 
    bigwig_file_7 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_7.bw" 
    bigwig_file_8 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_8.bw" 
    bigwig_file_9 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_9.bw" 
    bigwig_file_10 = "/home/ojpochemy/SamplerBigWig/foldbigwig1/interval.all.obs_10.bw" 
    ref_file = "/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    my_time = []
    selene_time = []
    window_sizes = []
    for i in range (0, 1):
        window_len = 4000
        batch_size = 1 
        window_sizes.append(batch_size)
        num_of_iterations = 4000
        num_worker = 1
        target = GenomicSignalFeatures(
            input_paths=[bigwig_file_1, bigwig_file_2, bigwig_file_3, bigwig_file_4, bigwig_file_5,
            bigwig_file_6, bigwig_file_7, bigwig_file_8, bigwig_file_9, bigwig_file_10],
            features=[], # doing nothing here
            shape=(window_len,),
        )

        ref_genome = Genome(
            input_path=ref_file,
            blacklist_regions="hg38",
        )

        sampler = RandomPositionsSampler(
            reference_sequence=ref_genome,
            target=target,
            features=['r1000'], # again doing nothing here
            test_holdout=['chr1', 'chr2', 'chr3'],
            validation_holdout=['chr9'],
            sequence_length=window_len,
            center_bin_to_predict=window_len,
            position_resolution=1,
            random_shift=0,
            random_strand=False,
            seed=239
        )

        sampler.mode = "train"

        loader = SamplerDataLoader(sampler, num_workers=num_worker, batch_size=batch_size)

        t1 = get_my_time(window_len, chr_big_wig, dna_seq, batch_size, num_of_iterations, num_worker)
        my_time.append(t1)
        t2 = get_selene_time(num_of_iterations)
        selene_time.append(t2)

    df = pd.DataFrame({
        'Batch size': window_sizes,
        'My loader': my_time,
        'Selene loader': selene_time
    })

    # # Устанавливаем стиль seaborn
    # sns.set(style="darkgrid")

    # # Создаем график с помощью Seaborn
    # sns.lineplot(data=df, x='Batch size', y='My loader', label='My loader')
    # sns.lineplot(data=df, x='Batch size', y='Selene loader', label='Selene loader')
    # plt.xlabel('Batch size')
    # plt.ylabel('Time in seconds')
    # plt.title('Time vs Batch size with window len 100000 and num of iteration 10')
    # plt.legend()


    # # Сохраняем график в файл "selenevsmy.png"
    # plt.savefig('selenevsmy4.png')
    # plt.show()
