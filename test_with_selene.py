from characteristics import СharacteristicBigWig, Characteristics
from sequences import DNASequence, BlankSequence, DNASequenceWithFasta
from samplers import _DnaDataset
from torch.utils.data import DataLoader
import time

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

if __name__ == '__main__':
    window_len = 100000
    batch_size = 16
    num_of_iterations = 10
    num_worker = 1

    chr_big_wig = СharacteristicBigWig("/home/ojpochemy/dnaloader/resfold", 1)
    dna_seq = DNASequenceWithFasta("/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai", 
                    "/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    dna_dataset = _DnaDataset(window_len, dna_seq, chr_big_wig, 4)
    datal = DataLoader(dna_dataset, num_workers = num_worker, batch_size = batch_size)

    j = 0
    start_time = time.time() 
    for cur_res in datal:
        j += 1
        if (j == num_of_iterations):
            break
    end_time = time.time() 
    print("end time my dna-loader ", end_time - start_time)


    bigwig_file_1 = "/home/ojpochemy/SamplerBigWig/foldbigwig/interval.all.obs_1.bw" 
    bigwig_file_2 = "/home/ojpochemy/SamplerBigWig/foldbigwig/interval.all.obs_2.bw" 
    bigwig_file_3 = "/home/ojpochemy/SamplerBigWig/foldbigwig/interval.all.obs_3.bw" 
    bigwig_file_4 = "/home/ojpochemy/SamplerBigWig/foldbigwig/interval.all.obs_4.bw" 
    bigwig_file_5 = "/home/ojpochemy/SamplerBigWig/foldbigwig/interval.all.obs_5.bw" 
    ref_file = "/home/ojpochemy/dnaloader/sequences/helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    target = GenomicSignalFeatures(
        input_paths=[bigwig_file_1, bigwig_file_2, bigwig_file_3, bigwig_file_4, bigwig_file_5],
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
        test_holdout=['chr8'],
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

    j = 0
    start_time = time.time() 
    for seq, x in loader:
        # print(seq.shape)
        # print(seq)
        # print(x.shape)
        j += 1
        if (j == num_of_iterations):
            break
    end_time = time.time() 
    print("end time selene ", end_time - start_time)