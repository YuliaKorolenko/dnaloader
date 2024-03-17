from common import Chromosome_Info, Chromosome_Info_For_Simple
from typing import List
import sys
import random

class CoordinateBase():

    # Constructor
    def __init__(self, chr_bounds : List[Chromosome_Info], window_size : int, 
                 chromosomes : List[int]):
        self.chr_bounds = chr_bounds
        self.window_size = window_size
        self.chromosomes = chromosomes

class HandMadeCoordinates(CoordinateBase):

    def __init__(self, chr_bounds : List[Chromosome_Info], window_size : int,
                  chromosomes : List[int], file_path : str):
        super().__init__(
            chr_bounds,
            window_size, 
            chromosomes
        )
        self.file_path = file_path
        self.cur_pos = -1
            
        try:
            with open(file_path, "r") as file:
                self.numbers_list = [(int(t.split()[0]), int(t.split()[1]), int(t.split()[2]))
                                  for t in file.readlines()]
        except FileNotFoundError:
            raise Exception("The file could not be opened: %s" % file_path)
            
        # check boarders of chromosome
        for elem_in_list in self.numbers_list:
            if (chr_bounds[elem_in_list[0]].lenght <= elem_in_list[2]):
                raise RuntimeError("Произошла ошибка!")

    def get_next_coord(self):
        self.cur_pos += 1
        if (self.cur_pos >= len(self.numbers_list)):
            self.cur_pos = 0
        return self.numbers_list[self.cur_pos]
    

class RandomCoordinates(CoordinateBase):

    def __init__(self, chr_bounds : List[Chromosome_Info], window_size : int,
                  chromosomes : List[int], seed):
        super().__init__(
            chr_bounds,
            window_size, 
            chromosomes
        )
        random.seed(seed)
        self.state = random.getstate()
        self.full_length = 0
        for cur_chr in chr_bounds:
            self.full_length += cur_chr.lenght

    def get_correct_chr(self, rnd_start : int):
        for pos_i in range(len(self.chr_bounds)):
            cur_bound = self.chr_bounds[pos_i]
            if cur_bound.start_pos <= rnd_start < cur_bound.start_pos + cur_bound.lenght:
                return cur_bound, pos_i
            pos_i += 1
        return Chromosome_Info(-1, -1, -1), -1

    def get_next_coord(self):
        random.setstate(self.state)

        find_chr = False
        while(not find_chr):
            rnd_start = random.randint(0, self.full_length)
            self.state = random.getstate()
            rnd_bounds, number = self.get_correct_chr(rnd_start)
            if (rnd_start + self.window_size < rnd_bounds.start_pos + rnd_bounds.lenght):
                # return `chromosome num`, `start position` in chromosome and `end position` in chromosome
                return (number,
                        rnd_start - rnd_bounds.start_pos,
                        rnd_start - rnd_bounds.start_pos + self.window_size)
        





        

    


    


