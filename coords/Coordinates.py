from typing import List
from common import Chromosome_Info
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

        with open(file_path, "r") as file:
            self.numbers_list = [(int(t.split()[0]), int(t.split()[1])) for t in file.readlines()]
            
        # TO DO: check format (правильные границы хромосом)

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

    def get_next_coord(self):
        random.setstate(self.state)
        chr_rnd_number = random.randint(0, len(self.chromosomes) - 1)
        # TO DO: если одна хромосома больше другой, должны ли мы в нее в два раза чаще попадать?
        chr_rnd = self.chromosomes[chr_rnd_number]
        start_rnd = random.randint(0, self.chr_bounds[chr_rnd].lenght - self.window_size)
        self.state = random.getstate()
        return (chr_rnd, start_rnd)



        

    


    


