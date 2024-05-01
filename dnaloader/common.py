
from dataclasses import dataclass


@dataclass
class Chromosome_Info:
    number: int
    start_pos: int
    lenght: int


@dataclass
class Chromosome_Info_For_Simple:
    number: int
    start_pos: int
    lenght: int
    wit_sep_lenght: int
    row_length: int
    with_sep_row_length: int
