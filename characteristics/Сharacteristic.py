import os
import pyBigWig
import math

class Сharacteristic():

    def __init__(self, path :  str):
        self.files = os.listdir(path)


class СharacteristicBigWig():

    def __init__(self, path :  str):
        super().__init__(
            path
        )
        # you can select any value from the segment -> [1...24]
        size_of_bw = 24

        bw_files = []
        for i in range(self.files):
            bw_files.append(pyBigWig.open("interval.all.obs.bw"))

    



