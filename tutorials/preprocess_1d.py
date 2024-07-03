from dnaloader.characteristics import Characteristic1dPreprocess
import h5py
import os

if __name__ == '__main__':
    preprocess_class = Characteristic1dPreprocess(prev_path="helper/filepaths.txt",
                                                  res_path="hi.hdf5")
    preprocess_class.preprocess()

    file = h5py.File('hi.hdf5', 'r')
    keys = list(file.keys())

    p = " "
    for key in keys:
        dataset = file[key] 
        print(key) 
        print(dataset)
        for i in dataset:
            print(p * 3,i)
            print(dataset[i][:10])
    file.close()