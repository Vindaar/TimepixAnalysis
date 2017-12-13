# a combination of some useful functions dealing with InGrid related helper functions
# mostly related to plotting data, which was produced by Nim

import numpy as np

def readFileColumn(filename):
    lines = open(filename, "r").readlines()
    data = []
    for line in lines:
        if "#" not in line and "nan" not in line:
            data.append(float(line))
    return np.asarray(data)
    
def readFadcSpectrumFile(filename):
    return readFileColumn(filename)

def readHitsFile(filename):
    return readFileColumn(filename)

def readTotsFile(filename):
    return readFileColumn(filename)
