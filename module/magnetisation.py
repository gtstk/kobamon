import numpy as np
import pandas as pd

def calc_mag(spin):
    mean = np.mean(spin, axis = 0)
    mag = np.linalg.norm(mean)
    return mag

def calc_mz(spin):
    mean = np.mean(spin, axis = 0)
    mz = mean[2]
    return mz