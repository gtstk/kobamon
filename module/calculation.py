import pandas as pd
import numpy as np

def calc_E0(spin, ID, NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5, J1, H):
	sum1 = spin[NN1_0]+spin[NN1_1]+spin[NN1_2]+spin[NN1_3]+spin[NN1_4]+spin[NN1_5]
	exc1 = J1*(np.dot(spin[ID], sum1))
	dem = np.dot(spin[ID], H)
	E0 = -exc1 -dem
	return E0

def calc_E1(next_spin, spin, NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5, J1, H):
    sum1 = spin[NN1_0]+spin[NN1_1]+spin[NN1_2]+spin[NN1_3]+spin[NN1_4]+spin[NN1_5]
    exc1 = J1*(np.dot(next_spin, sum1))
    dem = np.dot(next_spin, H)
    E1 = -exc1 -dem
    return E1

def calc_change(E0, E1, T):
	P0 = (np.exp(-E0/T))
	P1 = (np.exp(-E1/T))
	change = P1/(P0+P1)
	return change