import numpy as np
import pandas as pd

#get_neighber
def get_ID(i, j, k, L):
	ID = i + j*L + k*L**2
	return ID
 
def get_NN1(i, j, k, L):
	NN1 = np.array([[i-1, j, k],[i+1, j, k],[i, j-1,k],[i,j+1,k],[i,j,k-1],[i,j,k+1]])
	NN1 = NN1%L
	NN1_0 = NN1[0,0] + NN1[0,1]*L + NN1[0,2]*L**2
	NN1_1 = NN1[1,0] + NN1[1,1]*L + NN1[1,2]*L**2
	NN1_2 = NN1[2,0] + NN1[2,1]*L + NN1[2,2]*L**2
	NN1_3 = NN1[3,0] + NN1[3,1]*L + NN1[3,2]*L**2
	NN1_4 = NN1[4,0] + NN1[4,1]*L + NN1[4,2]*L**2
	NN1_5 = NN1[5,0] + NN1[5,1]*L + NN1[5,2]*L**2
	return NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5
