import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from module import neighber
from module import calculation
from module import magnetisation

#init_parameter
MaxT = 3.0
MinT = 0.0
numT = 31
H = np.array([0, 0, 1])
J1 = 1
step = 10000
L = 10
naccept = 0
MT = np.empty((0,2), float)

#init_spin
f = pd.read_csv('random_spin.csv')
kozo = f.values
pos = np.array(kozo[:, 0:3])
spin = np.array(kozo[:,3:6], dtype=float)

#main
for T in np.linspace(MaxT, MinT, numT):
    for loop in range(step):
        for k in range(L):
            for j in range(L):
                for i in range(L):
                    rand_spin = np.random.rand(3) - 0.5
                    rand_spin_norm = np.linalg.norm(rand_spin)
                    next_spin = rand_spin/rand_spin_norm # new spin

                    ID = neighber.get_ID(i, j, k, L)
                    NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5 = neighber.get_NN1(i, j, k, L)
                    E0 = calculation.calc_E0(spin, ID, NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5, J1, H)
                    E1 = calculation.calc_E1(next_spin, spin, NN1_0, NN1_1, NN1_2, NN1_3, NN1_4, NN1_5, J1, H)
                    change = calculation.calc_change(E0, E1, T)

                    metropolis = np.random.uniform(0, 1) #metropolis test
                    if(change > metropolis):
                        spin[ID, 0], spin[ID, 1], spin[ID, 2] = next_spin[0], next_spin[1], next_spin[2]
                        naccept=naccept+1;
    #mag = magnetisation.calc_mag(spin)
    mz = magnetisation.calc_mz(spin)
    MT = np.append(MT, np.array([[T, mz]]), axis=0)
    kozo = 	np.hstack([pos, spin])
    f = pd.DataFrame(kozo, columns=['i', 'j', 'k', 'Sx', 'Sy', 'Sz'])
    f.to_csv('spin' + str(T) + '.csv', index=False)

f = pd.DataFrame(MT, columns = ['T', 'Mz'])
f.to_csv('MT.csv', index=False)
