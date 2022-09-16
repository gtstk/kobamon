import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from numpy.random import *

#para
L = 5

#function
def set_new_spin():
    theta, fai = np.radians(randint(0, 181)), np.radians(randint(0, 361))
    x, y, z = np.sin(theta)*np.cos(fai), np.sin(theta)*np.sin(fai), np.cos(theta)
    new_spin = np.array([x, y, z])
    return new_spin

#init spin
init_spin = np.empty((L, L, L, 3))
for k in range(L):
    for j in range(L):
        for i in range(L):
            new_spin = set_new_spin()
            init_spin[k, j, i, :] = new_spin
print("L:", L)
print(init_spin)

#check equiliblium

#para
J1 = 1
H = 0.1*np.array([0, 0, 1])
kB = 1
T = 1.5
B = 1/(kB*T)
equi = 1000

#load spin
spin = np.copy(init_spin)

#function
def calc_sum_nn(spin, i, j, k, L):
    sum_nn1 = spin[k, j, (i+1)%L, :] + spin[k, j, (i-1)%L, :] + spin[k, (j+1)%L, i, :] + spin[k, (j-1)%L, i, :] + spin[(k+1)%L, j, i, :] + spin[(k-1)%L, j, i, :]
    sum_nn = sum_nn1
    return sum_nn


#main
size = L**3
naccept = 0
plt_E = np.empty(0)
plt_mcs = np.empty(0)
for loop in range(equi):
    plt_mcs = np.append(plt_mcs, loop)
    for mcs in range(size):
        new_spin = set_new_spin()
        i, j, k = randint(0, L, 3)
        sum_nn = calc_sum_nn(spin, i, j, k, L)
        E0 = -np.dot((J1*sum_nn+H), spin[k, j, i, :])
        E1 = -np.dot((J1*sum_nn+H), new_spin)
        dE = E1 - E0
        metro_p = np.exp(-B*dE)
        metropolis = uniform(0,1)
        #print("spin:", spin[k, j, i, :])
        #print("new_spin:", new_spin)
        #print("metro_p:", metro_p,"\n")
        if(metro_p>metropolis):
            spin[k, j, i, :] = new_spin
            naccept += 1
    #calc enegy
    E_total = 0
    for k in range(L):
        for j in range(L):
            for i in range(L):
                sum_nn = spin[k, j, (i+1)%L, :] + spin[k, j, (i-1)%L, :] + spin[k, (j+1)%L, i, :] + spin[k, (j-1)%L, i, :] + spin[(k+1)%L, j, i, :] + spin[(k-1)%L, j, i, :]
                E0 = -np.dot(J1*sum_nn, spin[k, j, i, :])
                E_total += E0
    plt_E = np.append(plt_E, E_total)
print("fin equi, accept ratio:", naccept/(equi*size))
#plot equi
plt.title("Equilibrium")
plt.xlabel("MCS")
plt.ylabel("Energy/site")
plt.plot(plt_mcs[:], plt_E[:])
plt.show()

#plot spin
#4D~2D
sm = np.empty((size, 6))
for k in range(L):
    for j in range(L):
        for i in range(L):
            id = L**2*k+L*j+i
            sm[id, :] = [i, j, k, spin[k, j, i, 0], spin[k,j,i,1], spin[k,j,i,2]]
# FigureとAxes
fig = plt.figure(figsize = (10, 10))
ax = fig.add_subplot(111, projection='3d')
#ax.grid()
ax.set_xlabel("x", fontsize = 16)
ax.set_ylabel("y", fontsize = 16)
ax.set_zlabel("z", fontsize = 16)
ax.set_xlim(0, L-1)
ax.set_ylim(0, L-1)
ax.set_zlim(0, L-1)

ax.quiver(sm[:, 0], sm[:, 1], sm[:, 2], sm[:,3], sm[:,4], sm[:,5],
          color = "red", length = 0.5,
          arrow_length_ratio = 0.5)