#!/usr/bin/env python
# coding: utf-8

# In[217]:


import numpy as np

#init parameter
L = 10
size = 10**3

#create random spin
rand_row = np.random.rand(3)-0.5
rand_norm = np.linalg.norm(rand_row)
rand_spin = rand_row/rand_norm

#init spin
init_spin = np.empty((L, L, L, 3))
str_spin = np.empty((0, 6))
for k in range(L):
    for j in range(L):
        for i in range(L):
            rand_row = np.random.rand(3)-0.5
            rand_norm = np.linalg.norm(rand_row)
            rand_spin = rand_row/rand_norm
            
            init_spin[k, j, i, :] = rand_spin #(L, L, L, 3)スピンのみ
            
            str0_spin = np.append([k, j, i], rand_spin)
            str_spin = np.append(str_spin, [str0_spin], axis=0) #スピンマップ
print(init_spin[0, 0, 0, :])
# got init_spin,str_spin


# In[222]:


# equi
naccept = 0
equi = 1000
them = 10
step = 1000
J1 = 1
mu = 1
H = np.array([0, 0, 1])
kB = 1
T = 0.1
B = 1/(kB*T)
E_total = np.empty(0)
E_exc = np.empty(0)
E_app = np.empty(0)

spin = np.copy(init_spin)

# equilizarion
for loop in range(equi):
    for k in range(L):
        for j in range(L):
            for i in range(L):
                rand_row = np.random.rand(3)-0.5
                rand_norm = np.linalg.norm(rand_row)
                # new spin
                new_spin = rand_row/rand_norm

                # culc enegy of origin spin
                sum_NN = spin[k, j, (i-1)%L, :]+spin[k, j, (i+1)%L, :]+spin[k, (j-1)%L, i, :]+spin[k, (j+1)%L, i, :]+spin[(k-1)%L, j, i, :]+spin[(k+1)%L, j, i, :]
                exc0 = -J1*(np.dot(spin[k, j, i, :], sum_NN))
                app0 = -mu*(np.dot(spin[k, j, i, :], H))
                E0 = exc0 + app0

                # culc enegy of new spin
                exc1 = -J1*(np.dot(new_spin, sum_NN))
                app1 = -mu*(np.dot(new_spin, H))
                E1 = exc1 + app1

                # dE, p
                dE = E1 - E0
                p = np.exp(-dE*B)

                #metropolis test
                metropolis = np.random.uniform(0, 1)
                if 1<p: #確率が1以上の場合承認
                    spin[k, j, i, :] = new_spin
                    naccept += 1;
                elif metropolis<p: #メトロポリス値より確率が高い場合承認
                    spin[k, j, i, :] = new_spin
                    naccept += 1;
                    
    
    # calc total energy
    calc_total, calc_exc, calc_app = 0, 0, 0
    for k in range(L):
        for j in range(L):
            for i in range(L):
                sum_NN = spin[k, j, (i-1)%L, :]+spin[k, j, (i+1)%L, :]+spin[k, (j-1)%L, i, :]+spin[k, (j+1)%L, i, :]+spin[(k-1)%L, j, i, :]+spin[(k+1)%L, j, i, :]
                exc0 = -J1*(np.dot(spin[k, j, i, :], sum_NN))
                app0 = -mu*(np.dot(spin[k, j, i, :], H))
                E0 = exc0 + app0
                calc_exc += exc0
                calc_app += app0
                calc_total += E0
    E_total = np.append(E_total, calc_total)
    E_exc = np.append(E_exc, calc_exc)
    E_app = np.append(E_app, calc_app)
    # エネルギーここまで
# equiここまで


# In[218]:


# MT曲線ここから
equi = 100
them = 100
MaxT = 3.0
MinT = 0.1
numT = 30
kB = 1
T = 0
MT = np.empty((0,2))

for T in np.linspace(MaxT, MinT, numT):
    print("温度変更", T)
    spin = np.copy(init_spin)
    B = 1/(kB*T)
    #　平衡化ここから
    for loop in range(equi):
        for k in range(L):
            for j in range(L):
                for i in range(L):
                    rand_row = np.random.rand(3)-0.5
                    rand_norm = np.linalg.norm(rand_row)
                    # new spin
                    new_spin = rand_row/rand_norm

                    # culc enegy of origin spin
                    sum_NN = spin[k, j, (i-1)%L, :]+spin[k, j, (i+1)%L, :]+spin[k, (j-1)%L, i, :]+spin[k, (j+1)%L, i, :]+spin[(k-1)%L, j, i, :]+spin[(k+1)%L, j, i, :]
                    exc0 = -J1*(np.dot(spin[k, j, i, :], sum_NN))
                    app0 = -mu*(np.dot(spin[k, j, i, :], H))
                    E0 = exc0 + app0

                    # culc enegy of new spin
                    exc1 = -J1*(np.dot(new_spin, sum_NN))
                    app1 = -mu*(np.dot(new_spin, H))
                    E1 = exc1 + app1

                    # dE, p
                    dE = E1 - E0
                    p = np.exp(-dE*B)

                    #metropolis test
                    metropolis = np.random.uniform(0, 1)
                    if 1<p: #確率が1以上の場合承認
                        spin[k, j, i, :] = new_spin
                        naccept += 1;
                    elif metropolis<p: #メトロポリス値より確率が高い場合承認
                        spin[k, j, i, :] = new_spin
                        naccept += 1;
    # 平衡化ここまで, spin(4D)
    
    # 熱平均ここから
    M_all = np.empty(0)
    
    for loop in range(them):
        for k in range(L):
            for j in range(L):
                for i in range(L):
                    rand_row = np.random.rand(3)-0.5
                    rand_norm = np.linalg.norm(rand_row)
                    # new spin
                    new_spin = rand_row/rand_norm

                    # culc enegy of origin spin
                    sum_NN = spin[k, j, (i-1)%L, :]+spin[k, j, (i+1)%L, :]+spin[k, (j-1)%L, i, :]+spin[k, (j+1)%L, i, :]+spin[(k-1)%L, j, i, :]+spin[(k+1)%L, j, i, :]
                    exc0 = -J1*(np.dot(spin[k, j, i, :], sum_NN))
                    app0 = -mu*(np.dot(spin[k, j, i, :], H))
                    E0 = exc0 + app0

                    # culc enegy of new spin
                    exc1 = -J1*(np.dot(new_spin, sum_NN))
                    app1 = -mu*(np.dot(new_spin, H))
                    E1 = exc1 + app1

                    # dE, p
                    dE = E1 - E0
                    p = np.exp(-dE*B)

                    #metropolis test
                    metropolis = np.random.uniform(0, 1)
                    if 1<p: #確率が1以上の場合承認
                        spin[k, j, i, :] = new_spin
                        naccept += 1;
                    elif metropolis<p: #メトロポリス値より確率が高い場合承認
                        spin[k, j, i, :] = new_spin
        calc_m = np.empty(0)
        calc_mz = np.empty(0)
        calc_m = np.mean(spin, axis=0) #3D
        calc_m = np.mean(calc_m, axis=0) #2D
        calc_m = np.mean(calc_m, axis=0) #1D
        calc_mz = calc_m[2] 
        M_all = np.append(M_all, calc_mz)
    M = np.mean(M_all)
    MT = np.append(MT, np.array([[T, M]]), axis=0)
print(MT)
        # 熱平均ここまで


# In[221]:


#データ加工、アウトプット
import pandas as pd
E_S = np.empty((0, 4), float)
S = np.arange(1, equi+1)
E_S = np.stack([S, E_total, E_exc, E_app], 1)

f = pd.DataFrame(E_S, columns=['step', 'E_total', 'E_exc', 'E_app'])
f.to_csv('ES.csv', index=False)

f = pd.DataFrame(MT, columns=['Temp', 'z-Mag'])
f.to_csv('MTcurve.csv', index=False)


# In[206]:


for T in np.linspace(3, 0.1, 30):
    B = 1/(T)
    p = np.exp(-dE*B)
    print(T)
    print(B)
    print(p)


# In[ ]:




