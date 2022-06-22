#!/usr/bin/env python
# coding: utf-8

# In[172]:


import numpy as np

#init parameter
L = 10
size = 10**3
step = 1000

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


# In[174]:


# main monte
naccept = 0
step = 1000
J1 = 1
mu = 1
H = np.array([0, 0, 1])
kB = 1.36*10**(-23)
T = 300
B = 1
E_total = np.empty(0)
E_exc = np.empty(0)
E_app = np.empty(0)

#確認
print(init_spin[0, 0, 0, :], "init_spin ここまで")
spin = np.copy(init_spin)
print(spin[0, 0, 0, :], "spin ここまで")

# main
for loop in range(step):
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
print(spin[0, 0, 0, :], "spin ここまで")
print("log", "accept_ratio:", naccept/(size*step), E_total, E_exc, E_app)


# In[180]:


#データ加工、アウトプット
import pandas as pd
E_S = np.empty((0, 4), float)
S = np.arange(1, 1001)
E_S = np.stack([S, E_total, E_exc, E_app], 1)

f = pd.DataFrame(E_S, columns=['step', 'E_total', 'E_exc', 'E_app'])
f.to_csv('ES.csv', index=False)


# In[ ]:




