#!/usr/bin/env python
# coding: utf-8

# In[59]:


import numpy as np

#init parameter
L = 10
size = 10**3
step = 10000
J1 = 1
mu = 1
app = np.array([0, 0, 1])
kB = 1.36*10**(-23)
T = 300
B = 1

#create random spin
rand_row = np.random.rand(3)-0.5
rand_norm = np.linalg.norm(rand_row)
rand_spin = rand_row/rand_norm
print(rand_spin)

#init spin
init_spin = np.empty((0, 3))
str_spin = np.empty((0, 6))
for i in range(L):
    for j in range(L):
        for k in range(L):
            rand_row = np.random.rand(3)-0.5
            rand_norm = np.linalg.norm(rand_row)
            rand_spin = rand_row/rand_norm
            
            init_spin = np.append(init_spin, [rand_spin], axis=0) #スピンのみ
            
            str0_spin = np.append([i, j, k], rand_spin)
            str_spin = np.append(str_spin, [str0_spin], axis=0) #スピンマップ
# got init_spin,str_spin


# In[67]:


# main monte
spin = init_spin
for k in range(L):
    for j in range(L):
        for i in range(L):
            rand_row = np.random.rand(3)-0.5
            rand_norm = np.linalg.norm(rand_row)
            # new spin
            new_spin = rand_row/rand_norm
            
            # origin spin ID
            ID = L*L*k + L*j + i
            # nearest neighbor spin ID
            NN1_0 = L*L*k + L*j + (i-1)%L
            NN1_1 = L*L*k + L*j + (i+1)%L
            NN1_2 = L*L*k + L*((j-1)%L) + i
            NN1_3 = L*L*k + L*((j+1)%L) + i
            NN1_4 = L*L*((k-1)%L) + L*j + i
            NN1_5 = L*L*((k+1)%L) + L*j + i
            
            # culc enegy of origin spin
            sum_NN = spin[NN1_0]+spin[NN1_1]+spin[NN1_2]+spin[NN1_3]+spin[NN1_4]+spin[NN1_5]
            exc0 = -J1*(np.dot(spin[ID], sum_NN))
            app0 = -mu*(np.dot(spin[ID], app))
            E0 = exc0 + app0
            
            # culc enegy of new spin
            exc1 = -J1*(np.dot(new_spin, sum_NN))
            app1 = -mu*(np.dot(new_spin, app))
            E1 = exc1 + app1
            
            # dE, p
            dE = E1 - E0
            p = np.exp(-dE*B)
            if 1<p:
                p=1
                print(E0, E1, dE, p)
                print("true")
            else:
                print(E0, E1, dE, p)
                print("false")


# In[ ]:




