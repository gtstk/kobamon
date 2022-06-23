#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas as pd
import numpy as np

df = pd.read_csv('sample_neuron.dat', skiprows=30, usecols=[2,3,4]) #.dat、上３０行無視、２、３、４列インプット
x = np.array(df)
m = x[:, 2].reshape(-1, 55) #行：温度、　列：磁場
H = x[0:55, 0] #磁場（データ）
T = np.arange(10, 60, 10) #温度（生成）

dS_all = np.empty((0, 4)) #dM/dTbox, columns:磁場
for i in range(H.size):
    if i==0:
        dH = 0
    else:
        dH = (H[i]-H[i-1])
    dS = np.empty(0) #dM/dTbox初期化, columns:磁場
    for j in range(T.size-1):
        #T[j+1]-T[j]
        y = (m[j+1, i]-m[j, i])*dH/10
        dS = np.append(dS, y)
    dS_all = np.append(dS_all, [dS], axis=0)
print(dS_all)
S = np.sum(dS_all, axis=0)
print(S)
T = np.arange(10, 50, 10) #温度（生成）
ST = np.stack([T, S], 1)
f = pd.DataFrame(ST, columns=['T', 'S'])
f.to_csv('ST.csv', index=False)


# In[ ]:




