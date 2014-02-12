# -*- coding: utf-8 -*-
"""
This program extracts the risk-neutral densities of inflation expectations using
data on interest rate caps and floors. The columns of the data matric represent strike prices 
and the rows represent different observation dates. We implement the procedure as in 
Breeden & Litzenberger (2013) and plot the market expectations.

@author: Hamed Faquiryan
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

b = pd.io.parsers.read_csv(r'XXXXX\ircaps_2y.csv', index_col=[0], parse_dates=True)

plt.figure();
b.plot();
plt.legend(('1%','2%','3%','4%','5%'), loc='best'); 
plt.title('2-Year IR Cap Prices by Strike'); 
plt.xlabel('Date'); 
plt.ylabel('Price');
plt.savefig(r'C:\Users\l1hsf01\Desktop\ircaps_2y.png', dpi=300)


h = b.values

g = np.empty((np.size(h,0),np.size(h,1)))

#Procedure for RND extraction a la Breeden & Litzenberger (2013)
for jj in range(np.size(h,1)):
    for ii in range(np.size(h,0)):
        if (jj == 0):
            g[ii,jj] = h[ii,jj]
        elif (jj == 3):
            g[ii,jj] = h[ii,jj] - h[ii,jj+1]
        elif (jj >= 4):
            g[ii,jj] = g[ii,jj-1]
        else:
            g[ii,jj] = h[ii,jj-1] + h[ii,jj+1] - 2*h[ii,jj]


for ii in range(np.size(h,0)):
    q = np.sum(g[ii,:])
    for jj in range(np.size(h,1)):
        g[ii,jj] = g[ii,jj]/q


final = pd.DataFrame(g[:,0:4], index=b.index)
            
plt.figure();
final.plot();
plt.legend(('1%','2%','3%','4%'), loc='best'); 
plt.title('RNDs Extracted from 2-Year IR Caps'); 
plt.xlabel('Date'); 
plt.ylabel('Density');
plt.savefig(r'XXXX\ircaps_2y_Density.png', dpi=300)


g[:,0] = 1*g[:,0]  
g[:,1] = 2*g[:,1]  
g[:,2] = 3*g[:,2]  
g[:,3] = 4*g[:,3]  
g[:,4] = 5*g[:,4]

temp = np.sum(g, axis=1)
sd_up = np.std(g, axis=1) + temp
sd_down = (-1)*np.std(g, axis=1) + temp

temper = np.array([temp, sd_up, sd_down])
temper = temper.T

avg = pd.DataFrame(temper, index=b.index)



plt.figure();
avg.plot();
plt.legend(('E[LIBOR]','StD[+]','StD[-]')); 
plt.title('E[LIBOR] in 2 Years based on IR Caps'); 
plt.xlabel('Date'); 
plt.ylabel('Interest Rate');
plt.savefig(r'XXXX\ircaps_2y_mean.png', dpi=300) 

print 'This program has run.'     
    