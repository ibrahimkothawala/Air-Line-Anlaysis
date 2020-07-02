# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:17:00 2020

@author: cemam
"""
import numpy as np
from scipy.optimize import fsolve
#%%Functions
#dh hydraulic diameter (m)
#e roughness (m)
#Re reynolds number
def f(dh,e,Re): #turbulent
    if Re > 2300:
        f = (-1.8*np.log10(((e/dh)/3.7)**1.11+(6.9/Re)))**(-2)
    else: #laminar
        f = 64/Re
    return f

def f_colebrook(dh,e,Re):
    if Re > 2300:
        def func(f_init):
            return 1/np.sqrt(f_init) + 2*np.log(e/3.7/dh +2.51/Re/np.sqrt(f_init))
        f_val = fsolve(func,f(dh,e,Re))
    else: #laminar
        f_val = 64/Re
    return f_val
    
if __name__ == "__main__":

    import matplotlib.pyplot as plt

    e_range = np.logspace(-2,-6,7)
    Re_range = np.logspace(1,8,100)
    d = 1
    funcs = [f,f_colebrook]
    soln = np.zeros((len(funcs)+1,len(e_range),len(Re_range)))
    for ii in range(len(funcs)):
        for jj in range(len(e_range)):
            for kk in range(len(Re_range)):
                soln[ii,jj,kk] = funcs[ii](d,e_range[jj],Re_range[kk])
#%%  
    for ii in range(len(funcs)):
        for jj in range(len(e_range)):
            plt.loglog(Re_range,soln[ii,jj,:],'x',color = 'tab:olive')
   
    for jj in range(len(e_range)):
        for kk in range(len(Re_range)):
                soln[len(funcs),jj,kk] = np.average(soln[:,jj,kk])
    for jj in range(len(e_range)):
            plt.loglog(Re_range,soln[0,jj,:],'-', color ='black')
    plt.title('Moody Diagram Generated with different friction factor calculators')
    plt.xlabel('Reynolds Number')
    plt.ylabel('Friction Factor')
    plt.savefig('friction factor.png')


    



# %%
