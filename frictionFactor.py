# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:17:00 2020

@author: cemam
"""
import numpy as np
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