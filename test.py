#%%
import numpy as np
from scipy.optimize import fsolve
import frictionFactor as ff
import cantera as ct
import matplotlib.pyplot as plt
import DifferentialFannoPressureDrop as DFP


#&&Cantera
gas = ct.Solution('air.cti') 
gas.transport_model = 'Mix'
#%%Inputs
mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation
Tstag = T0
P1 =10.4321e6#MPa #pressure at pipe inlet [Pa] 1500 PSI 
PInit = P1
Rspec = 287#J/kg-K
gas.TPX = T0,P1,{'O2':1,'N2':3.76}
rho = gas.density_mass
rhoInit = rho
gamma = 1.4 #ratio of specific heats [-]
#visc 
visc = gas.viscosity
#step size for analysis
step = 0.001 #(m) break tube up into 1mm increments
#for loss calculations https://docs.google.com/spreadsheets/d/1gn8gLvc1uLY0blPmABOoTxjrHLD8UaT0t1OvKjpqaj4/edit#gid=0
d = 0.004572
l = 3
eps = 0.00001
sumk = 0.750

DFP.Fanno(mdot,T0,P1,gas, gamma, d, l,step,eps,sumk)