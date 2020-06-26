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
gas.TPX = T0,P1,{'O2':1,'N2':3.76}
rho = gas.density_mass
MW = gas.mean_molecular_weight
#%%Inputs
mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation

