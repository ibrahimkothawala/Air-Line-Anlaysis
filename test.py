#%%
import numpy as np
from scipy.optimize import fsolve
import frictionFactor as ff
import cantera as ct
import matplotlib.pyplot as plt
import DifferentialFannoPressureDrop as DFP
import Orifice_size as ORS



#%% Inputs
pts = 100
mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation
P1_range = np.linspace(1000,2000,pts)
R = 8.314
D_1 = 0.004572 #tube diameter (m) out of pressure regulator is 1/4" flex tube

#%% Calculated Inputs
Ainit = np.pi*(D_1/2)**2 #area of the tube after the regulator



#%% Cantera
gas = ct.Solution('air.cti') 
gas.transport_model = 'Mix'
MW = gas.mean_molecular_weight
Rspec = R/MW*1000 #J/kg-K
gamma = gas.cp/gas.cv


#%% Stagnation Properties
#stagnation properties calcs

for ii in range(pts):
    #gas initial properties
    gas.TPX = T0,P1_range[ii],{'O2':1,'N2':3.76}
    rho = gas.density_mass
    u = mdot/(rho*Ainit)
    dynamicP = 0.5*rho*u**2
    
    #stagnation properties
    P0 =  P1_range[ii] + dynamicP #Pascals stagnation pressure
    a_o = np.sqrt(gamma*Rspec*T0) #stagnation speed of sound, stagnation properties
    rho_o = P0/(Rspec*T0) #stagnation density

    

