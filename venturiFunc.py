#author: Ibrahim Kothawala 

import numpy as np

#Gas Constant
R = 8.31446261815324 #j/K/mol

#inputs: gamma
#outputs: Pstar/P0
#References: https://en.wikipedia.org/wiki/Choked_flow
#pg 601 eqn 15.46 Gas Dynamics James John 
def PchokeCondition(gamma): #pressure ratio of Pthroat/P_UpstreamStagnation to velocity choke flow
    return (2/(gamma + 1))**(gamma/(gamma-1))

#inputs: throat area, gamma, stagnation density, stagnation pressure, discharge coefficient
#outputs: choked mass flow rate through throat area.
#References: https://en.wikipedia.org/wiki/Choked_flow
def chokedMassFlow(throatArea,gamma,rho_0,P0,c_d):
    return c_d*throatArea*np.sqrt((gamma*rho_0*P0*(2/(gamma+1)))**((gamma+1)/(gamma-1)))

#inputs: gamma,specific gas constant,absolute temperature
#outputs: local speed of sound
#References: https://en.wikipedia.org/wiki/Speed_of_sound
def speedOfSound(gamma,R_spec,T):
    return np.sqrt(gamma*R_spec*T)

#inputs: velocity, local speed of sound
#outputs: local mach number
#refernces: https://en.wikipedia.org/wiki/Mach_number
def mach(u,c):
    return u/c

#inputs: gamma mach number
#outputs: P0/pStatic
#reference: pg 75 eqn 3.15 Gas Dynamics James John
def stagnationPressureRatio(gamma,mach):
    return (1+(gamma-1)/2*mach*mach)**(gamma/(gamma-1))

#inputs: stagnationPressureRatio, gamma
#outputs: rhoStagnation/rhoStatic
#reference: pg 76 right above eqn 3.17 Gas Dynamics James John
def stagnationDensityRatio(stagnationPressureRatio,gamma):
    return stagnationPressureRatio**(1/gamma)

#inputs: gamma, upstream density, upstream static pressure, throat static pressure
# throat area, upstream area, discharge coefficient
#outputs: the mass flowrate through a venturi meter
#references: pg 598 eqn 15.33 Gas Dynamics James John
#additional information: can be used to determine the mass flowrate thru a venturi
# if the throat pressure, upstream pressure, upstream density are known.
def venturiMassflow(gamma,rho_ups,pStaticUps,pStatThroat,A_throat,A_ups,c_d):
    return (c_d*A_throat*np.sqrt((2*gamma*rho_ups*pStaticUps)/(gamma-1)*((pStatThroat/pStaticUps)**(2/gamma)-(pStatThroat/pStaticUps)**((gamma+1)/gamma))))     \
    /np.sqrt(1-(A_throat/A_ups)**2*(pStatThroat/pStaticUps)**(2/gamma))

if __name__ == '__main__':
    import cantera as ct
    
#%% Setup Cantera and input variables
    gas = ct.Solution('air.cti')
    gasMolar = {'O2':1,'N2':3.76} #define gas and molar ratio
    gas.basis = 'mass'
    T0 = 300.0 # stagnation temperature [K]
    pStaticUps = 3e6 # upstream static pressure [Pa] 
    mDot = 0.3 # mass flowrate [kg/s]
    A_ups = np.pi/4*(12.7e-3)**2 # upstream area [m^2] (based of 1/2 in internal diameter pipe)
    gas.TPX = T0, pStaticUps, gasMolar #setting up gas object
    R_spec = 1000*R/gas.mean_molecular_weight # specific gas constant [j/kg/K] (need to multiply by 100 because cantera use kmole)
    gamma = gas.cp/gas.cv # upstream gamma [unitless] (assume gamma doesn't change across the venturi)

#%% Calculated Inputs
    rho_ups = gas.density_mass # density of the gas upstream of the venturi [kg/m^3] 
    u_Ups = mDot/A_ups/rho_ups # average speed of gas upstream [m/s]
    P0_ups = pStaticUps*stagnationPressureRatio(gamma,mach(u_Ups,speedOfSound(gamma,R_spec,T0))) # upstream stagnation pressure [Pa]
    Pthroat = P0_ups*PchokeCondition(gamma) # static pressure at throat at choking conditions [Pa]

#%% Testing venturi mass flowrate function
    pts = 100
    pStaticUps_range = np.linspace(2e6,5e6,pts)
    