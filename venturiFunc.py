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
    return 0#c_d*throatArea*np.sqrt((gamma*rho_0*P0*(2/(gamma+1)))**((gamma+1)/(gamma-1))) needs to be rewritten 

#inputs: mass flowrate,gamma,upstream stagnation density, upstream stagnation pressure, discharge coefficient
#outputs: choked area
#same equation as above did some algebra (note: throat area arbitrarily 1)
def chokedArea(mdot,gamma,rho_0,P0,c_d):
    return 0#mdot/(c_d*np.sqrt((gamma*rho_0*P0*(2/(gamma+1)))**((gamma+1)/(gamma-1)))) needs to be rewritten 

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
# if the throat pressure, upstream pressure, upstream density are known. (not reliable as a design equation)
def venturiMassflow(gamma,rho_ups,pStaticUps,pStatThroat,A_throat,A_ups,c_d):
    return (c_d*A_throat*np.sqrt((2*gamma*rho_ups*pStaticUps)/(gamma-1)*((pStatThroat/pStaticUps)**(2/gamma)-(pStatThroat/pStaticUps)**((gamma+1)/gamma))))     \
    /np.sqrt(1-(A_throat/A_ups)**2*(pStatThroat/pStaticUps)**(2/gamma))

def areaToDiam(area):
    return np.sqrt(4*area/np.pi)

def diamToarea(D):
    return D**2/4*np.pi

if __name__ == '__main__':
    import cantera as ct
    import matplotlib.pyplot as plt

#%% Setup Cantera and input variables
    gas = ct.Solution('air.cti')
    gasMolar = {'O2':1,'N2':3.76} #define gas and molar ratio
    gas.basis = 'mass'
    T0 = 300.0 # stagnation temperature [K]
    pStaticUps = 3e6 # upstream static pressure [Pa] 
    mDot = 0.3 # mass flowrate [kg/s]
    A_ups = diamToarea(12.7e-3) # upstream area [m^2] (based of 1/2 in internal diameter pipe)
    gas.TPX = T0, pStaticUps, gasMolar #setting up gas object
    R_spec = 1000*R/gas.mean_molecular_weight # specific gas constant [j/kg/K] (need to multiply by 100 because cantera use kmole)
    gamma = gas.cp/gas.cv # upstream gamma [unitless] (assume gamma doesn't change across the venturi)

#%% Calculated Inputs
    rho_ups = gas.density_mass # density of the gas upstream of the venturi [kg/m^3] 
    u_Ups = mDot/A_ups/rho_ups # average speed of gas upstream [m/s]
    P0_ups = pStaticUps*stagnationPressureRatio(gamma,mach(u_Ups,speedOfSound(gamma,R_spec,T0))) # upstream stagnation pressure [Pa]
    Pthroat_choked = P0_ups*PchokeCondition(gamma) # static pressure at throat at choking conditions [Pa]

#%% Testing venturi mass flowrate function
    pts = 100
    pStaticUps_range = np.linspace(2e6,10e6,pts) # range of upstream pressures to test
    Pthroat = 2e6 # single downstream pressure [Pa]
    mdot = 0.3 # goal mass flow rate [kg/s]
    soln = np.zeros((pts,3))

    for ii in range(pts):
        gas.TP = T0, pStaticUps_range[ii]
        rho = gas.density_mass
        soln[ii,0] = venturiMassflow(gamma,rho,pStaticUps_range[ii],Pthroat,A_ups/2,A_ups,0.95)
        soln[ii,1] = chokedArea(mdot,gamma,rho,pStaticUps_range[ii],0.95)
        soln[ii,2] = chokedMassFlow(A_ups/2,gamma,rho,pStaticUps_range[ii],0.95)

#%% Plotting
    fig, ax1 = plt.subplots(constrained_layout =True)

    ax1.plot(pStaticUps_range,soln[:,0],label = "mass flowrate against upstream pressure")
    ax1.set_xlabel("upstream pressure [pa]")
    ax1.set_ylabel("mass flowrate [kg/s]")
    ax1.legend()
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("throat diameter [mm]")
    ax2.plot(pStaticUps_range,areaToDiam(soln[:,1])*1e3,label = "area against upstream stagnation pressure")
    ax2.legend()
    plt.legend()
    plt.show()

# %%
