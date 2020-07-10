#author: Ibrahim Kothawala 

import numpy as np
from scipy.optimize import fsolve
#%% Functions and Global Constants
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
def chokedMassFlow_wiki(throatArea,gamma,rho_0,P0,c_d):
    return c_d*throatArea*np.sqrt(gamma*rho_0*P0*(2/(gamma+1)**((gamma+1)/(gamma-1))))

#inputs stagnation pressure, throat area, stagnation temperature, specific gas constant, gamma
#outputs the maximum mass flowrate through that area
#References: pg 79 Gas Dynamics James John
def chokedMassFlow_gasDyn(P0,throatArea,T0,Rspec,gamma):
    return P0*throatArea/np.sqrt(Rspec*T0)*np.sqrt(gamma)*((gamma+1)/2)**((gamma+1)/(-2*(gamma-1)))

#inputs: mass flowrate,gamma,upstream stagnation density, upstream stagnation pressure, discharge coefficient
#outputs: choked area
#same equation as chokedMassFlow_wiki did some algebra (note: throat area arbitrarily 1)
def chokedArea_wiki(mdot,gamma,rho_0,P0,c_d):
    return mdot/chokedMassFlow_wiki(1,gamma,rho_0,P0,c_d)                                           

#inputs: mass flowrate, stagnation pressure, stagnation temperature, specific gas constant, gamma
#outputs: choked area
#same equation as chokedMassFlow_gasDyn did some algebra (note: throat area arbitrarily 1)
def chokedArea_gasDyn(mdot,P0,T0,Rspec, gamma):
    return mdot/chokedMassFlow_gasDyn(P0,1,T0,Rspec,gamma)

#inputs: gamma,specific gas constant,absolute temperature
#outputs: local speed of sound
#References: https://en.wikipedia.org/wiki/Speed_of_sound
def speedOfSound(gamma,Rspec,T):
    return np.sqrt(gamma*Rspec*T)

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

#inputs: initial mach number guess,area ratio of a star to a actual(aActual/Astar), gamma
#outputs: mach number at the area specified
#references: pg 79 eqn 3.23 and pg 80 for derivative function Gas dynamics james john 
#note: since two solutions for isentropic flow through a variable area channel one supersonic one subsonic if initial
# guess is subsonic will return subsonic root if initial guess is supersonic will return supersonic root.
def machAreaRatio(machGuess,areaRatio,gamma):
    def f(mach):
        return 1/mach*((2/(gamma+1))*(1+mach**2/2*(gamma-1)))**((gamma+1)/2/(gamma-1))-areaRatio
    def fprime(mach):
        b = (gamma+1)/2/(gamma-1); c = 2/(gamma+1)
        return areaRatio - c**(b-1)*mach*(1+mach**2/2/b/c)**(b-1)
    return fsolve(f,machGuess,fprime=fprime)

#inputs: mach number before shock, gamma
#outputs: ratio of stagnation pressures before and after shock p_o2/p_o1
#refernces: pg 119 eqn 4.15 Gas Dynamics james john 
def stagPratioAcrossNormShock(mach1, gamma):
    return (((gamma+1)/2*mach1**2)/(1+(gamma-1)/2*mach1**2))**(gamma/(gamma-1))*(1/(2*gamma*mach1**2/(gamma+1)-(gamma-1)/(gamma+1)))**(1/(gamma-1))

#inputs: area of throat, exit area, shock area guess, gamma, upstream stagnation pressure
#outputs: ratio of upstream stagnation pressure to exit pressure for a shock that occurs at that area (Pb/P01)
#references: pg 137 Gas Dynamics James John
#note this is used in other functions not useful by itself
# currently does not work
def BackPressureAfterNormalShock(throatArea,exitArea,shockArea, gamma, P0):
    mach1 = machAreaRatio(5,shockArea/throatArea,gamma) #calculates the mach number right before shock
    P0ratioShock = stagPratioAcrossNormShock(mach1,gamma) #note this is equivalent to p_02/p_01 = A*_1/A*_2 eqn 4.21 Gas Dynamamics 
    AexitA2starRatio = exitArea*P0ratioShock/exitArea #
    machExit = machAreaRatio(0.4,AexitA2starRatio,gamma)
    PbP02Ratio = 1/stagnationPressureRatio(gamma,machExit) #p_staticExit/p02
    return PbP02Ratio*P0ratioShock*P0

#inputs: area of throat, exit area, gamma, upstream stagnation pressure, exit static pressure
#outputs: if normal shock occurs: location of normal shock
# if normal shock doesn't occur: returns 0 and prints message normal shock doesn't occur within diverging section
# if flow is not choked: message returns 0 and prints flow is not choked
# currently does not work 
def normalShockLoop(throatArea,exitArea,gamma, P0, Pb):
    # chokedCondition = 1/PchokeCondition(gamma)
    # if P0/Pb < chokedCondition: 
    #     print("flow is not choked")
    #     return 0
    # else:
    #     PbP01ratioCalcExit = normalShockCalc(throatArea,exitArea,exitArea,gamma,P0) # calculates the back pressure ratio to upstream stagnation pressure ratio for a shock that occurs at the exit
    #     if PbP01ratioCalcExit < Pb/P0:
    #         print("normal shock does not occur in diverging section")
    #         return 0
    #     else:
    #         shockArea = np.average([throatArea,exitArea])
    #         PbP01ratioCalcLoop = normalShockCalc(throatArea,exitArea,shockArea,gamma,P0)
    #         ctr = 0; ctrMax = 100; tol = 1
    #         while ctr < ctrMax and np.abs(PbP01ratioCalcLoop - Pb/P0) > tol:
    #             shockArea = ((PbP01ratioCalcLoop)/(Pb/P0) -1)*shockArea
    #             PbP01ratioCalcLoop = normalShockCalc(throatArea,exitArea,shockArea,gamma,P0)
    #             ctr+=1
    #         if ctr == ctrMax:
    #             print("convergence failed")
    #         return shockArea
    return 0 

def areaToDiam(area):
    return np.sqrt(4*area/np.pi)

def diamToarea(D):
    return D**2/4*np.pi

def mmToIn(mm):
    return mm/25.4

def inTomm(inches):
    return inches*25.4

def paToPsi(pa):
    return pa/6895

def psiToPa(psi):
    return psi*6895

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
    soln = np.zeros((pts,5))

    for ii in range(pts):
        gas.TP = T0, pStaticUps_range[ii]
        rho = gas.density_mass
        gamma = gas.cp/gas.cv
        soln[ii,0] = venturiMassflow(gamma,rho,pStaticUps_range[ii],Pthroat,A_ups/2,A_ups,1)
        soln[ii,1] = chokedArea_wiki(mdot,gamma,rho,pStaticUps_range[ii],1)
        soln[ii,2] = chokedMassFlow_wiki(A_ups/2,gamma,rho,pStaticUps_range[ii],1)
        soln[ii,3] = chokedMassFlow_gasDyn(pStaticUps_range[ii],A_ups/2,T0,R_spec,gamma)
        soln[ii,4] = chokedArea_gasDyn(mdot,pStaticUps_range[ii],T0,R_spec,gamma)

#%% Testing machAreaRatio function
    ptsMach = 100
    gamma = 1.4
    areaRatioRange = np.linspace(1,10,ptsMach)
    solnMach = np.zeros(ptsMach*2)
    subSonicGuess = 0.4; supSonicGuess = 5
    for ii in range(ptsMach):
            solnMach[ii] = machAreaRatio(subSonicGuess,areaRatioRange[ii],gamma)
            solnMach[ii+ptsMach] = machAreaRatio(supSonicGuess,areaRatioRange[ii],gamma)
     

#%% Testing Normal Shock Functions
    ptsShock = 100
    gamma = 1.4
    exitArea = diamToarea(1e-3*inTomm(0.5))
    solnShock = np.zeros((ptsShock,2))
    P0_range = np.linspace(2e6,10e6,ptsShock)

    for ii in range(ptsShock):
        throatArea = chokedArea_gasDyn(mdot,P0_range[ii],T0,R_spec,gamma)
        solnShock[ii,0] = BackPressureAfterNormalShock(throatArea,exitArea,throatArea,gamma,P0_range[ii])
        solnShock[ii,1] = BackPressureAfterNormalShock(throatArea,exitArea,exitArea,gamma,P0_range[ii])
   
#%% Plotting
    fig, ax1 = plt.subplots(3,constrained_layout =True)
    ax1[0].set_title("Mass flow rate against upstream pressure")
    ax1[0].plot(pStaticUps_range,soln[:,0],label ="Throat Pressure: "+"{:.3f}".format(Pthroat/1e6)+" [MPa]"+" Throat diameter: "+"{:.3f}".format(1e3*areaToDiam(A_ups/2))+" [mm]")
    ax1[0].set_xlabel("upstream static pressure [Pa]")
    ax1[0].set_ylabel("mass flowrate [kg/s]")
    secaxx = ax1[0].secondary_xaxis('top',functions = (paToPsi,psiToPa))
    secaxx.set_xlabel("upstream static pressure [psi]")
    ax1[0].legend()
    
    ax1[1].set_title("Choked diameter against upstream stagnation pressure for "+"{:.2f}".format(mdot)+" [kg/s] mass flowrate")
    ax1[1].set_ylabel("throat diameter [mm]")
    ax1[1].set_xlabel("upstream stagnation pressure [Pa]")
    ax1[1].plot(pStaticUps_range,areaToDiam(soln[:,1])*1e3,label = "(wiki func)")
    ax1[1].plot(pStaticUps_range,areaToDiam(soln[:,4])*1e3,label = "(gas dyn func)")
    secaxy = ax1[1].secondary_yaxis('right',functions = (mmToIn,inTomm))
    secaxy.set_ylabel("throat diameter [in]")
    secaxx1 = ax1[1].secondary_xaxis('top',functions = (paToPsi,psiToPa))
    secaxx1.set_xlabel("upstream stagnation pressure [psi]")
    ax1[1].legend()
    
    ax1[2].set_title("Mass flow rate against upstream pressure")
    ax1[2].set_ylabel("Mass flow rate [kg/s]")
    ax1[2].set_xlabel("upstream stagnation pressure [Pa]")
    ax1[2].plot(pStaticUps_range,soln[:,2],label = "(Wiki func) Throat diameter: "+"{:.3f}".format(1e3*areaToDiam(A_ups/2))+" [mm]")
    ax1[2].plot(pStaticUps_range,soln[:,3],label = "(Gas Dyn func) Throat diameter: "+"{:.3f}".format(1e3*areaToDiam(A_ups/2))+" [mm]")
    ax1[2].legend()
    secaxx2 = ax1[2].secondary_xaxis('top',functions = (paToPsi,psiToPa))
    secaxx2.set_xlabel("upstream stagnation pressure [psi]")
    plt.show()

#%% Plotting Mach number tests
    plt.figure(2)
    plt.title("Area Ratio against Mach Number")
    plt.ylabel("Area Ratio")
    plt.xlabel("Mach Number")
    plt.plot(solnMach[0:ptsMach],areaRatioRange,color = "blue")
    plt.plot(solnMach[ptsMach:2*ptsMach],areaRatioRange,color = "blue")
    plt.show()

#%% Plotting Normal Shock Area tests
    plt.figure(3)
    plt.title("Back pressure required normal shock at throat or exit")
    plt.ylabel("Back Pressure [Pa]")
    plt.xlabel("Upstream stagnation pressure [Pa]")
    plt.plot(P0_range,solnShock[:,0], label = "Back Pressure for Normal Shock at Throat")
    plt.plot(P0_range,solnShock[:,1], label = "Back Pressure for Normal Shock at Exit")
    plt.legend()
    plt.show()
# %%
