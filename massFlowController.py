# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:20:32 2020

@author: Ibrahim Kothawala and Christopher Emami
"""


import numpy as np
import LineLossAnalysis as pDropCalc
import cantera as ct
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import venturiFunc as VF

#&&Cantera
gas = ct.Solution('air.cti') 
gasMolar = {'O2':1,'N2':3.76} #define gas and molar ratio
#Inputs
dPurchase = 0.25 #[inches] the size of the orifice based off of Mccmaster availability or other suppliers
C_d = 0.7 #drilled orifice Cd
D_1 = 0.004572 #Initial upstream tube diameter (m) (right after the regulator)

mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation
P2Goal = 300 #[PSI] #goal Pressure (chamber pressure)
step = 0.001 #(m) #step size for line loss analysis
npts = 100 #number of pressure points for analysis
pressureRange = np.linspace(1200,2500,npts)
soln = np.zeros((7,npts))

Kvec = [0.75, 0.00001006424457,0.2,1.816115216,7.5, 0.5, 3.30000078] #valve and fitting coefficients
Lvec = [1e-16, 1e-16, 1.2192, 0.5461, 1.5367, 1.2192, 0.0635] #lengths of pipes in between valves and fittings
Dvec = [0.004572, 0.009398, 0.008636, 0.009398, 0.009398, 0.008636,0.009398] #diameters of pipes in between valves and fittings
epsvec = [0.015e-3, 0.015e-3, 3e-6, 0.015e-3, 0.001e-3, 3e-6, 0.015e-3]  #mm #https://www.engineeringtoolbox.com/lined-pipe-pressure-loss-d_1178.html ptfe roughness for flex hose
pipeSysProp = [Kvec,Lvec,Dvec,epsvec]

#orifice sizing functions

#equation 15.44 inputs: stagnation density, stagnation speed of sound, stagnation pressure, pressure after orifice. 
#returns the value of the fractional section in eqn 15.44 including the stagnation density
def eqn1544Helper(gamma,rho_o,a_o,P0,P4):
    return ((np.sqrt(((2/(gamma - 1)*(P4/P0)**(2/gamma)*((P0/P4)**((gamma - 1)/gamma) - 1))))/(1+(2/np.pi)*(P4/P0)**(1/gamma)))*rho_o*a_o)

#equation 15.44 inputs: stagnation density, stagnation speed of sound, stagnation pressure, pressure after orifice, orifice throat area.
#returns the mass flow rate through an orifice of the given dimensions at the given pressures. 
def massFlowRateThruOrifice(gamma,rho_o,a_o,P0,P4,orificeThroatArea):
    return eqn1544Helper(gamma,rho_o,a_o,P0,P4)*orificeThroatArea

#equation 15.44 inputs: stagnation density, stagnation speed of sound,stagnation pressure, pressure after orifice, mass flow rate thru orifice. 
#returns the area of the orifice that outputs the given mass flowrate at the given inputs.
def orificeThroatArea(gamma,rho_o,a_o,P0,P4,mdot):
    return mdot/eqn1544Helper(gamma,rho_o,a_o,P0,P4)

#equation 15.44 inputs: stagnation density, stagnation speed of sound,stagnation pressure, mass flow rate thru orifice, orifice throat area, and a guess for the pressure after the orifice. 
#returns the pressure after the orifice for the given mass flow rate and the inputs above.
def pressureAfterOrifice(gamma,rho_o,a_o,P0,mdot,orificeThroatArea,P4_guess):
    return fsolve(lambda Pout: - mdot + eqn1544Helper(gamma,rho_o,a_o,P0,Pout)*orificeThroatArea,P4_guess) #eqn 15.44

#equation 15.47 inputs: stagnation density, stagnation speed of sound, orifice throat area
#returns the maximum massflowrate possible thru an orifice of that area. 
def maxMassFlowRate(gamma,rho_o,a_o,orificeThroatArea):
    return (((2/(gamma + 1))**((gamma + 1)/(2*(gamma -1)))))/(1 + (2/np.pi)*(gamma/(gamma + 1))**(1/(gamma - 1))) * rho_o* a_o*orificeThroatArea

def PchokeCondition(gamma): #pressure ratio of Pthroat/Pdownstream to velocity choke flow
    return (2/(gamma + 1))**(gamma/(gamma-1))

# def venturiConditions(mdot, T0,)

#Evaluate different upstream pressures from 450 to 1500 PSI then graph results
P2Goal = pDropCalc.psi_to_MPa(P2Goal)*1e6 #[Pa]
dPurchase = dPurchase*0.0254 #purchase orifice size from mcmaster
P02vec = np.zeros(npts)


orifice = False
venturi = True

if venturi: 
    solnVenturi = np.zeros((npts,5))
    exitArea = VF.diamToarea(1e-3*VF.inTomm(0.5))

for ii in range(npts):

    InitialPressure = pressureRange[ii] #Insert pressure at pressure regulator here in PSI
    InitialPressure = pDropCalc.psi_to_MPa(InitialPressure)*1e6 #Pascals
    
    #define gas
    gas.transport_model = 'Mix'
    gas.TPX = T0,InitialPressure, gasMolar
    R = 8.314
    gamma = gas.cp/gas.cv
    MW = gas.mean_molecular_weight
    gas.basis = 'mass'
    rhoInit = gas.density #initial density at the beginning of the system
    Rspec = R/MW*1000 #J/kg-K
    #for loss calculations https://docs.google.com/spreadsheets/d/1gn8gLvc1uLY0blPmABOoTxjrHLD8UaT0t1OvKjpqaj4/edit#gid=0
    
    #try catch statement for when the mach number exceeds 1
    try:
        LineLoss = pDropCalc.LineAnalysis(mdot, T0, InitialPressure, gas, gamma, step, pipeSysProp,plots=False)
    except:
        print("Upstream pressure data point thrown out as flow Fanno choked at this mass flowrate: "+str(mdot)+" kg/s and this pressure "+str('%.2f'%pressureRange[ii])+" PSI")
        continue     
        
    P_static = (LineLoss[0][1,-1])# pressure at exit of piping system
    PdropLines = LineLoss[2]#Pressure drop across the lines 
    mach = LineLoss[1][4,-1]
    P0 = P_static*VF.stagnationPressureRatio(gamma,mach)
    P02vec[ii] = P0

    
    if venturi:
        solnVenturi[ii,0] = LineLoss[0][1,0] #static upstream pressure
        solnVenturi[ii,1] = P0 #stagnation downstream pressure
        solnVenturi[ii,2] = VF.chokedArea_gasDyn(mdot,P0,T0,Rspec,gamma) #choked area that that gives specified mass flowrate 
        solnVenturi[ii,3] = VF.BackPressureAfterNormalShock(solnVenturi[ii,2],exitArea,exitArea,gamma,P0) #back pressure needed to produce shock at exit area
        solnVenturi[ii,4] = VF.BackPressureAfterNormalShock(solnVenturi[ii,2],exitArea,solnVenturi[ii,2],gamma,P0) # back pressure needed to produce shock at throat
        

    if orifice:

    
        #minimum downstream pressure to choke flow
        MinimumPressure = P_static/(1+PchokeCondition(gamma)) # 
    
        #stagnation properties calcs
        Ainit = np.pi*(D_1/2)**2
        u = mdot/(rhoInit*Ainit)
        dynamicP = 0.5*rhoInit*u**2
        #P0 = InitialPressure + dynamicP #Pascals stagnation pressure
        a_o = np.sqrt(gamma*Rspec*T0) #stagnation speed of sound, stagnation properties
        rho_o = P0/(Rspec*T0) #stagnation density

        #evaluate pressure drop across orifice, add a factor of safety to not experience pressure fluctuations and
        #to achieve 300 psi chamber pressure, pg. 466 in fluid mechanics textbook
        #%%
        #mdot check
        #pg 602 gas dynamics James E. John
        #given flow rate, what diameter is required, and outlet pressure associated:
        orificeAmin = orificeThroatArea(gamma,rho_o,a_o,P0,P2Goal,mdot) #15.44
        d_orificeMin = np.sqrt((orificeAmin*4)/np.pi) 
        d_inchesMin = d_orificeMin/0.0254
        P4_real = pressureAfterOrifice(gamma,rho_o, a_o, P0, mdot, orificeAmin, P2Goal) #15.44
    
        #given diameter, what mdot is achieved and downstream pressure:
        #observe any corrections due to standard orifice sizing
        
        Apurchase = np.pi*(dPurchase/2)**2 #m^2 
        #theoretical maximum from purchased orifice
        mdotPurchase = maxMassFlowRate(gamma,rho_o,a_o,Apurchase) # eqn 15.47
        #solve for downstream pressure after orifice, want to equal goal P
        P4 = pressureAfterOrifice(gamma, rho_o, a_o, P0, mdot, Apurchase,P2Goal) #eqn 15.44
        #if the following downstream pressure is observed in experimentation, this is the flow rate 
        mdotExperimental = massFlowRateThruOrifice(gamma,rho_o,a_o,P0,P4,Apurchase) #eqn 15.44

        PdropOrifice = (P_static - P4)*1e-6
        if PdropOrifice < 0.0:
            soln[:,ii] = np.zeros((len(soln[:,ii])))
            print("Upstream pressure data point thrown of "+str('%.2f'%pressureRange[ii])+" PSI out as pressure drop across orifice is too high, need more upstream pressure.")
        else:
            PdropSystem = (PdropLines + PdropOrifice)
            print('\nOrifice Results:')
            print('Orifice diameter:', dPurchase/0.0254,'(in)')
            print('Upstream Pressure:', InitialPressure*1e-6,'(MPa)')
            print('Upstream Pressure:', InitialPressure*1e-6*145.038,'(PSI)')
            print('Pressure before orifice:', P_static*1e-6,'(MPa)')
            print('Pressure entering chamber',P4*1e-6,'(MPa)')
            print('Pressure entering chamber',P4*1e-6*145.038,'(PSI)')
            print('Minimum downstream pressure to choke', MinimumPressure*1e-6,'(MPa)')
            print('Pressure drop across lines', PdropLines,'(MPa)')
            print('Pressure drop across orifice ',PdropOrifice,'(MPa)')
            print('Total Pressure drop across system', PdropSystem,'(MPa)')
            print('mdot from experimentation from purchased orifice', mdotExperimental, '(kg/s)')
                
            soln[:,ii] =[P4, PdropLines, PdropOrifice, PdropSystem, d_inchesMin, mdotExperimental, pressureRange[ii]] 

            #%%
        xlabel = 'Upstream Pressure (MPa)'
        #hides all zeros in soln caused by non real results.
        soln2 = np.ma.masked_equal(soln,0)
#%% Plotting Venturi
if venturi:
    estimatedBackPressure = 1.2*2e6
    solnVenturi = np.ma.masked_equal(solnVenturi,0)
    plt.figure(1)
    plt.title("Static Upstream Pressure Against Downstream Stagnation Pressure")
    plt.ylabel("Downstream Pressure")
    plt.xlabel("Upstream Static Pressure")
    plt.plot(solnVenturi[:,0],solnVenturi[:,1],label = "Stagnation Pressure at Pipe System Outlet")
    plt.plot(solnVenturi[:,0],solnVenturi[:,3],label = "Downstream back pressure for normal shock at venturi exit")
    plt.plot(solnVenturi[:,0],solnVenturi[:,4],label = "Downstream back pressure for normal shock at venturi throat")
    plt.plot(solnVenturi[:,0],estimatedBackPressure*np.ones(len(solnVenturi[:,0])),'-.',label = "Estimated back pressure")
    plt.legend()
    
#%% Plotting Orifice






if orifice:
        
    #plotting
    fig, axs = plt.subplots(3, 1, constrained_layout=True)
    axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln2[0,:]/1e6,'-',label = 'Pressure after Orifice')
    axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln2[1,:],'-', label = 'Delta P Across Lines')
    axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln2[2,:],'-', label = 'Delta P Across Orifice')
    axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln2[3,:],'-', label = 'Total Delta P')
    axs[0].legend()
    axs[0].set_title('Changes in System Pressures with Upstream Pressure')
    axs[0].set_xlabel(xlabel)
    axs[0].set_ylabel('Pressure (MPa)')
    #fig.suptitle(' Flow Results', fontsize=12)
    axs[1].set_ylabel('Orifice Diameter [in]')

    axs[1].set_ylabel('Mass Flowrate [kg/s]')
    axs[1].plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[5,:], '-',label = 'Mdot Experimental [kg/s]', color = 'tab:blue')
    axs[1].set_xlabel(xlabel)
    axs[1].set_title('Orifice Properties Against Upstream Pressure')
    axs[1].legend(loc='upper left')

    axs[2].plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[4,:],label = "theoretical orifice area")
    axs[2].set_title("Theoretical Orifice size against upstream pressure \n")
    axs[2].set_ylabel("orifice area [in^2]")
    axs[2].set_xlabel(xlabel)


