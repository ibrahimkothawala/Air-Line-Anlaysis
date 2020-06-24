# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:20:32 2020

@author: Ibrahim Kothawala and Christopher Emami
"""
#testing github stuff
#testing again for fun
#https://guides.github.com/introduction/flow/

import numpy as np
import LineLossAnalysis as pDropCalc
import cantera as ct
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
#Orifice Plate try 300
#&&Cantera
#Inputs
Kvec = [0.75, 0.00001006424457,0.2,1.816115216,7.5, 0.5, 3.30000078]
Lvec = [1e-16, 1e-16, 1.2192, 0.5461, 1.5367, 1.2192, 0.0635]
Dvec = [0.004572, 0.009398, 0.008636, 0.009398, 0.009398, 0.008636,0.009398]
epsvec = [0.015e-3, 0.015e-3, 3e-6, 0.015e-3, 0.001e-3, 3e-6, 0.015e-3]  #mm #https://www.engineeringtoolbox.com/lined-pipe-pressure-loss-d_1178.html ptfe roughness for flex hose
dPurchase = 0.18 #[inches] the size of the orifice based off of Mccmaster availabality or other suppliers
mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation
npts = 5
Tstag = T0
gas = ct.Solution('air.cti') 
pressureRange = np.linspace(1200,1515,npts)
soln = np.zeros((7,npts))
pipeSysProp = [Kvec,Lvec,Dvec,epsvec]
#Evaluate different upstream pressures from 450 to 1500 PSI then graph results
for ii in range(npts):

    InitialPressure = pressureRange[ii] #Insert pressure at pressure regulator here in PSI
    InitialPressure = pDropCalc.psi_to_MPa(InitialPressure)*1e6
    #define gas
    
    gas.transport_model = 'Mix'
    gas.TPX = T0,InitialPressure,{'O2':1,'N2':3.76}
    R = 8.314
    gamma = gas.cp/gas.cv
    MW = gas.mean_molecular_weight
    gas.basis = 'mass'
    rhoInit = gas.density #initial density at the beginning of the system
    #step size for analysis
    step = 0.001 #(m) break tube up into 1mm increments
    #for loss calculations https://docs.google.com/spreadsheets/d/1gn8gLvc1uLY0blPmABOoTxjrHLD8UaT0t1OvKjpqaj4/edit#gid=0
    try:
        LineLoss = pDropCalc.LineAnalysis(mdot, T0, InitialPressure, gas, gamma, step, pipeSysProp,plots=False)
    except:
        print("Upstream pressure data point thrown out as flow Fanno choked at this mass flowrate: "+str(mdot)+"kg/s and this pressure "+str(pressureRange[ii])+"psi")
        continue     
        
    

    #print(LineLoss.)
    T_static = LineLoss[0][2,0] #initial temperature (this is the temperature after one step)
    P_static = (LineLoss[0][1,-1])# pressure upstream of the orifice 
    PdropLines = LineLoss[2]#Pressure drop across the lines 



    #%%
    Rspec = R/MW*1000 #J/kg-K
    C_d = 0.99 #machined venturi #https://neutrium.net/fluid_flow/discharge-coefficient-for-nozzles-and-orifices/
    #goal Downstream pressure
    P2Goal = P_static/1.528 # the 1.528 is the pressure ratio required for choked flow (find ref) #MPa

    #stagnation properties calcs
    D_1 = 0.004572 #tube diameter (m) out of pressure regulator is 1/4" flex tube
    Ainit = np.pi*(D_1/2)**2
    u = mdot/(rhoInit*Ainit)
    dynamicP = 0.5*rhoInit*u**2
    P0 = InitialPressure + dynamicP #Pascals
    a_o = np.sqrt(gamma*Rspec*T0) #stagnation speed of sound, stagnation properties
    rho_o = P0/(Rspec*T0) #stagnation density

    #evaluate pressure drop across orifice, add a factor of safety to not experience pressure fluctuations and
    #to achieve 300 psi chamber pressure, pg. 466 in fluid mechanics textbook
    #%%
    #mdot check
    #pg 602 gas dynamics James E. John

    orificeAmin = mdot/((np.sqrt(((2/(gamma - 1)*(P2Goal/P0)**(2/gamma)*((P0/P2Goal)**((gamma - 1)/gamma) - 1))))/(1+(2/np.pi)*(P2Goal/P0)**(1/gamma)))*rho_o*a_o)
    d_orificeMin = np.sqrt((orificeAmin*4)/np.pi)
    d_inchesMin = d_orificeMin/0.0254

    P4_real = fsolve(lambda Pout: - mdot + (np.sqrt(((2/(gamma - 1)*(Pout/P0)**(2/gamma)*((P0/Pout)**((gamma - 1)/gamma) - 1))))/(1+(2/np.pi)*(Pout/P0)**(1/gamma)))*rho_o*a_o*orificeAmin, 2.0684e6)

    #observe any corrections due to standard orifice sizing
    dPurchase = dPurchase*0.0254 #purchase orifice size from mcmaster
    dPurchase = d_orificeMin
    Apurchase = np.pi*(dPurchase/2)**2 #m^2 
    #theoretical maximum from purchased orifice
    mdotPurchase = (((2/(gamma + 1))**((gamma + 1)/(2*(gamma -1)))))/(1 + (2/np.pi)*(gamma/(gamma + 1))**(1/(gamma - 1))) * rho_o* a_o*Apurchase
    #solve for downstream pressure after orifice, want to equal goal P
    P4 = fsolve(lambda P: - mdotPurchase + (np.sqrt(((2/(gamma - 1)*(P/P0)**(2/gamma)*((P0/P)**((gamma - 1)/gamma) - 1))))/(1+(2/np.pi)*(P/P0)**(1/gamma)))*rho_o*a_o*Apurchase, 2.0684e6)
    #if the following downstream pressure is observed in experimentation, this is the
    #flow rate 
    mdotExperimental = (np.sqrt(((2/(gamma - 1)*(P4/P0)**(2/gamma)*((P0/P4)**((gamma - 1)/gamma) - 1))))/(1+(2/np.pi)*(P4/P0)**(1/gamma)))*rho_o*a_o*Apurchase

    #mdot check other #pg. 79 gas dynamics
    #theoretical max based off of stagnation properties and calculated orifice diameter
    mdotmaxO = 0.040418*P0*orificeAmin/np.sqrt(T0) 
    Acrit = mdot*np.sqrt(T0)/(0.040418*P0)
    dcrit = np.sqrt((4*Acrit)/(np.pi))/0.0254 #in

    #calculation of P* for choking condition (relevance unknown)
    Pstar = (2/(gamma + 1))**(gamma/(gamma - 1))*P0
    PdropOrifice = (P_static - P4_real)*1e-6
    if PdropOrifice < 0.0:
        soln[:,ii] = np.zeros((len(soln[:,ii])))
    else:
        PdropSystem = (PdropLines + PdropOrifice)
        print('\nOrifice Results:')
        print('Upstream Pressure:', InitialPressure*1e-6,'(MPa)')
        print('Solved orifice diameter (purchase one close to this)',d_inchesMin,'(in)')
        print('critical orifice diameter',dcrit,'(in)')
        print('Pressure entering chamber',P4_real*1e-6,'(MPa)')
        print('Pressure entering chamber',P4_real*1e-6*145.038,'(PSI)')
        print('Pressure drop across lines', PdropLines,'(MPa)')
        print('Pressure drop across orifice ',PdropOrifice,'(MPa)')
        print('Total Pressure drop across system', PdropSystem,'(MPa)')
        print('mdot from experimentation from purchased orifice', mdotExperimental, '(kg/s)')
        print('mdot theoretical max from calculated downstream pressure', mdotPurchase, '(kg/s)')
        print('mdot theoretical from stagnation properties', mdotmaxO, '(kg/s)\n')
        
        soln[:,ii] =[P4_real,PdropLines,PdropOrifice,PdropSystem,d_inchesMin,mdotExperimental,pressureRange[ii]] 
#write helper function to delete zeros of soln vector
#%%
xlabel = 'Upstream Pressure (MPa)'
soln2 = np.ma.masked_equal(soln,0)

fig, axs = plt.subplots(2, 1, constrained_layout=True)
axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln2[0,:]/1e6,'-',label = 'Pressure after Orifice')
axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln[1,:],'-', label = 'Delta P Across Lines')
axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln[2,:],'-', label = 'Delta P Across Orifice')
axs[0].plot(pDropCalc.psi_to_MPa(soln2[6,:]), soln[3,:],'-', label = 'Total Delta P')
axs[0].legend()
axs[0].set_title('Changes in System Pressures with Upstream Pressure')
axs[0].set_xlabel(xlabel)
axs[0].set_ylabel('Pressure (MPa)')
#fig.suptitle(' Flow Results', fontsize=12)
axs[1].set_ylabel('Orifice Diameter [in]')
lns1 = axs[1].plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[4,:], '-',label = 'Orifice Diameter [in]', color = 'tab:orange')
ax2 = axs[1].twinx()
ax2.set_ylabel('Mass Flowrate [kg/s]')
lns2 = ax2.plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[5,:], 'o',label = 'Mdot Experimental [kg/s]', color = 'tab:blue')
axs[1].set_xlabel(xlabel)
axs[1].set_title('Orifice Properties Against Upstream Pressure')
axs[1].legend(loc='upper right')
ax2.legend(loc='right')

plt.show()
#Explanation of Results

#mdotmax0 is a theoretical maximum flow rate based off a given orifice diameter
#in this case the inner diameter of a 1/4" steel tube

#dcrit is the diameter required to achieve this theoretical flow rate ofmdotmaxO 

#orificeAmin is the minimum area required to choke the flow based on the
#the maximum flow rate calculated from the AirLineAnalysis

#Solved orifice diameter (d_inchesMin) is the diameter required to choke the flow
#at orificeAmin

#P4_real is based from the maximum flow rate allowed in the system from 
#Airline Analysis and predicts the pressure entering the chamber after an orifice
#which is d_inchesMin

#mdotExperimental is the theoretical flow rate achieved based on stagnation
#pressure, static pressure after the orifice (set as ideal) and dPurchase.
#if that downstream pressure is measured it will produce mdotExperimental mdot
#in the current case it is the goal pressure (upstream Pressure/1.528), the
#minimum required to choke the flow across the orifice


#mdotPurchase is a theoretical flow rate ahiceved based off stagnation conditions
#and given orifice diameter (dPurchase)

#P4 is a calculation to determine the downstream pressure in the theoretical case of 
#mdotPurchase instead of setting the downstream pressure




