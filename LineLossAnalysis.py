# -*- coding: utf-8 -*-
"""
Created on Sun May 10 16:11:11 2020

@author: cemam
"""
#%%Imports
import numpy as np
from scipy.optimize import fsolve
import frictionFactor as ff
import cantera as ct
import matplotlib.pyplot as plt
import DifferentialFannoPressureDrop as FannoF

#%%
#Inputting Minor/Major Losses
#NOTE:any value that is zero should be written as 1e-16 to avoid errors

#Minor/Major Loss calculations:
#https://docs.google.com/spreadsheets/d/1gn8gLvc1uLY0blPmABOoTxjrHLD8UaT0t1OvKjpqaj4/edit#gid=89965144
Kvec = [0.75, 0.00001006424457,0.2,1.816115216,7.5, 0.5, 3.30000078]
Lvec = [1e-16, 1e-16, 1.2192, 0.5461, 1.5367, 1.2192, 0.0635]
Dvec = [0.004572, 0.009398, 0.008636, 0.009398, 0.009398, 0.008636,0.009398]
epsvec = [0.015e-3, 0.015e-3, 3e-6, 0.015e-3, 0.001e-3, 3e-6, 0.015e-3]  #mm#https://www.engineeringtoolbox.com/lined-pipe-pressure-loss-d_1178.html ptfe roughness for flex hose
#https://www.enggcyclopedia.com/2011/09/absolute-roughness/ aluminum/steel roughness
def psi_to_MPa(psi):
    return  0.00689*psi


def LineAnalysis(mdot, T0,P1,gas,gamma, step,intermediateOutputs=False,plots=True):
    
    pipeSegments = len(Kvec)
    soln = np.zeros((5,pipeSegments))
    pDropvec = np.zeros(pipeSegments)
    Pvec = np.zeros(pipeSegments)
    Tvec = np.zeros(pipeSegments)
    rhovec = np.zeros(pipeSegments)
    viscvec = np.zeros(pipeSegments)

    for ii in range(0,pipeSegments):
        flowData = FannoF.Fanno(mdot,T0,P1,gas, gamma, Dvec[ii], Lvec[ii], step, epsvec[ii], Kvec[ii],plots = False) #choose step size in this function
        if ii == 0:
            fullArrayFlowData = flowData
            Lvector = np.arange(0,Lvec[ii],step)
        else: 
            fullArrayFlowData = np.append(fullArrayFlowData,flowData,axis=1)
            Lvector =  np.append(Lvector,np.arange(0,Lvec[ii],step)+np.add(Lvector[-1],0))
            
        
        #obtain conditions at the end of the tube    
        P1 = flowData[1,int(Lvec[ii]/step-1)]    
        pDrop1 = flowData[7,int(Lvec[ii]/step-1)]
        T0 = flowData[0,int(Lvec[ii]/step-1)]
        rho = flowData[2,int(Lvec[ii]/step-1)]
        visc = flowData[8,int(Lvec[ii]/step-1)]
        M2 = flowData[5,int(Lvec[ii]/step-1)]
        kReal = flowData[9,int(Lvec[ii]/step-1)]
        k2 = flowData[10,int(Lvec[ii]/step-1)]
               
        
        pDropvec[ii] = pDrop1
        Pvec[ii] = P1
        Tvec[ii] = T0
        rhovec[ii] = rho
        viscvec[ii] = visc
        if intermediateOutputs:
            print('\nSection of analysis',ii)
            print('kReal',kReal)
            print('k2', k2)
            print('M2',M2)
            print('pDrop',pDrop1/1e6)
            print('P1',P1)
            print('T0',T0)
            print('rho',rho)
            print('visc',visc)
    
          
    soln[0,:] = pDropvec
    soln[1,:] = Pvec
    soln[2,:] = Tvec
    soln[3,:] = rhovec
    soln[4,:] = viscvec
    PdropTotal = sum(soln[0,:])
    print('\nTotal Pressure Drop Across System',PdropTotal/1e6,'(MPa)\n')
    
    
    if plots:
     
     xlabel = 'Length of tube [m]'
     fig, axs = plt.subplots(4, 1, constrained_layout=True)
     axs[0].plot(Lvector, fullArrayFlowData[4,:],'-')
     axs[0].set_title('Mach Number Along Tube')
     axs[0].set_xlabel(xlabel)
     axs[0].set_ylabel('Mach Number')
     fig.suptitle('Fanno Flow Results', fontsize=12)
    
     axs[1].plot(Lvector,fullArrayFlowData[6,:], '-')
     axs[1].set_xlabel(xlabel)
     axs[1].set_title('Friction Factor Along Tube')
     axs[1].set_ylabel('Friction Factor')
    
     axs[2].plot(Lvector,fullArrayFlowData[1,:]/1e6, '-')
     axs[2].set_xlabel(xlabel)
     axs[2].set_title('Pressure Along Tube')
     axs[2].set_ylabel('P [mPa]')

     axs[3].plot(Lvector,fullArrayFlowData[0,:], '-')
     axs[3].set_xlabel(xlabel)
     axs[3].set_title('Temperature Along Tube')
     axs[3].set_ylabel('T [K]')
     plt.show()
    return soln,fullArrayFlowData, PdropTotal/1e6
#%%
if __name__ == '__main__':
    mdot = 0.3 #mass flow rate [kg/s]
    T0 = 300# #ambient temperature [K]/stagnation
    stag = T0
    P1 = psi_to_MPa(1500)*1e6 #(100/145.038)*1e6#MPa #pressure at pipe inlet [Pa] 1500 PSI 
    PInit = P1
    #define gas
    gas = ct.Solution('air.cti') 
    gas.transport_model = 'Mix'
    gas.TPX = T0,P1,{'O2':1,'N2':3.76}
    R = 8.314
    gamma = gas.cp/gas.cv #ratio of specific heats [-]
    MW = gas.mean_molecular_weight
    step = 0.001 #(m) break tube up into 1mm increments
    #for loss calculations https://docs.google.com/spreadsheets/d/1gn8gLvc1uLY0blPmABOoTxjrHLD8UaT0t1OvKjpqaj4/edit#gid=0

    LineLoss = LineAnalysis(mdot, T0, P1, gas, gamma, step)
    #print(LineLoss.)
    T = LineLoss[0][2,0] #initial temperature (this is the temperature after one step)
    P = (LineLoss[0][0,-1])/1e6
    PdropTotal = LineLoss[2]
