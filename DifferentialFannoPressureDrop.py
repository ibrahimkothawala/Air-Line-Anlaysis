"""
Created on: 11/16/2019 at 12:08pm
Author: Christopher Emami & Ibrahim Kothawala
Last Updated: 11/16/2019
"""
#%%Imports
import numpy as np
from scipy.optimize import fsolve
import frictionFactor as ff
import matplotlib.pyplot as plt

#%%Inputs
def maxfLoverD(mach,gamma):
    return ((1-(mach**2))/(gamma*(mach**2)))+(((gamma+1)/(2*gamma))*np.log(((gamma+1)*(mach**2))/(2+((gamma-1)*(mach**2)))))

def mach2(mach1,gamma,fLD2):
    return fsolve(lambda mach1: ((1-(mach1**2))/(gamma*(mach1**2)))+(((gamma+1)/(2*gamma))*np.log(((gamma+1)*(mach1**2))/(2+((gamma-1)*(mach1**2)))))-fLD2,0.01)

def Pratio(mach, gamma):
    return (1/mach)*(((gamma+1)/(2+((gamma-1)*(mach**2))))**0.5)

def Pdrop(M1, M2, gamma, P1):
    pRatio1 = Pratio(M1, gamma)
    pRatio2 = Pratio(M2, gamma)
    P2 = (pRatio2/pRatio1)*P1
    return P1-P2

def Tratio(mach1, mach2, gamma):
    tRatio1 = (gamma +1)/(2+(gamma-1)*mach1**2)
    tRatio2 = (gamma +1)/(2+(gamma-1)*mach2**2)
    return tRatio2/tRatio1

def rhoRatio(mach1, mach2, gamma):
    rhoRatio1 = (1/mach1)*(((2+((gamma-1)*(mach1**2)/(gamma+1))))**0.5)
    rhoRatio2 = (1/mach2)*(((2+((gamma-1)*(mach2**2)/(gamma+1))))**0.5)
    return rhoRatio2/rhoRatio1



def Fanno(mdot, T0, P1, gas, gamma, d, L, step, eps, sumK,plots=True):
    
    #define gas    
    #define properties of gas
    MW = gas.mean_molecular_weight
    rho = gas.density_mass
    visc = gas.viscosity
    R = 8.314
    Rspec = R/MW*1000
    
    #Fanno Flow
    u1 = mdot/(rho*np.pi*((d/2)**2)) #initial velocity
    fF = ff.f(d,eps,(rho*u1*d)/visc)/4
    Lvector = np.arange(0.0,L,step)
    #end conditions
    Tvec = np.zeros(len(Lvector))
    P1vec = np.zeros(len(Lvector))
    P2vec = np.zeros(len(Lvector))
    Uvec = np.zeros(len(Lvector))
    rhovec = np.zeros(len(Lvector))
    M1vec =  np.zeros(len(Lvector))
    M2vec =  np.zeros(len(Lvector))
    fLDmaxvec =  np.zeros(len(Lvector))
    fFvec =  np.zeros(len(Lvector))
    pDropvec = np.zeros(len(Lvector))
    fLD1Total = np.zeros(len(Lvector))
    viscvec = np.zeros(len(Lvector))
    fLD2vec =  np.zeros(len(Lvector))
    soln = np.zeros((11,len(Lvector)))
    
    pDropTotal = 0
    T1 = T0
    M1 = u1/np.sqrt(gamma*Rspec*T1)
    
    if M1 > 0.98: #catch when the flow chokes and do not calculate that solution for faster run time
        raise Exception("Flow is Fanno choked, upstream pressure and mass flow rate combo invalid")
    
    fLDmax = maxfLoverD(M1,gamma)  #constant fL/D based on gamma and M1
    rho1 = rho
    sumKparsed = 0
    for i in range(0,len(Lvector)):
        
        if len(Lvector) == 1: # if there is no length of tube for example when analyzing a valve
            fLD1 = sumK
        else: 
            fLD1 = ((4*fF*Lvector[i])/d) + sumKparsed
        
        fLD1Total[i] = fLD1
        M1 = u1/np.sqrt(gamma*Rspec*T1)
        fLD2 = (fLDmax-fLD1)
        M2 = mach2(M1,gamma,fLD2)
        pDrop = Pdrop(M1, M2, gamma, P1)
        pDropTotal = pDropTotal + pDrop
        T2 = (Tratio(M1, M2, gamma))*T1
        rho2 = rhoRatio(M1,M2,gamma)*rho1

        #reset values
        u1 = mdot/(rho2*np.pi*((d/2)**2))
        rho1 = rho2
        T1 = T2
        P1vec[i] = P1
        P1 -=pDrop #P2 at the end of the tube being fed into the next iteration
        gas.TP = T1,P1
        visc1 = gas.viscosity
        fF = ff.f(d,eps,(rho2*u1*d)/visc1)/4
        sumKparsed = sumKparsed + sumK/len(Lvector)
        
        Tvec[i] = T2
        P2vec[i] = P1
        rhovec[i] = rho2
        Uvec[i] = u1
        M1vec[i] = M1
        M2vec[i] = M2
        fLDmaxvec[i] = fLDmax
        fFvec[i] = fF
        pDropvec[i] = pDropTotal
        viscvec[i] = visc1
        fLD2vec[i] = fLD2

        soln[0,:] = Tvec
        soln[1,:] = P2vec
        soln[2,:] = rhovec
        soln[3,:] = Uvec
        soln[4,:] = M1vec
        soln[5,:] = M2vec
        soln[6,:] = fFvec
        soln[7,:] = pDropvec
        soln[8,:] = viscvec
        soln[9,:] = fLD1Total
        soln[10,:] = fLD2vec
        
        
    if plots:    
        xlabel = 'Length of tube [m]'
        fig, axs = plt.subplots(3, 1, constrained_layout=True)
        axs[0].plot(Lvector, M1vec,'-')
        axs[0].set_title('Mach Number Along Tube')
        axs[0].set_xlabel(xlabel)
        axs[0].set_ylabel('Mach Number')
        fig.suptitle('Fanno Flow Results', fontsize=12)
    
        axs[1].plot(Lvector,fFvec, '-')
        axs[1].set_xlabel(xlabel)
        axs[1].set_title('Friction Factor Along Tube')
        axs[1].set_ylabel('Friction Factor')
    
        axs[2].plot(Lvector,pDropvec/1e-6, '-')
        axs[2].set_xlabel(xlabel)
        axs[2].set_title('Pressure Drop Along Tube')
        axs[2].set_ylabel('Delta P [mPa]')
        plt.show()
    return soln



