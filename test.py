  
#%%
import numpy as np
from scipy.optimize import fsolve
import frictionFactor as ff
import cantera as ct
import matplotlib.pyplot as plt
import DifferentialFannoPressureDrop as DFP
import Orifice_size as ORS
import LineLossAnalysis as LL

p1 = LL.psi_to_MPa(300)*1e6
p2 = LL.psi_to_MPa(2500)*1e6
#%% Inputs
pts = 30
mdot = 0.3 #mass flow rate [kg/s]
T0 = 300# #ambient temperature [K]/stagnation (stagnation temperature remains constant in fanno flow see https://en.wikipedia.org/wiki/Fanno_flow)
P1_range = np.linspace(p1,p2,pts) #static pressure before orifice [Pa]
R = 8.314
D_1 = 0.004572 #tube diameter (m) out of pressure regulator is 1/4" flex tube

#%% Calculated Inputs
Ainit = np.pi*(D_1/2)**2 #area of the tube after the regulator
soln = np.zeros((pts,6))


#%% Cantera
gas = ct.Solution('air.cti') 
gas.transport_model = 'Mix'
MW = gas.mean_molecular_weight
Rspec = R/MW*1000 #J/kg-K



#%% Main Loop
for ii in range(pts):
    #gas initial properties
    gas.TPX = T0,P1_range[ii],{'O2':1,'N2':3.76}
    gamma = gas.cp/gas.cv
    rho = gas.density_mass
    u = mdot/(rho*Ainit)
    dynamicP = 0.5*rho*u**2
  
    #stagnation properties
    #P0 =  P1_range[ii] + dynamicP #Pascals stagnation pressure
    P0 = ORS.P02vec[ii]
    a_o = np.sqrt(gamma*Rspec*T0) #stagnation speed of sound, stagnation properties
    rho_o = P0/(Rspec*T0) #stagnation density

    #Testing theoretical orifice sizing funcs
    det_P4 = P1_range[ii]*ORS.PchokeCondition(gamma) #determined pressure after orifice from 
    theo_area = ORS.orificeThroatArea(gamma,rho_o,a_o,P0,det_P4,mdot)
    calc_P4 = ORS.pressureAfterOrifice(gamma, rho_o, a_o, P0, mdot, theo_area, det_P4)
    
    #testing maximum mass flowrate functions
    mdot_thru_theo = ORS.massFlowRateThruOrifice(gamma,rho_o, a_o,P0,det_P4,theo_area) #massflowrate thru theoreical orifice area at choking conditions
    mdot_max_thru_theo = ORS.maxMassFlowRate(gamma,rho_o,a_o,theo_area)
    
    #putting parameters into solution array
    soln[ii] = theo_area,det_P4,calc_P4,mdot_thru_theo,mdot_max_thru_theo, P0

#%% Plotting Helper Variables
det_Pratio = np.divide(P1_range,soln[:,1])
calc_Pratio = np.divide(P1_range,soln[:,2])
mdotRatio_act = soln[:,3]/mdot
mdotRatio_max = soln[:,4]/mdot
P1_range_mPa = P1_range/1e6
#%% Plotting results from testing theoretical orifice sizing funcs
fig, axs = plt.subplots(3, 1,constrained_layout=True)
#aspectRatio = 4
xlabel = "Upstream Pressure [mPa]"

##Sub-plot of downstream pressure ratios
#axs[0].plot(P1_range_mPa,det_Pratio,'o',label = "Determined")
#axs[0].plot(P1_range_mPa,calc_Pratio,'x', label = "Calculated")
#axs[0].set_title('Comparison of 2 different methods of calculating downstream pressure \n and showing they are intrinsically the same')
#axs[0].set_ylabel('Upstream to Downstream \n Pressure ratio')
#axs[0].set_xlabel(xlabel)
##axs[0].set_aspect(aspectRatio)
#axs[0].legend()

#sub-plot of orifice area
axs[0].plot(P1_range_mPa,soln[:,0],label = "theoretical orifice area")
axs[0].set_title("Theoretical Orifice size against upstream pressure \n")
axs[0].set_ylabel("orifice area [m^2]")
axs[0].set_xlabel(xlabel)
#axs[1].set_aspect(aspectRatio)

#sub-plot of mdot ratios
axs[1].plot(P1_range_mPa,mdotRatio_act, label = "ratio mdot specified to calculated from pressures")
axs[1].plot(P1_range_mPa,mdotRatio_max, label = "ratio mdot specified to calculated maximum")
axs[1].set_title("Ratio of M dot specified to calculated \n Proves they are the same")
axs[1].set_xlabel(xlabel)
axs[1].set_ylabel("mdot ratio")
axs[1].legend()

#sub-plot of orifice area
axs[2].plot(P1_range_mPa,soln[:,5],label = "stagnation pressure")
axs[2].set_title("Stagnation pressure against upstream pressure \n")
axs[2].set_ylabel("orifice area [m^2]")
axs[2].set_xlabel(xlabel)
#axs[2].set_aspect(aspectRatio)
plt.tight_layout()
plt.show()




# fig, axs = plt.subplots(1, 1, constrained_layout=True)
# axs.plot(P1_range,det_Pratio,'o',label = 'Determined')
# axs.plot(P1_range,calc_Pratio,'x',label = 'Calculated')
# axs.legend()
# axs.set_title('Comparing results from determined Pressure after orifice and calculated')
# axs.set_xlabel('Orifice Upstream Pressure [Pa] ')
# axs.set_ylabel('Pressure ratio')

# #fig.suptitle(' Flow Results', fontsize=12)
# axs[1].set_ylabel('Orifice Diameter [m]')
# lns1 = axs[1].plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[4,:], '-',label = 'Theoretical Orifice Diameter', color = 'tab:orange')
# ax2 = axs[1].twinx()
# ax2.set_ylabel('Orifice Area')
# lns2 = ax2.plot(pDropCalc.psi_to_MPa(soln2[6,:]),soln2[5,:], 'o',label = 'Mdot Experimental [kg/s]', color = 'tab:blue')
# axs[1].set_xlabel(xlabel)
# axs[1].set_title('Orifice Properties Against Upstream Pressure')
# axs[1].legend(loc='upper right')
# ax2.legend(loc='right')



    



# %%