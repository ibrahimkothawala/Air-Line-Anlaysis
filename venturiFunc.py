#author: Ibrahim Kothawala 

import numpy as np

#Gas Constant
R = 8.31446261815324 #j/K/mol

#inputs: gamma
#outputs: Pstar/P0
#References: https://en.wikipedia.org/wiki/Choked_flow
#pg 601 eqn 15.46 Gas Dynamics James John 
def PchokeCondition(gamma): #pressure ratio of Pthroat/PdownstreamStagnation to velocity choke flow
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
    return (A_throat*np.sqrt((2*gamma*rho_ups*pStaticUps)/(gamma-1)*((pStatThroat/pStaticUps)**(2/gamma)-(pStatThroat/pStaticUps)**((gamma+1)/gamma))))/np.sqrt(1-(A_throat/A_ups)**2*(pStatThroat/pStaticUps)**(2/gamma))
    