#author: Ibrahim Kothawala 

import numpy as np

#inputs: gamma
#outputs: Pstar/P0
#References: https://en.wikipedia.org/wiki/Choked_flow
def PchokeCondition(gamma): #pressure ratio of Pthroat/PdownstreamStagnation to velocity choke flow
    return (2/(gamma + 1))**(gamma/(gamma-1))

#Outputs
#inputs
