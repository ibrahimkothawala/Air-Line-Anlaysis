#author: Ibrahim Kothawala 

import numpy as np

#inputs 
#outputs
def PchokeCondition(gamma): #pressure ratio of Pthroat/PdownstreamStagnation to velocity choke flow
    return (2/(gamma + 1))**(gamma/(gamma-1))

#Outputs
#inputs
