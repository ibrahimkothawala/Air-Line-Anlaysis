# Air Line Analysis
This calculator is used to determine the pressure drop of the gaseous oxidizer in our propellent delivery system and to size the orifice which acts as the oxidizer injector. 

# Required Libraries
-Numpy
-Cantera
-Matplotlib
-Scipy

## How does Air Line Analysis calculate pressure drop?
It utilizes the fanno flow model for compressible gasses in pipes. A compressible model was utilized instead of an incompressible model in order to capture the density changes that would occur as the pressure within the system dropped due to friction. In addition Fanno Flow assumes the flow to be adiabatic. This model assumption is the most likely to fail. To learn more about fanno flow and the governing equations read chapter 6 of The Dynamics and Thermodynamics of COMPRESSIBLE FLUID FLOW by Ascher Shapiro. The pressure drop is calculated using DifferentialFannoPressureDrop.py. LineLossAnalysis utilizes DifferentialFannoPressureDrop.py to calculate the pressure drop across the entire pipe system. The pressure drop across the orifice utilizes a different model which will be discussed in a further section. 

### Overall Algorithm
The algorithm for determining a pressure drop requires both an initial mass flow rate and pressure. As the code marches down the pipe system it adds up the pressure drops across each element and at the end determines the pressure drop across the entire system. However, sometimes the mass flow rate and initial pressure result in the Mach number reaching or exceeding 1 before the exit of the system. In that case the flow is considered fanno choked and any result is deemed non real and an exception is raised for that mass flow rate and initial pressure combination. Either the initial pressure needs to be increased to achieve that mass flow rate or the mass flow rate needs to be decreased. 

### Calculating gas thermodynamic and transport properties:
Density, viscosity, and molecular weight are needed in order to perform the fanno flow analysis. Since the density and viscosity are not constant along the pipe, Cantera is utilized to calculate those properties at each element in the system. 

### Calculating friction factor
In order to calculate the pressure drop in pipes due to skin friction the friction factor must be calculated. This is done in frictionFactor.py where both an approximation of the Colebrook equation and the Colebrook equation itself are used to determine the friction factor. A log-log plot of the friction factor at different reynolds numbers and pipe roughness is done to show model agreement and comparison to the moody diagram.

### Pressure drop in the pipe segments
The pipe segment is broken into 1 mm long pipe segments. This is done because as the gas moves down the pipe the pressure increases causing the velocity and the friction factor to increase. 1 mm was chosen because at that length the pressure drop and resultant friction factor increase was deemed negligible. First fL/D max is calculated then the fL/D actual is calculated, for the pipe segment, the difference between them is used to calculate the pressure, velocity, and density at the end of the segment.

### Pressure drop in a valve or fitting
Since fL/D is equivalent to the a valve or fittings coefficient K, that is substituted into where fL/D actual goes. Then the algorithm continues as it would in the pressure drop across a pipe segment. 

## How does Air Line Analysis size an Orifice? 
First a theoretical orifice area is determined by setting the ratio of the upstream pressure to the downstream pressure to be 1.528 using eqn 15.47 in Gas Dynamics by James E. John. Pressure ratios higher than this will create choked flow (find ref). The pressure upstream of the orifice is determined using the fanno flow analysis and is the exit pressure of the piping system. Once the area is determined, the pressure downstream of the theoretical orifice is solved for using eqn 15.45 in Gas Dynamics. Since the theoretical area may not be available to purchase the closest size smaller than the theoretical area that can be purchased should be inputted manually. It is already done for our case. The maximum mass flow rate that can be achieved by using the purchased orifice is then calculated using eqn 15.47 in Gas Dynamics and the pressure downstream of the orifice is calculated using eqn 15.45 using that maximum mass flow rate. Then using the downstream pressure calculated in the previous step a mass flow rate is determined. 
