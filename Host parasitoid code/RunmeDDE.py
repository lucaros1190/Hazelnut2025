
# Script to solve the delayed host-parasitoid model
#
# Created by Luca Rossini on 27 January 2025
# Last update 9 April 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

from Parameters import *
from Functions import *
import numpy as np
import scipy.integrate as sc 
from ddeint import ddeint
import matplotlib.pyplot as plt
import time


# Timer start

t_start = time.time()


# Acquisition of the time span: consider that the Delta t subdivision affects the precision of the algorithm. Make some attempts to find a good compromise. The Delta t is the variable IntegrationStep, to define according to the needs of the user.

IntegrationStep = 10000

t_span = np.linspace(0, len(DailyTemp), IntegrationStep)
print(t_span)

# Define the initial history for the DDE: change the array inside the function to modify the initial condition

def history(t):
    #s = 2 * np.heaviside(t, 0.5) * np.array([100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 10])
    s = 2 * np.heaviside(t, 0.5) * np.array([10000, 10000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100])
    return s


# Compute the normalisation coefficients for the development rates and adult survival

[NormCoeff_Egg_H, NormCoeff_N1_H, NormCoeff_N2_H, NormCoeff_N3_H, NormCoeff_N4_H, NormCoeff_N5_H, NormCoeff_AM_H, NormCoeff_AF_H, NormCoeff_Egg_P, NormCoeff_L1_P, NormCoeff_L2_P, NormCoeff_L3_P, NormCoeff_P_P, NormCoeff_AF_P, NormCoeff_AM_P] = NormCoefficients(DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P)


# Solve the equation using 'ddeint' - Normalised rates

print('\n1. Solving DDE model WITH normalisation coefficients \n')

#SolDDE_Normalised = ddeint(DDE_FuncNormalised, history, t_span, fargs = (FertPar_H, MortPar_Egg_H, MortPar_N1_H, MortPar_N2_H, MortPar_N3_H, MortPar_N4_H, MortPar_N5_H, MortPar_Egg_P, MortPar_L1_P, MortPar_L2_P, MortPar_L3_P, DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P, NormCoeff_Egg_H, NormCoeff_N1_H, NormCoeff_N2_H, NormCoeff_N3_H, NormCoeff_N4_H, NormCoeff_N5_H, NormCoeff_AM_H, NormCoeff_AF_H, NormCoeff_Egg_P, NormCoeff_L1_P, NormCoeff_L2_P, NormCoeff_L3_P, NormCoeff_P_P, NormCoeff_AF_P, NormCoeff_AM_P, DailyTemp, SR_H, SR_P, a, T_h, b))
SolDDE_Complete = ddeint(DDE_FuncComplete, history, t_span, fargs = (FertPar_H, MortPar_Egg_H, MortPar_N1_H, MortPar_N2_H, MortPar_N3_H, MortPar_N4_H, MortPar_N5_H, MortPar_Egg_P, MortPar_L1_P, MortPar_L2_P, MortPar_L3_P, DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P, NormCoeff_Egg_H, NormCoeff_N1_H, NormCoeff_N2_H, NormCoeff_N3_H, NormCoeff_N4_H, NormCoeff_N5_H, NormCoeff_AM_H, NormCoeff_AF_H, NormCoeff_Egg_P, NormCoeff_L1_P, NormCoeff_L2_P, NormCoeff_L3_P, NormCoeff_P_P, NormCoeff_AF_P, NormCoeff_AM_P, DailyTemp, SR_H, SR_P, a, T_h, b))


# Timer stop

t_f = time.time() - t_start

print('\nExecution time:', t_f, 's\n')
    

# Plot the results

    # Egg stage
'''   
plt.figure(1)

plt.plot(t_span, SolDDE_Normalised[:,0], color = 'orange', label = 'Egg - Host')
plt.plot(t_span, SolDDE_Normalised[:,1], color = 'blue', label = 'Egg parasitised - Host')
plt.plot(t_span, SolDDE_Normalised[:,10], color = 'red', label = 'Egg - Parasitoid')
plt.legend()


    # Overall solutions overlapped - Model WITH normalised development rates

plt.figure(2)

plt.plot(t_span, SolDDE_Normalised[:, 0], label = 'Egg - Host')
plt.plot(t_span, SolDDE_Normalised[:, 1], label = 'Parasitised eggs')
plt.plot(t_span, SolDDE_Normalised[:, 2], label = 'N1 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 3], label = 'N2 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 4], label = 'N3 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 5], label = 'N4 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 6], label = 'N5 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 7], label = 'Ad male - Host')
plt.plot(t_span, SolDDE_Normalised[:, 8], label = 'Ad Female 1 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 9], label = 'Ad Female 2 - Host')
plt.plot(t_span, SolDDE_Normalised[:, 10], linestyle = 'dashed', label = 'Egg - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 11], linestyle = 'dashed', label = 'L1 - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 12], linestyle = 'dashed', label = 'L2 - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 13], linestyle = 'dashed', label = 'L3 - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 14], linestyle = 'dashed', label = 'Pupa - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 15], linestyle = 'dashed', label = 'Ad male - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 16], linestyle = 'dashed', label = 'Ad Female 1 - Parasitoid')
plt.plot(t_span, SolDDE_Normalised[:, 17], linestyle = 'dashed', label = 'Ad Female 2 - Parasitoid')
plt.legend()


plt.show()
'''

    # Egg stage

plt.figure(1)

plt.plot(t_span, SolDDE_Complete[:,0] + SolDDE_Complete[:,1], color = 'orange', label = 'Egg - Host')
plt.plot(t_span, SolDDE_Complete[:,0], color = 'yellow', label = 'Egg - Host no transient')
plt.plot(t_span, SolDDE_Complete[:,15] + SolDDE_Complete[:,16], color = 'red', label = 'Egg - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:,15], color = 'blue', label = 'Egg - Parasitoid no transient')
plt.legend()
#plt.ylim(0, 10000)


    # Overall solutions overlapped - Model WITH normalised development rates

plt.figure(2)

plt.plot(t_span, SolDDE_Complete[:, 0] + SolDDE_Complete[:, 1], label = 'Egg - Host')
plt.plot(t_span, SolDDE_Complete[:, 2] + SolDDE_Complete[:, 3], label = 'N1 - Host')
plt.plot(t_span, SolDDE_Complete[:, 4] + SolDDE_Complete[:, 5], label = 'N2 - Host')
plt.plot(t_span, SolDDE_Complete[:, 6] + SolDDE_Complete[:, 7], label = 'N3 - Host')
plt.plot(t_span, SolDDE_Complete[:, 8] + SolDDE_Complete[:, 9], label = 'N4 - Host')
plt.plot(t_span, SolDDE_Complete[:, 10] + SolDDE_Complete[:, 11], label = 'N5 - Host')
plt.plot(t_span, SolDDE_Complete[:, 12], label = 'Ad male - Host')
plt.plot(t_span, SolDDE_Complete[:, 13], label = 'Ad Female 1 - Host')
plt.plot(t_span, SolDDE_Complete[:, 14], label = 'Ad Female 2 - Host')
plt.plot(t_span, SolDDE_Complete[:, 15] + SolDDE_Complete[:, 16], linestyle = 'dashed', label = 'Egg - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 17] + SolDDE_Complete[:, 18], linestyle = 'dashed', label = 'L1 - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 19] + SolDDE_Complete[:, 20], linestyle = 'dashed', label = 'L2 - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 21] + SolDDE_Complete[:, 22], linestyle = 'dashed', label = 'L3 - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 23] + SolDDE_Complete[:, 24], linestyle = 'dashed', label = 'Pupa - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 25], linestyle = 'dashed', label = 'Ad male - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 26], linestyle = 'dashed', label = 'Ad Female 1 - Parasitoid')
plt.plot(t_span, SolDDE_Complete[:, 27], linestyle = 'dashed', label = 'Ad Female 2 - Parasitoid')
plt.legend()
#plt.ylim(0, 10000)


    # Summary plot for the manuscript

fig, graphres = plt.subplots(nrows = 2)
fig.suptitle('Host-parasitoid dynamics: Simulation 2')

graphres[0].plot(t_span, SolDDE_Complete[:,0] + SolDDE_Complete[:,1], color = 'orange', label = 'Host - Egg')
graphres[0].plot(t_span, SolDDE_Complete[:,15] + SolDDE_Complete[:,16], color = 'red', linestyle = 'dashed', label = 'Parasitoid - Egg')
graphres[0].set(ylabel = 'Number of individuals')
graphres[0].legend()

graphres[1].plot(t_span, SolDDE_Complete[:,13] + SolDDE_Complete[:,14], color = 'purple', label = 'Host - Adult females')
graphres[1].plot(t_span, SolDDE_Complete[:,26] + SolDDE_Complete[:,27], color = 'grey', linestyle = 'dashed', label = 'Parasitoid - Adult Females')
graphres[1].set(xlabel = 'Time (days)')
graphres[1].set(ylabel = 'Number of individuals')
graphres[1].legend()

plt.show()