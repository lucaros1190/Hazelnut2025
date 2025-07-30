
# List of parameters needed by RunmeDDE.py to solve the model.
#
# Created by Luca Rossini on 17 January 2025
# Last update 10 April 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

import pandas as pd
import numpy as np


# Absorb the temperature array needed by the model

InputFile = pd.read_excel('TemperatureInput.xlsx')
DailyTemp = InputFile['Average_Temperature']


# Parasitisation parameters - Functional response

    # Maximum predation rate
    
a = 0.01#10**(-12)

    # Handling time
    
T_h = 0.05

    # Shape parameter of the type 3 response
    
b = 3


###################################################
#                                                 #
#       HOST PARAMETERS - Halyomorpha halys       #
#                                                 #
###################################################


    # Sex ratio - Host

SR_H = 0.5


    # Fertility - Host

alpha_H = 74968.6851351
gamma_H = 10.0000002
Lambda_H = 21.1615365
delta_H = 5.0865985
tau_H = 25.0653308

FertPar_H = [alpha_H, gamma_H, Lambda_H, delta_H, tau_H]


    # Mortality - Host EGGS

a_e_H = 0.0001667
b_e_H = -0.0162653
c_e_H = 0.5844314
d_e_H = -9.163731
e_e_H = 52.9093736

MortPar_Egg_H = [a_e_H, b_e_H, c_e_H, d_e_H, e_e_H]


    # Mortality - Host N1

a_N1_H = 3.01e-05
b_N1_H = -0.003706
c_N1_H = 0.1668317
d_N1_H = -3.2539274
e_N1_H = 23.2205487

MortPar_N1_H = [a_N1_H, b_N1_H, c_N1_H, d_N1_H, e_N1_H]


    # Mortality - Host N2

a_N2_H = 5.9e-05
b_N2_H = -0.0057541
c_N2_H = 0.2135545
d_N2_H = -3.5937791
e_N2_H = 23.3620976

MortPar_N2_H = [a_N2_H, b_N2_H, c_N2_H, d_N2_H, e_N2_H]


    # Mortality - Host N3

a_N3_H = 0.0001129
b_N3_H = -0.0113512
c_N3_H = 0.4254478
d_N3_H = -7.0417152
e_N3_H = 43.5214056

MortPar_N3_H = [a_N3_H, b_N3_H, c_N3_H, d_N3_H, e_N3_H]


    # Mortality - Host N4

a_N4_H = 0.000156
b_N4_H = -0.015777
c_N4_H = 0.591451
d_N4_H = -9.7212598
e_N4_H = 59.0601027

MortPar_N4_H = [a_N4_H, b_N4_H, c_N4_H, d_N4_H, e_N4_H]


    # Mortality - Host N5

a_N5_H = 6.24e-05
b_N5_H = -0.006405
c_N5_H = 0.2495236
d_N5_H = -4.3471306
e_N5_H = 28.4991245

MortPar_N5_H = [a_N5_H, b_N5_H, c_N5_H, d_N5_H, e_N5_H]


    # Development - Host EGGS
    
a_Egg_H = 0.000238
T_L_Egg_H = 11.03
T_M_Egg_H = 38
m_Egg_H = 2.812

DevPar_Egg_H = [a_Egg_H, T_L_Egg_H, T_M_Egg_H, m_Egg_H]


    # Development - Host N1

a_N1_H = 0.000198
T_L_N1_H = 13.79
T_M_N1_H = 39.41
m_N1_H = 1.892

DevPar_N1_H = [a_N1_H, T_L_N1_H, T_M_N1_H, m_N1_H]


    # Development - Host N2

a_N2_H = 0.0002009
T_L_N2_H = 7.299
T_M_N2_H = 36.01
m_N2_H = 15.32

DevPar_N2_H = [a_N2_H, T_L_N2_H, T_M_N2_H, m_N2_H]


    # Development - Host N3

a_N3_H = 0.00004836
T_L_N3_H = 14.28
T_M_N3_H = 39.93
m_N3_H = 1.106

DevPar_N3_H = [a_N3_H, T_L_N3_H, T_M_N3_H, m_N3_H]


    # Development - Host N4

a_N4_H = 0.00008691
T_L_N4_H = 7.464
T_M_N4_H = 38.98
m_N4_H = 2.011

DevPar_N4_H = [a_N4_H, T_L_N4_H, T_M_N4_H, m_N4_H]


    # Development - Host N5

a_N5_H = 0.00004418
T_L_N5_H = 9.265
T_M_N5_H = 39.99
m_N5_H = 1.49

DevPar_N5_H = [a_N5_H, T_L_N5_H, T_M_N5_H, m_N5_H]


    # Development - Host Adult males and females

a_A_H = 0.00002012
T_L_A_H = 11.0
T_M_A_H = 39.97
m_A_H = 2.097

DevPar_AM_H = [a_A_H, T_L_A_H, T_M_A_H, m_A_H]
DevPar_AF_H = [a_A_H, T_L_A_H, T_M_A_H, m_A_H]


# Lag functions parameters

# ATTENTION PLEASE: if you change the parameters because you want to apply
# this code to a different species, check the if/else filters into the
# specific lag functions!


    # Host Egg lag

a_Egg_H = 0.917942
b_Egg_H = -57.757927
c_Egg_H = 984.505134
d_Egg_H = 1.711403

LagPar_Egg_H = [a_Egg_H, b_Egg_H, c_Egg_H, d_Egg_H]


    # Host N1 lag

a_N1_H = -0.048521
b_N1_H = 2.064456
c_N1_H = 28.760370
d_N1_H = -12.939449

LagPar_N1_H = [a_N1_H, b_N1_H, c_N1_H, d_N1_H]


    # Host N2 lag

a_N2_H = 64.011770
b_N2_H = -0.086601

LagPar_N2_H = [a_N2_H, b_N2_H]


    # Host N3 lag

a_N3_H = -0.044249
b_N3_H = 2.357864
c_N3_H = -27.126639

LagPar_N3_H = [a_N3_H, b_N3_H, c_N3_H]


    # Host N4 lag

a_N4_H = 0.002183
b_N4_H = -0.179449
c_N4_H = 4.514058
d_N4_H = -29.120932

LagPar_N4_H = [a_N4_H, b_N4_H, c_N4_H, d_N4_H]


    # Host N5 lag

a_N5_H = -2.834640
b_N5_H = 250.842032
c_N5_H = 1.626221

LagPar_N5_H = [a_N5_H, b_N5_H, c_N5_H]

    
    # Adult non mated females lag - Preoviposition period!

a_PreOvi_H = 7.0e+15
b_PreOvi_H = -11.147985
c_PreOvi_H = 10.613228

LagPar_PreOvi_H = [a_PreOvi_H, b_PreOvi_H, c_PreOvi_H]


#############################################################
#                                                           #
#       PARASITOID PARAMETERS - Trissolcus japonicus        #
#                                                           #
#############################################################


    # Sex ratio - Parasitoid

SR_P = 0.5


    # Fertility rate - Parasitoid

alpha_P = 37000.000000
gamma_P = 17.761392
Lambda_P = 26.036039
tau_P = 4.000000
delta_P = 22.944989

FertPar_P = [alpha_P, gamma_P, Lambda_P, tau_P, delta_P]


    # Development Egg - Lactin-2 Parasitoid

a_Egg_P = 0.138374
TM_Egg_P = 33.944648
DeltaT_Egg_P = 6.954508
Lambda_Egg_P = -0.289941

DevPar_Egg_P = [a_Egg_P, TM_Egg_P, DeltaT_Egg_P, Lambda_Egg_P]


    # Development L1 - Logan Parasitoid
    
psi_L1_P = 0.085063
rho_L1_P = 0.092881
TM_L1_P = 32.619411
DeltaT_L1_P = 1.553932

DevPar_L1_P = [psi_L1_P, rho_L1_P, TM_L1_P, DeltaT_L1_P]


    # Development L2 - Logan Parasitoid
    
psi_L2_P = 0.043368
rho_L2_P = 0.116276
TM_L2_P = 32.329243
DeltaT_L2_P = 2.217094

DevPar_L2_P = [psi_L2_P, rho_L2_P, TM_L2_P, DeltaT_L2_P]


    # Development L3 - Sharpe and De Michele Parasitoid

A_L3_P = 4.186140
B_L3_P = -130.733316
C_L3_P = 35.404626
D_L3_P = 713.394861
E_L3_P = 7.183079
F_L3_P = -148.630299

DevPar_L3_P = [A_L3_P, B_L3_P, C_L3_P, D_L3_P, E_L3_P, F_L3_P]
    
    
    # Development P - Lactin-2 Parasitoid

a_P_P = 0.032619
TM_P_P = 45.052977
DeltaT_P_P = 14.029446
Lambda_P_P = -1.001966

DevPar_P_P = [a_P_P, TM_P_P, DeltaT_P_P, Lambda_P_P]

    
    # Survival Ad Males - Lactin-1 Parasitoid

a_AM_P = 0.106800
TM_AM_P = 34.999343
DeltaT_AM_P = 9.333671
    
DevPar_AM_P = [a_AM_P, TM_AM_P, DeltaT_AM_P]

    
    # Survival Ad Females - Lactin-1 Parasitoid

a_AF_P = 0.097196
TM_AF_P = 35.000000
DeltaT_AF_P = 10.243270

DevPar_AF_P = [a_AF_P, TM_AF_P, DeltaT_AF_P]


    # Mortality Egg - Parasitoid

a_Mort_Egg_P = 0.000012
b_Mort_Egg_P = -0.000801
c_Mort_Egg_P = 0.021179
d_Mort_Egg_P = -0.299753
e_Mort_Egg_P = 2.000001

MortPar_Egg_P = [a_Mort_Egg_P, b_Mort_Egg_P, c_Mort_Egg_P, d_Mort_Egg_P, e_Mort_Egg_P]

    
    # Mortality L1 - Parasitoid

a_Mort_L1_P = 0.000032
b_Mort_L1_P = -0.002357
c_Mort_L1_P = 0.063728
d_Mort_L1_P = -0.749858
e_Mort_L1_P = 3.450317

MortPar_L1_P = [a_Mort_L1_P, b_Mort_L1_P, c_Mort_L1_P, d_Mort_L1_P, e_Mort_L1_P]
    
    
    # Mortality L2 - Parasitoid
    
a_Mort_L2_P = 0.000046
b_Mort_L2_P = -0.003554
c_Mort_L2_P = 0.100508
d_Mort_L2_P = -1.217428
e_Mort_L2_P = 5.363510

MortPar_L2_P = [a_Mort_L2_P, b_Mort_L2_P, c_Mort_L2_P, d_Mort_L2_P, e_Mort_L2_P]
    
    
    # Mortality L3 - Parasitoid

a_Mort_L3_P = 0.000047
b_Mort_L3_P = -0.003659
c_Mort_L3_P = 0.103070
d_Mort_L3_P = -1.247460
e_Mort_L3_P = 5.457718
    
MortPar_L3_P = [a_Mort_L3_P, b_Mort_L3_P, c_Mort_L3_P, d_Mort_L3_P, e_Mort_L3_P]


# Lag functions parameters - Parasitoid

    # Egg lag - Parasitoid

a_Egg_P = 0.7677908
b_Egg_P = 23.2738017
c_Egg_P = -0.2085237

LagPar_Egg_P = [a_Egg_P, b_Egg_P, c_Egg_P]

    # L1 lag - Parasitoid
        
a_L1_P = 0.6654668
b_L1_P = 15.1767956
c_L1_P = -0.1329154

LagPar_L1_P = [a_L1_P, b_L1_P, c_L1_P]

    # L2 lag - Parasitoid
        
a_L2_P = 0.3328681
b_L2_P = 11.250843
c_L2_P = -0.0902007

LagPar_L2_P = [a_L2_P, b_L2_P, c_L2_P]

    # L3 lag - Parasitoid

a_L3_P = -4.7934135
b_L3_P =  16.917582
c_L3_P = -0.0414033

LagPar_L3_P = [a_L3_P, b_L3_P, c_L3_P]

    # Lag pupa - Parasitoid
        
a_p_P = 4.2511857
b_p_P = 1281.5844394
c_p_P = -0.391008

LagPar_P_P = [a_p_P, b_p_P, c_p_P]

    # Lag AF2 - Parasitoid
        
a_AF2_P = 0.0271894
b_AF2_P = 0.5603925

LagPar_PreOvi_P = [a_AF2_P, b_AF2_P]
