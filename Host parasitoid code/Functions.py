
# List of the functions that support the script RunmeDDE.py
#
# Created by Luca Rossini on 17 January 2025
# Last update 9 April 2025
#
# E-mail: luca.rossini@ulb.be


# Import list

from Parameters import *
import matplotlib.pyplot as plt


# Temperature reader

def TempFunction(time, TempArray):

    if time >= len(TempArray):
        Temp = TempArray[0]
    else:
        Temp = TempArray[time]

    return Temp
    

 # Delay function 1 - Host Eggs and N1

def DelayFun_EggN1_H(T, Parameters):

    np.array(Parameters)
    
    Delay = (Parameters[0] * T**2 + Parameters[1] * T +
                Parameters[2]) / (T + Parameters[3])

    if Delay <= 0:
        Delay = 1
            
    elif np.any(T <= 15): # If out of the lower threshold
        T = 15
        Delay = (Parameters[0] * T**2 + Parameters[1] * T +
                    Parameters[2]) / (T + Parameters[3])

    else:
        Delay = (Parameters[0] * T**2 + Parameters[1] * T +
                    Parameters[2]) / (T + Parameters[3])
    
    return Delay


# Delay function 2 - Host N2

def DelayFun_N2_H(T, Parameters):

    Delay = Parameters[0] * np.exp(Parameters[1] * T)

    if Delay <= 0:
        Delay = 1

    elif np.any(T <= 15): # If out of the lower threshold
        T = 15
        Delay = Parameters[0] * np.exp(Parameters[1] * T)

    else:
        Delay = Parameters[0] * np.exp(Parameters[1] * T)
    
    return Delay


# Delay function 3 - Host N3

def DelayFun_N3_H(T, Parameters):

    Delay = Parameters[0] * T**2 + Parameters[1] * T + Parameters[2]

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] * T**2 + Parameters[1] * T + Parameters[2];
    
    return Delay


# Delay function 4 - Host N4

def DelayFun_N4_H(T, Parameters):

    Delay = Parameters[0] * T**3 + Parameters[1] * T**2 + Parameters[2] * T + Parameters[3]

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] * T**3 + Parameters[1] * T**2 + Parameters[2] * T + Parameters[3]
    
    return Delay


# Delay function 5 - Host N5

def DelayFun_N5_H(T, Parameters):

    Delay = (Parameters[0] * T + Parameters[1]) / (T + Parameters[2])

    if Delay <= 0:
        Delay = 1
    else:
        Delay = (Parameters[0] * T + Parameters[1]) / (T + Parameters[2])
    
    return Delay


# Delay function 6 - Preoviposition time of the host

def DelayFun_PreOvi_H(T, Parameters):

    Delay = Parameters[0] * T**Parameters[1] + Parameters[2]

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] * T**Parameters[1] + Parameters[2]
    
    return Delay


# Delay function 7 - Parasitoid Eggs

def DelayFun_Egg_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    if Delay <= 0:
        Delay = 1

    else:
        Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    return Delay


# Delay function 8 - Parasitoid L1

def DelayFun_L1_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    if Delay <= 0:
        Delay = 1

    else:
        Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))
    
    return Delay


# Delay function 9 - Parasitoid L2

def DelayFun_L2_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))
    
    return Delay


# Delay function 10 - Parasitoid L3

def DelayFun_L3_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))
    
    return Delay


# Delay function 11 - Parasitoid Pupa

def DelayFun_P_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]
    c = Parameters[2]

    Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))

    if Delay <= 0:
        Delay = 1

    else:
        Delay = Parameters[0] + Parameters[1] * (np.exp(Parameters[2] * T))
    
    return Delay


# Delay function 12 - Parasitoid preoviposition time

def DelayFun_PreOvi_P(T, Parameters):

    a = Parameters[0]
    b = Parameters[1]

    Delay = Parameters[0] * T + Parameters[1]

    if Delay <= 0:
        Delay = 1
    else:
        Delay = Parameters[0] * T + Parameters[1]
    
    return Delay
    

# Temperature-dependent delays - To feed the DDE

def Delays(t, TempArray, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P):
    
        # Manage time and temperatures

    time = int(t)

    DayTemp = TempFunction(time, TempArray)

    # Host
    
        # Egg lag - Host
    tau_E_H = DelayFun_EggN1_H(DayTemp, LagPar_Egg_H)
    
        # Parasitised eggs lag - Host
    tau_EPar_H = 1
    
        # N1 lag - Host
    tau_N1_H = DelayFun_EggN1_H(DayTemp, LagPar_N1_H)
    
        # N2 lag - Host
    tau_N2_H = DelayFun_N2_H(DayTemp, LagPar_N2_H)
    
        # N3 lag - Host
    tau_N3_H = DelayFun_N3_H(DayTemp, LagPar_N3_H)
    
        # N4 lag - Host
    tau_N4_H = DelayFun_N4_H(DayTemp, LagPar_N4_H)
    
        # N5 lag - Host
    tau_N5_H = DelayFun_N5_H(DayTemp, LagPar_N5_H)
    
        # Male lag - Host
    tau_AM_H = 1
    
        # Female 1 lag - Host
    tau_AF1_H = 1
                               
        # Female 2 lag - Host
    tau_AF2_H = DelayFun_PreOvi_H(DayTemp, LagPar_PreOvi_H)
    
    # Parasitoid

        # Egg lag - Parasitoid
    tau_E_P = DelayFun_Egg_P(DayTemp, LagPar_Egg_P)

        # L1 lag  - Parasitoid
    tau_L1_P = DelayFun_L1_P(DayTemp, LagPar_L1_P)

        # L2 lag  - Parasitoid
    tau_L2_P = DelayFun_L2_P(DayTemp, LagPar_L2_P)

        # L3 lag  - Parasitoid
    tau_L3_P = DelayFun_L3_P(DayTemp, LagPar_L3_P)

        # P_lag  - Parasitoid
    tau_P_P = DelayFun_P_P(DayTemp, LagPar_P_P)

        # Male lag - Parasitoid
    tau_AM_P = 1

        # Female 1 lag - Parasitoid
    tau_AF1_P = 1

        # Female 2 lag - Parasitoid
    tau_AF2_P = DelayFun_PreOvi_P(DayTemp, LagPar_PreOvi_P)
    
        # Store the results
    Del = [tau_E_H, tau_EPar_H, tau_N1_H, tau_N2_H, tau_N3_H, tau_N4_H, tau_N5_H, tau_AM_H, tau_AF1_H, tau_AF2_H, tau_E_P, tau_L1_P, tau_L2_P, tau_L3_P, tau_P_P, tau_AM_P, tau_AF1_P, tau_AF2_P]
    
    return np.array(Del)


# Definition of Lactin-2 rate function

def LacTwo(Par, T):

    G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - (Par[1] - T) / Par[2]) + Par[3]

    if G > 0:
        G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - (Par[1] - T) / Par[2]) + Par[3]
    else:
        G = 0.001
    
    return G
    

# Definition of Logan rate function
    
def Log(Par, T):

    G = Par[0] * (np.exp(Par[1] * T) - np.exp(Par[1] * Par[2] - ((Par[2] - T) / Par[3])))
    
    if G > 0:
        G = Par[0] * (np.exp(Par[1] * T) - np.exp(Par[1] * Par[2] - ((Par[2] - T) / Par[3])))
    else:
        G = 0.001
    
    return G


# Definition of Sharpe and De Michele rate function

def SDM(Par, Temp):

    T = Temp # Temperature in Kelvin for this function!
    
    if T <= 0:
        T = 1
    else:
        T = Temp

    G = (T * (np.exp(Par[0] - (Par[1] / T)))) / ((1 + np.exp(Par[2] - (Par[3] / T)) + np.exp(Par[4] - (Par[5] / T))))

    if G > 0:
        G = (T * (np.exp(Par[0] - (Par[1] / T)))) / ((1 + np.exp(Par[2] - (Par[3] / T)) + np.exp(Par[4] - (Par[5] / T))))
    else:
        G = 0.001
    
    return G


# Definition of Lactin-1 rate function

def LacOne(Par, T):

    G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - ( (Par[1] - T) / Par[2] ) )

    if G > 0:
        G = np.exp(Par[0] * T) - np.exp(Par[0] * Par[1] - ( (Par[1] - T) / Par[2] ) )
    else:
        G = 0.001
        
    return G
    
    
# Definition of the Briere rate function

def BriRate(Par, T):

    a = Par[0]
    T_L = Par[1]
    T_M = Par[2]
    m = Par[3]

    G = a * T * (T - T_L) * (T_M - T)**(1/m)
            
    if G > 0:
    
            G = a * T * (T - T_L) * (T_M - T)**(1/m)
    else:
        G = 0.001

    return G
    

# Definition of the fertility rate function

def FertFunc(Par, T):

    FP = Par[0] * ( ((Par[1] + 1) / (np.pi * (Par[2] ** (2 * Par[1] + 2)))) * ((Par[2] ** 2) - ( ((T - Par[4]) ** 2) + (Par[3] ** 2)) ) ** Par[1] )

    if FP > 0:
        FP = Par[0] * ( ((Par[1] + 1) / (np.pi * (Par[2] ** (2 * Par[1] + 2)))) * ((Par[2] ** 2) - ( ((T - Par[4]) ** 2) + (Par[3] ** 2)) ) ** Par[1] )
    else:
        FP = 0

    return FP


# Definition of the Mortality rate function

    # Host

def MortFunc_H(Par, T):

    MP1 = Par[0] * (T ** 4)
    MP2 = Par[1] * (T ** 3)
    MP3 = Par[2] * (T ** 2)
    MP4 = Par[3] * T
    MP5 = Par[4]

    MP = MP1 + MP2 + MP3 + MP4 + MP5
    
    if MP >= 0 and MP <=1:

        MP = MP1 + MP2 + MP3 + MP4 + MP5
        
    elif MP >= 1:
    
        MP = 1
        
    else:
        MP = 0
    
    return MP / 3

    # Parasitoid

def MortFunc_P(Par, T):

    MP1 = Par[0] * (T ** 4)
    MP2 = Par[1] * (T ** 3)
    MP3 = Par[2] * (T ** 2)
    MP4 = Par[3] * T
    MP5 = Par[4]

    MP = MP1 + MP2 + MP3 + MP4 + MP5
    
    if MP >= 0 and MP <=1:

        MP = MP1 + MP2 + MP3 + MP4 + MP5
        
    elif MP >= 1:
    
        MP = 1
        
    else:
        MP = 0

    return 0.1#MP


# Compute the normalisation coefficients for the development rates – Including adult survival

def NormCoefficients(DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P):
                    
    # Array to store solutions

    DevRate_Egg_H = []
    DevRate_N1_H = []
    DevRate_N2_H = []
    DevRate_N3_H = []
    DevRate_N4_H = []
    DevRate_N5_H = []
    DevRate_AF_H = []
    DevRate_AM_H = []
    
    DevRate_Egg_P = []
    DevRate_L1_P = []
    DevRate_L2_P = []
    DevRate_L3_P = []
    DevRate_P_P = []
    DevRate_AF_P = []
    DevRate_AM_P = []
    
    # Temperature series to evaluate the functions in their domain of existance
    
    Temp = np.linspace(10, 50, 1000)
    
    # Compute development rates in their domain of existance
    
        # Host
        
    DevRate_Egg_H = np.array([BriRate(DevPar_Egg_H, T) for T in Temp])
    DevRate_N1_H = np.array([BriRate(DevPar_N1_H, T) for T in Temp])
    DevRate_N2_H = np.array([BriRate(DevPar_N2_H, T) for T in Temp])
    DevRate_N3_H = np.array([BriRate(DevPar_N3_H, T) for T in Temp])
    DevRate_N4_H = np.array([BriRate(DevPar_N4_H, T) for T in Temp])
    DevRate_N5_H = np.array([BriRate(DevPar_N5_H, T) for T in Temp])
    DevRate_AF_H = np.array([BriRate(DevPar_AF_H, T) for T in Temp])
    DevRate_AM_H = np.array([BriRate(DevPar_AM_H, T) for T in Temp])
    
        # Parasitoid
        
    DevRate_Egg_P = np.array([LacTwo(DevPar_Egg_P, T) for T in Temp])
    DevRate_L1_P = np.array([Log(DevPar_L1_P, T) for T in Temp])
    DevRate_L2_P = np.array([Log(DevPar_L2_P, T) for T in Temp])
    DevRate_L3_P = np.array([SDM(DevPar_L3_P, T) for T in Temp])
    DevRate_P_P = np.array([LacTwo(DevPar_P_P, T) for T in Temp])
    DevRate_AF_P = np.array([LacOne(DevPar_AF_P, T) for T in Temp])
    DevRate_AM_P = np.array([LacOne(DevPar_AM_P, T) for T in Temp])
    
    # Compute maximum values
    
        # Host

    Egg_Norm_H = np.max(DevRate_Egg_H)
    N1_Norm_H = np.max(DevRate_N1_H)
    N2_Norm_H = np.max(DevRate_N2_H)
    N3_Norm_H = np.max(DevRate_N3_H)
    N4_Norm_H = np.max(DevRate_N4_H)
    N5_Norm_H = np.max(DevRate_N5_H)
    AF_Norm_H = np.max(DevRate_AF_H)
    AM_Norm_H = np.max(DevRate_AM_H)
    
        # Parasitoid
        
    Egg_Norm_P = np.max(DevRate_Egg_P)
    L1_Norm_P = np.max(DevRate_L1_P)
    L2_Norm_P = np.max(DevRate_L2_P)
    L3_Norm_P = np.max(DevRate_L3_P)
    P_Norm_P = np.max(DevRate_P_P)
    AF_Norm_P = np.max(DevRate_AF_P)
    AM_Norm_P = np.max(DevRate_AM_P)
    
    return [Egg_Norm_H, N1_Norm_H, N2_Norm_H, N3_Norm_H, N4_Norm_H, N5_Norm_H, AM_Norm_H, AF_Norm_H, Egg_Norm_P, L1_Norm_P, L2_Norm_P, L3_Norm_P, P_Norm_P, AF_Norm_P, AM_Norm_P]
    
    
# Compute HOST Egg daily development rate normalised

def Egg_NormDevRate_H(Par, T, NormCoeff):

    Egg_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return Egg_NormDevRate_Value / 2
    

# Compute HOST N1 daily development rate normalised

def N1_NormDevRate_H(Par, T, NormCoeff):

    N1_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return N1_NormDevRate_Value
    
    
# Compute HOST N2 daily development rate normalised

def N2_NormDevRate_H(Par, T, NormCoeff):

    N2_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return N2_NormDevRate_Value
    

# Compute HOST N3 daily development rate normalised

def N3_NormDevRate_H(Par, T, NormCoeff):

    N3_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return N3_NormDevRate_Value
    
    
# Compute HOST N4 daily development rate normalised

def N4_NormDevRate_H(Par, T, NormCoeff):

    N4_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff

    return N4_NormDevRate_Value
    

# Compute HOST N5 daily development rate normalised

def N5_NormDevRate_H(Par, T, NormCoeff):

    N5_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return N5_NormDevRate_Value
    
    
# Compute HOST Adult male daily development rate normalised

def AM_NormDevRate_H(Par, T, NormCoeff):

    AM_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return AM_NormDevRate_Value
    
    
# Compute HOST Adult female daily development rate normalised

def AF_NormDevRate_H(Par, T, NormCoeff):

    AF_NormDevRate_Value = BriRate(Par, T) #/ NormCoeff
    
    return AF_NormDevRate_Value
    

# Compute PARASITOID Egg daily development rate normalised

def Egg_NormDevRate_P(Par, T, NormCoeff):

    Egg_NormDevRate_Value = LacTwo(Par, T) #/ NormCoeff
    
    return Egg_NormDevRate_Value / 20


# Compute PARASITOID L1 daily development rate normalised

def L1_NormDevRate_P(Par, T, NormCoeff):

    L1_NormDevRate_Value = Log(Par, T) #/ NormCoeff
    
    return L1_NormDevRate_Value / 15
    

# Compute PARASITOID L2 daily development rate normalised

def L2_NormDevRate_P(Par, T, NormCoeff):

    L2_NormDevRate_Value = Log(Par, T) #/ NormCoeff
    
    return L2_NormDevRate_Value / 15
    

# Compute PARASITOID L3 daily development rate normalised

def L3_NormDevRate_P(Par, T, NormCoeff):

    L3_NormDevRate_Value = SDM(Par, T) #/ NormCoeff
    
    return L3_NormDevRate_Value / 15
    
    
# Compute PARASITOID Pupa daily development rate normalised

def P_NormDevRate_P(Par, T, NormCoeff):

    PNormDevRate_Value = LacTwo(Par, T) #/ NormCoeff
    
    return PNormDevRate_Value / 15


# Compute PARASITOID Adult Male daily development rate normalised

def AM_NormDevRate_P(Par, T, NormCoeff):

    AM_NormDevRate_Value = LacOne(Par, T) #/ NormCoeff
    
    return AM_NormDevRate_Value / 15
    

# Compute PARASITOID Adult Female daily development rate normalised

def AF_NormDevRate_P(Par, T, NormCoeff):

    AF_NormDevRate_Value = LacOne(Par, T) #/ NormCoeff
    
    return AF_NormDevRate_Value / 15
    
    
    
# DDE system with development rates NORMALISED - To feed the DDE solver!

def DDE_FuncNormalised(Y, t, FertPar_H, MortPar_Egg_H, MortPar_N1_H, MortPar_N2_H, MortPar_N3_H, MortPar_N4_H, MortPar_N5_H, MortPar_Egg_P, MortPar_L1_P, MortPar_L2_P, MortPar_L3_P, DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P, NormCoeff_Egg_H, NormCoeff_N1_H, NormCoeff_N2_H, NormCoeff_N3_H, NormCoeff_N4_H, NormCoeff_N5_H, NormCoeff_AM_H, NormCoeff_AF_H, NormCoeff_Egg_P, NormCoeff_L1_P, NormCoeff_L2_P, NormCoeff_L3_P, NormCoeff_P_P, NormCoeff_AF_P, NormCoeff_AM_P, TempArray, SR_H, SR_P, a, T_h, b):
                   
    # Calculate the daily temperature from the 'TemperatureInput.xlsx' file in 'Parameters.py'

    temp = TempFunction(int(t), TempArray)
    
    Z = Delays(t, TempArray, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P)
    
    # Approximation of the lags - Needed by ddeint in Y(t - Z)
    
        # Host
        
    # Z[0] = Egg lag
    # Z[1] = Egg parasitised lag
    # Z[2] = N1 lag
    # Z[3] = N2 lag
    # Z[4] = N3 lag
    # Z[5] = N4 lag
    # Z[6] = N5 lag
    # Z[7] = AM lag
    # Z[8] = AF1 lag
    # Z[9] = AF2 lag
    
        # Parasitoid
    
    # Z[10] = Egg lag
    # Z[11] = L1 lag
    # Z[12] = L2 lag
    # Z[13] = L3 lag
    # Z[14] = P lag
    # Z[15] = AM lag
    # Z[16] = AF1 lag
    # Z[17] = AF2 lag


    # DDE system state variables

        # Host
        
    # Y(t)[0] = Egg without lag
    # Y(t)[1] = Egg parasitised without lag
    # Y(t)[2] = N1 without lag
    # Y(t)[3] = N2 without lag
    # Y(t)[4] = N3 without lag
    # Y(t)[5] = N4 without lag
    # Y(t)[6] = N5 without lag
    # Y(t)[7] = AM without lag
    # Y(t)[8] = AF1 without lag
    # Y(t)[9] = AF2 without lag

        # Parasitoid

    # Y(t)[10] = Egg without lag
    # Y(t)[11] = L1 without lag
    # Y(t)[12] = L2 without lag
    # Y(t)[13] = L3 without lag
    # Y(t)[14] = P without lag
    # Y(t)[15] = AM without lag
    # Y(t)[16] = AF1 without lag
    # Y(t)[17] = AF2 without lag


        # Define the functional responses to use into the model
        
        # Type 1
        
    # alpha = LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp) * Y(t)[0]
    
        # Type 2
        
    alpha = ((a * Y(t)[0]) / (1 + a * T_h * Y(t)[0])) * LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp)
    
        # Type 3
        
    # alpha = FertFunc(FertPar_P, temp) * ((a * Y(t)[0] ** 2) / (1 + b * Y(t)[0] **2))
    
        # Custom type (look at notes of 11/04/2025)
    
    # alpha = LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp) * 0.5
    
        # DDE Model to solve
    
    # Host - Halyomorpha halys
    
        # Egg

    dE_h = BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp) * Y(t - Z[9])[9] - Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * Y(t)[0] - MortFunc_H(MortPar_Egg_H, temp) * Y(t)[0] - alpha * Y(t)[0] * Y(t)[17]
    
        # Parasitised eggs

    dEPar_h = alpha * Y(t)[0] * Y(t)[17] - MortFunc_H(MortPar_Egg_P, temp) * Y(t)[10] - MortFunc_H(MortPar_L1_P, temp) * Y(t)[11] - MortFunc_H(MortPar_L2_P, temp) * Y(t)[12] - MortFunc_H(MortPar_L3_P, temp) * Y(t)[13] - P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * Y(t)[14]

        # N1

    dN1_h = Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * Y(t - Z[0])[0] - N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * Y(t)[2] - MortFunc_H(MortPar_N1_H, temp) * Y(t)[2]

        # N2

    dN2_h = N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * Y(t - Z[2])[2] - N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * Y(t)[3] -  MortFunc_H(MortPar_N2_H, temp) * Y(t)[3]

        # N3

    dN3_h = N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * Y(t - Z[3])[3] - N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * Y(t)[4] -  MortFunc_H(MortPar_N3_H, temp) * Y(t)[4]
    
        # N4

    dN4_h = N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * Y(t - Z[4])[4] - N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * Y(t)[5] -  MortFunc_H(MortPar_N4_H, temp) * Y(t)[5]
    
        # N5

    dN5_h = N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * Y(t - Z[5])[5] - N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * Y(t)[6] -  MortFunc_H(MortPar_N5_H, temp) * Y(t)[6]
    
        # Adult males

    dAM_h = (1 - SR_H) * N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * Y(t - Z[6])[6] - (1 - AM_NormDevRate_H(DevPar_AM_H, temp, NormCoeff_AM_H)) * Y(t)[7]

        # Female-1 (Non-mated)

    dAF1_h = SR_H * N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * Y(t - Z[6])[6] - Y(t)[8]
    
        # Female-2 (Mated)
    dAF2_h = Y(t)[8] - (1 - AF_NormDevRate_H(DevPar_AF_H, temp, NormCoeff_AF_H)) * Y(t)[9] - (1 - AF_NormDevRate_H(DevPar_AF_H, temp, NormCoeff_AF_H)) * Y(t)[9]
    
    # Parasitoid - Trissulcus japonicus
    
        # Egg

    dE_p = alpha * Y(t)[0] * Y(t - Z[17])[17] - Egg_NormDevRate_P(DevPar_Egg_P, temp, NormCoeff_Egg_P) * Y(t)[10] - MortFunc_P(MortPar_Egg_P, temp) * Y(t)[10]

        # L1
    
    dL1_p = Egg_NormDevRate_P(DevPar_Egg_P, temp, NormCoeff_Egg_P) * Y(t - Z[10])[10] - L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * Y(t)[11] - MortFunc_P(MortPar_L1_P, temp) * Y(t)[11]

        # L2
    
    dL2_p = L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * Y(t - Z[11])[11] - L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * Y(t)[12] - MortFunc_P(MortPar_L2_P, temp) * Y(t)[12]

        # L3
    
    dL3_p = L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * Y(t - Z[12])[12] - L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * Y(t)[13] - MortFunc_P(MortPar_L3_P, temp) * Y(t)[13]

        # P
    
    dP_p = L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * Y(t - Z[13])[13] - P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * Y(t)[14]

        # AM
    
    dAM_p = (1 - SR_P) * P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * Y(t - Z[14])[14] - (1 - AM_NormDevRate_P(DevPar_AM_P, temp, NormCoeff_AM_P)) * Y(t)[15]

        # Female-1 (Non mated)
    
    dAF1_p = SR_P * P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * Y(t - Z[14])[14] - Y(t)[16]

        # Female-2 (Mated)
    
    dAF2_p = Y(t)[16] - (1 - AF_NormDevRate_P(DevPar_AF_P, temp, NormCoeff_AF_P)) * Y(t)[16] - (1 - AF_NormDevRate_P(DevPar_AF_P, temp, NormCoeff_AF_P)) * Y(t)[17]
        
    
    PartialSolution = [dE_h, dEPar_h, dN1_h, dN2_h, dN3_h, dN4_h, dN5_h, dAM_h, dAF1_h, dAF2_h, dE_p, dL1_p, dL2_p, dL3_p, dP_p, dAM_p, dAF1_p,dAF2_p]
    
    return np.array(PartialSolution)



# DDE system COMPLETE with transition stages - To feed the DDE solver!

def DDE_FuncComplete(Y, t, FertPar_H, MortPar_Egg_H, MortPar_N1_H, MortPar_N2_H, MortPar_N3_H, MortPar_N4_H, MortPar_N5_H, MortPar_Egg_P, MortPar_L1_P, MortPar_L2_P, MortPar_L3_P, DevPar_Egg_H, DevPar_N1_H, DevPar_N2_H, DevPar_N3_H, DevPar_N4_H, DevPar_N5_H, DevPar_AM_H, DevPar_AF_H, DevPar_Egg_P, DevPar_L1_P, DevPar_L2_P, DevPar_L3_P, DevPar_P_P, DevPar_AF_P, DevPar_AM_P, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P, NormCoeff_Egg_H, NormCoeff_N1_H, NormCoeff_N2_H, NormCoeff_N3_H, NormCoeff_N4_H, NormCoeff_N5_H, NormCoeff_AM_H, NormCoeff_AF_H, NormCoeff_Egg_P, NormCoeff_L1_P, NormCoeff_L2_P, NormCoeff_L3_P, NormCoeff_P_P, NormCoeff_AF_P, NormCoeff_AM_P, TempArray, SR_H, SR_P, a, T_h, b):
                   
    # Calculate the daily temperature from the 'TemperatureInput.xlsx' file in 'Parameters.py'

    MortPupa_P = 0.001

    temp = TempFunction(int(t), TempArray)
    
    Z = Delays(t, TempArray, LagPar_Egg_H, LagPar_N1_H, LagPar_N2_H, LagPar_N3_H, LagPar_N4_H, LagPar_N5_H, LagPar_PreOvi_H, LagPar_Egg_P, LagPar_L1_P, LagPar_L2_P, LagPar_L3_P, LagPar_P_P, LagPar_PreOvi_P)
    assert all([zi > 0 for zi in Z]), f"Negative or zero delay found: {Z}"
    tau_Egg_H = Z[0]
    tau_N1_H = Z[1]
    tau_N2_H = Z[2]
    tau_N3_H = Z[3]
    tau_N4_H = Z[4]
    tau_N5_H = Z[5]
    tau_AM_H = Z[6]
    tau_AF1_H = Z[7]
    tau_AF2_H = Z[8]

    tau_Egg_P = Z[9]
    tau_L1_P = Z[10]
    tau_L2_P = Z[11]
    tau_L3_P = Z[12]
    tau_P_P = Z[13]
    tau_AM_P = Z[14]
    tau_AF1_P = Z[15]
    tau_AF2_P = Z[16]


    # DDE system state variables

        # Host
        
    Egg_H = Y(t)[0] # Egg without lag
    Eggtr_H = Y(t)[1] # Egg transient without lag
    N1_H = Y(t)[2] # N1 without lag
    N1tr_H = Y(t)[3] # N1 transient without lag
    N2_H = Y(t)[4] # N2 without lag
    N2tr_H = Y(t)[5] # N2 transient without lag
    N3_H = Y(t)[6] # N3 without lag
    N3tr_H = Y(t)[7] # N3 transient without lag
    N4_H = Y(t)[8] # N4 without lag
    N4tr_H = Y(t)[9] # N4 transient without lag
    N5_H = Y(t)[10] # N5 without lag
    N5tr_H = Y(t)[11] # N5 transient without lag
    AM_H = Y(t)[12] # Adult males without lag
    AF1_H = Y(t)[13] # Adult non mated females without lag
    AF2_H = Y(t)[14] # Adult mated females without lag

        # Parasitoid

    Egg_P = Y(t)[15] # Egg without lag
    Eggtr_P = Y(t)[16] # Egg transient without lag
    L1_P = Y(t)[17] # L1 without lag
    L1tr_P = Y(t)[18] # L1 transient without lag
    L2_P = Y(t)[19] # L2 without lag
    L2tr_P = Y(t)[20] # L2 transient without lag
    L3_P = Y(t)[21] # L3 without lag
    L3tr_P = Y(t)[22] # L3 transient without lag
    P_P = Y(t)[23] # Pupa without lag
    Ptr_P = Y(t)[24] # Pupa transient without lag
    AM_P = Y(t)[25] # Adult males without lag
    AF1_P = Y(t)[26] # Adult non mated females without lag
    AF2_P = Y(t)[27] # Adult mated females without lag


# Approximation of the lags - Needed by ddeint in Y(t - Z)
    
        # Host
    
    Egg_HLag = Y(t - Z[0])[1] # Egg lag
    N1_HLag = Y(t - Z[1])[3] # N1 lag
    N2_HLag = Y(t - Z[2])[5] # N2 lag
    N3_HLag = Y(t - Z[3])[7] # N3 lag
    N4_HLag = Y(t - Z[4])[9] # N4 lag
    N5_HLag = Y(t - Z[5])[11] # N5 lag
    AM_HLag = Y(t - Z[6])[12] # AM lag
    AF1_HLag = Y(t - Z[7])[13] # AF1 lag
    AF2_HLag = Y(t - Z[8])[14] # AF2 lag
    
        # Parasitoid
    
    Egg_PLag = Y(t - Z[9])[16] # Egg lag
    L1_PLag = Y(t - Z[10])[18] # L1 lag
    L2_PLag = Y(t - Z[11])[20] # L2 lag
    L3_PLag = Y(t - Z[12])[22] # L3 lag
    P_PLag = Y(t - Z[13])[24] # P lag
    AM_PLag = Y(t - Z[14])[25] # AM lag
    AF1_PLag = Y(t - Z[15])[26] # AF1 lag
    AF2_PLag = Y(t - Z[16])[27] # AF2 lag


        # Define the functional responses to use into the model
        
        # Type 1
        
    # alpha_1 = a * LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp) #* Egg_H
    # alpha_2 = a * LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp) #* Eggtr_H
    #print(alpha_1, alpha_2, BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp))

        # Type 2
        
    alpha_1 = (a / (1 + a * T_h * Egg_H)) #* LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp)
    #if alpha_1 < 10:
    #    alpha_1 = (a * Egg_H / (1 + a * T_h * Egg_H)) #* LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp)
    #else:
    #    alpha_1 = 0.001
    
    alpha_2 = (a / (1 + a * T_h * Eggtr_H)) #* LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp)
    #if alpha_2 < 10:
    #    alpha_2 = (a * Eggtr_H / (1 + a * T_h * Eggtr_H)) #* LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp)
    #else:
    #    alpha_2 = 0.001
    
    #print(alpha_1, alpha_2, BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp), LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp))
    
        # Type 3
        
    # alpha_1 = FertFunc(FertPar_P, temp) * ((a * Y(t)[0] ** 2) / (1 + b * Y(t)[0] **2))
    # alpha_2 = FertFunc(FertPar_P, temp) * ((a * Y(t)[1] ** 2) / (1 + b * Y(t)[1] **2))
    # print(alpha_1, alpha_2, BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp) * AF2_H)
    
        # Custom type (look at notes of 11/04/2025)
    
    # alpha = LacOne(DevPar_AF_P, temp) * FertFunc(FertPar_P, temp) * 0.5
    
        # DDE Model to solve
    
    # Host - Halyomorpha halys
    
        # Egg

    dE_h = BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp) * AF2_H - Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * Egg_H - MortFunc_H(MortPar_Egg_H, temp) * Egg_H - alpha_1 * Egg_H * AF2_P
    #print(BriRate(DevPar_AF_H, temp) * FertFunc(FertPar_H, temp), LacOne(DevPar_AF_P, temp)*FertFunc(FertPar_P, temp))
        # Egg transient

    dEtr_h = Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * Egg_H - Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * np.exp(-tau_Egg_H * MortFunc_H(MortPar_Egg_H, temp)) * Egg_HLag - MortFunc_H(MortPar_Egg_H, temp) * Eggtr_H - alpha_2 * Eggtr_H * AF2_P
    
        # N1

    dN1_h = Egg_NormDevRate_H(DevPar_Egg_H, temp, NormCoeff_Egg_H) * np.exp(-tau_Egg_H * MortFunc_H(MortPar_Egg_H, temp)) * Egg_HLag - N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * N1_H - MortFunc_H(MortPar_N1_H, temp) * N1_H
    
        # N1 transient

    dN1tr_h = N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * N1_H - N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * np.exp(-tau_N1_H * MortFunc_H(MortPar_N1_H, temp)) * N1_HLag - MortFunc_H(MortPar_N1_H, temp) * N1tr_H
    
        # N2

    dN2_h = N1_NormDevRate_H(DevPar_N1_H, temp, NormCoeff_N1_H) * np.exp(-tau_N1_H * MortFunc_H(MortPar_N1_H, temp)) * N1_HLag - N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * N2_H - MortFunc_H(MortPar_N2_H, temp) * N2_H
    
        # N2 transient

    dN2tr_h = N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * N2_H - N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * np.exp(-tau_N2_H * MortFunc_H(MortPar_N2_H, temp)) * N2_HLag - MortFunc_H(MortPar_N2_H, temp) * N2tr_H
    
        # N3

    dN3_h = N2_NormDevRate_H(DevPar_N2_H, temp, NormCoeff_N2_H) * np.exp(-tau_N2_H * MortFunc_H(MortPar_N2_H, temp)) * N2_HLag - N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * N3_H - MortFunc_H(MortPar_N3_H, temp) * N3_H
    
        # N3 transient

    dN3tr_h = N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * N3_H - N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * np.exp(-tau_N3_H * MortFunc_H(MortPar_N3_H, temp)) * N3_HLag - MortFunc_H(MortPar_N3_H, temp) * N3tr_H
     
        # N4

    dN4_h = N3_NormDevRate_H(DevPar_N3_H, temp, NormCoeff_N3_H) * np.exp(-tau_N3_H * MortFunc_H(MortPar_N3_H, temp)) * N3_HLag - N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * N4_H - MortFunc_H(MortPar_N4_H, temp) * N4_H
    
        # N4 transient

    dN4tr_h = N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * N4_H - N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * np.exp(-tau_N4_H * MortFunc_H(MortPar_N4_H, temp)) * N4_HLag - MortFunc_H(MortPar_N4_H, temp) * N4tr_H
    
        # N5

    dN5_h = N4_NormDevRate_H(DevPar_N4_H, temp, NormCoeff_N4_H) * np.exp(-tau_N4_H * MortFunc_H(MortPar_N4_H, temp)) * N4_HLag - N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * N5_H - MortFunc_H(MortPar_N5_H, temp) * N5_H
    
        # N5 transient

    dN5tr_h = N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * N5_H - N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * np.exp(-tau_N5_H * MortFunc_H(MortPar_N5_H, temp)) * N5_HLag - MortFunc_H(MortPar_N5_H, temp) * N5tr_H
    
        # Adult males

    dAM_h = (1 - SR_H) * N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * np.exp(-tau_N5_H * MortFunc_H(MortPar_N5_H, temp)) * N5_HLag - AM_NormDevRate_H(DevPar_AM_H, temp, NormCoeff_AM_H) * AM_H

        # Female-1 (Non-mated)

    dAF1_h = SR_H * N5_NormDevRate_H(DevPar_N5_H, temp, NormCoeff_N5_H) * np.exp(-tau_N5_H * MortFunc_H(MortPar_N5_H, temp)) * N5_HLag - AF1_H
    
        # Female-2 (Mated)

    dAF2_h = AF1_H - AF_NormDevRate_H(DevPar_AF_H, temp, NormCoeff_AF_H) * AF1_H - (AF_NormDevRate_H(DevPar_AF_H, temp, NormCoeff_AF_H)) * AF2_H
    
    # Parasitoid - Trissulcus japonicus
    
        # Egg

    dE_p = alpha_1 * Egg_H * AF2_P + alpha_2 * Eggtr_H * AF2_P - Egg_NormDevRate_P(DevPar_Egg_H, temp, NormCoeff_Egg_P) * Egg_P - MortFunc_P(MortPar_Egg_H, temp) * Egg_P 

        # Egg transient

    dEtr_p = Egg_NormDevRate_P(DevPar_Egg_P, temp, NormCoeff_Egg_P) * Egg_P - Egg_NormDevRate_P(DevPar_Egg_P, temp, NormCoeff_Egg_P) * np.exp(-tau_Egg_P * MortFunc_P(MortPar_Egg_P, temp)) * Egg_PLag - MortFunc_P(MortPar_Egg_P, temp) * Eggtr_P

        # L1

    dL1_p = Egg_NormDevRate_P(DevPar_Egg_P, temp, NormCoeff_Egg_P) * np.exp(-tau_Egg_P * MortFunc_P(MortPar_Egg_P, temp)) * Egg_PLag - L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * L1_P - MortFunc_P(MortPar_L1_P, temp) * L1_P

        # L1 transient

    dL1tr_p = L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * L1_P - L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * np.exp(-tau_L1_P * MortFunc_P(MortPar_L1_P, temp)) * L1_PLag - MortFunc_P(MortPar_L1_P, temp) * L1tr_P

        # L2

    dL2_p = L1_NormDevRate_P(DevPar_L1_P, temp, NormCoeff_L1_P) * np.exp(-tau_L1_P * MortFunc_P(MortPar_L1_P, temp)) * L1_PLag - L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * L2_P - MortFunc_P(MortPar_L2_P, temp) * L2_P

        # L2 transient

    dL2tr_p = L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * L2_P - L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * np.exp(-tau_L2_P * MortFunc_P(MortPar_L2_P, temp)) * L2_PLag - MortFunc_P(MortPar_L2_P, temp) * L2tr_P

        # L3

    dL3_p = L2_NormDevRate_P(DevPar_L2_P, temp, NormCoeff_L2_P) * np.exp(-tau_L2_P * MortFunc_P(MortPar_L2_P, temp)) * L2_PLag - L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * L3_P - MortFunc_P(MortPar_L3_P, temp) * L3_P

        # L3 transient

    dL3tr_p = L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * L3_P - L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * np.exp(-tau_L3_P * MortFunc_P(MortPar_L3_P, temp)) * L3_PLag - MortFunc_P(MortPar_L3_P, temp) * L3tr_P

        # P

    dP_p = L3_NormDevRate_P(DevPar_L3_P, temp, NormCoeff_L3_P) * np.exp(-tau_L3_P * MortFunc_P(MortPar_L3_P, temp)) * L3_PLag - P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * P_P - MortPupa_P * P_P

        # P transient

    dPtr_p = P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * P_P - P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * np.exp(-tau_P_P * MortPupa_P) * P_PLag - MortPupa_P * Ptr_P

        # AM
    
    dAM_p = (1 - SR_P) * P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * np.exp(-tau_P_P * MortPupa_P) * P_PLag - (1 - AM_NormDevRate_P(DevPar_AM_P, temp, NormCoeff_AM_P)) * AM_P

        # Female-1 (Non mated)
    
    dAF1_p = SR_P * P_NormDevRate_P(DevPar_P_P, temp, NormCoeff_P_P) * np.exp(-tau_P_P * MortPupa_P) * P_PLag - AF1_P

        # Female-2 (Mated)
    
    dAF2_p = AF1_P - AF_NormDevRate_P(DevPar_AF_P, temp, NormCoeff_AF_P) * AF1_P - AF_NormDevRate_P(DevPar_AF_P, temp, NormCoeff_AF_P) * AF2_P
        
    
    PartialSolution = [dE_h, dEtr_h, dN1_h, dN1tr_h, dN2_h, dN2tr_h, dN3_h, dN3tr_h, dN4_h, dN4tr_h, dN5_h, dN5tr_h, dAM_h, dAF1_h, dAF2_h, dE_p, dEtr_p, dL1_p, dL1tr_p, dL2_p, dL2tr_p, dL3_p, dL3tr_p, dP_p, dPtr_p, dAM_p, dAF1_p, dAF2_p]
    
    return np.array(PartialSolution)