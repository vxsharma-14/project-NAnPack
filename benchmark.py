# coding: utf-8
import numpy as np
import postprocess as post
import math

def HeatConduction(Length, Height, iMax, jMax, T1, T2, T3, T4):
    '''Obtain the analytical solution of the heat conduction within the 
       rectangular domain. Use the results to validate the numerical solution of
       - the parabolic 2D unsteady heat conduction PDE, or
       - the elliptic 2D steady state heat conduction PDE

       Case description: A rectangular bar is initially heated to a temperature T0
       (or initial conditions are guessed for steady state simulation).
       Subsequently, its surfaces are subject to the constant temperatures of T1, T2
       T3 and T4 as depicted below.
       
                            Y-axis
                            |
                            |          T3
                             _______________________
                            |                       |
                            |                       |
                            |                       |
                            |                       |
                        T2  |                       | T4
                            |                       |
                            |                       |
                            |                       |
                            |                       |
                            |_______________________|   _____ X-axis
                                       T1
    
    Call signature:

        HeatConduction(Length, Height, iMax, jMax)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, iMax, jMax, L, H, ConvCriteria, nMax,\
                   nWrite, OutFileName, nDisplay, HistFileName, FrameOpt,\
                   FrameWrite = InputSettings
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    Inf = 20
    
    L = Length
    W = Height
    
    X = np.zeros((iMax,jMax))
    Y = np.zeros((iMax,jMax))
    
    # Calculate grid steps in X and Y directions
    dX = L/(iMax - 1.0)
    dY = W/(jMax - 1.0)
    
    # Define X(i,j) and Y(i,j)
    for i in range (0,iMax):
        for j in range(0,jMax):
            X[i][j] = dX*i
            Y[i][j] = dY*j
            
    Ta = np.zeros((iMax,jMax))
    Tb = np.zeros((iMax,jMax))
    Tc = np.zeros((iMax,jMax))
    Td = np.zeros((iMax,jMax))
    
    Ta1 = np.zeros((iMax,jMax))
    Tb1 = np.zeros((iMax,jMax))
    Tc1 = np.zeros((iMax,jMax))
    Td1 = np.zeros((iMax,jMax))
    
    Ta2 = np.zeros((iMax,jMax))
    Tb2 = np.zeros((iMax,jMax))
    Tc2 = np.zeros((iMax,jMax))
    Td2 = np.zeros((iMax,jMax))
    
    Ta3 = np.zeros((iMax,jMax))
    Tb3 = np.zeros((iMax,jMax))
    Tc3 = np.zeros((iMax,jMax))
    Td3 = np.zeros((iMax,jMax))
    
    TA = np.zeros((iMax,jMax))
    TB = np.zeros((iMax,jMax))
    TC = np.zeros((iMax,jMax))
    TD = np.zeros((iMax,jMax))
    TAnaly = np.zeros((iMax,jMax))
    
    m = 0
    while m <= Inf:
        
        m = m + 1
        mPi = m*math.pi
        
        for i in range (1,iMax-1):
            for j in range (1,jMax-1):

                #**************************************************************
                # Compute Ta
                #**************************************************************
                Ta1[i][j] = (1.0 - math.cos(mPi))/mPi
                Ta2[i][j] = math.sinh(mPi*(W - Y[i][j])/L)/math.sinh(mPi*W/L)
                Ta3[i][j] = math.sin(mPi*X[i][j]/L)

                Ta[i][j] = Ta[i][j] + Ta1[i][j]*Ta2[i][j]*Ta3[i][j]

                #**************************************************************
                # Compute Tb
                #**************************************************************
                Tb1[i][j] = (1.0 - math.cos(mPi))/mPi
                Tb2[i][j] = math.sinh(mPi*Y[i][j]/L)/math.sinh(mPi*W/L)
                Tb3[i][j] = math.sin(mPi*X[i][j]/L)

                Tb[i][j] = Tb[i][j] + Tb1[i][j]*Tb2[i][j]*Tb3[i][j]

                #**************************************************************
                # Compute Tc
                #**************************************************************
                Tc1[i][j] = (1.0 - math.cos(mPi))/mPi
                Tc2[i][j] = math.sinh(mPi*(L - X[i][j])/W)/math.sinh(mPi*L/W)
                Tc3[i][j] = math.sin(mPi*Y[i][j]/W)

                Tc[i][j] = Tc[i][j] + Tc1[i][j]*Tc2[i][j]*Tc3[i][j]

                #**************************************************************
                # Compute Td
                #**************************************************************
                Td1[i][j] = (1.0 - math.cos(mPi))/mPi
                Td2[i][j] = math.sinh(mPi*X[i][j]/W)/math.sinh(mPi*L/W)
                Td3[i][j] = math.sin(mPi*Y[i][j]/W)

    
                Td[i][j] = Td[i][j] + Td1[i][j]*Td2[i][j]*Td3[i][j]
                
    TA[1:-1,1:-1] = T1*2.0*Ta[1:-1,1:-1]
    TB[1:-1,1:-1] = T3*2.0*Tb[1:-1,1:-1]
    TC[1:-1,1:-1] = T2*2.0*Tc[1:-1,1:-1]
    TD[1:-1,1:-1] = T4*2.0*Td[1:-1,1:-1]

    TAnaly[1:-1,1:-1] = TA[1:-1,1:-1] + TB[1:-1,1:-1] +\
                        TC[1:-1,1:-1] + TD[1:-1,1:-1]

    # At boundaries
    TAnaly[:,0] = T1  # At Y=0
    TAnaly[:,-1] = T3 # At Y=H
    TAnaly[0,:] = T2  # At X=0
    TAnaly[-1,:] = T4 # At X=L

    # At corners
    TAnaly[0,0] = T1   # At X=0, Y=0
    TAnaly[-1,0] = T1  # At X=L, Y=0
    TAnaly[0,-1] = T3  # At X=0, Y=H
    TAnaly[-1,-1] = T3 # At x=L, Y=H

    # Write to file
    
    post.WriteSolutionToFile('./Output/Analytic.dat', iMax, jMax, dX, dY, TAnaly )
  
    return TAnaly
