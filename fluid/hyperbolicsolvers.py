# coding: utf-8
import numpy as np
#********************************************************************************
def FirstOrderUpwind1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit first upwind differencing method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    Uold = U.copy()
    
    if Linear == 'Yes':

        U[1:-1] = U[1:-1]\
                  - 0.5*convX*(1 + np.sign(convX))*(U[1:-1] - U[0:-2])\
                  - 0.5*convX*(1 - np.sign(convX))*(U[2:] - U[1:-1])

    elif Linear == 'No':
        E = U*U/2
        U[1:-1] = U[1:-1] - convX*(E[1:-1] - E[0:-2])
    
    else:
        print('Use either Yes or No for argument "Linear"')
    
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def Lax1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit Lax method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    Uold = U.copy()
    if Linear == 'Yes':
        U[1:-1] = 0.5*(U[2:] + U[0:-2]) - 0.5*convX*(U[2:] - U[0:-2])
    
    elif Linear == 'No':
        U[1:-1] = 0.5*(U[2:] + U[0:-2]) - 0.25*(convX*\
                                                (U[2:]**2 - U[0:-2]**2))
    
    else:
        print('Use either Yes or No for argument "Linear"')

    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def MidpointLeapfrog1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit Midpoint Leapfrog method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
    if Linear == 'Yes':
        print('try this later')
    
    elif Linear == 'No':
        print('try this later')
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    return U

#********************************************************************************
def LaxWendroff1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit Lax-Wendroff method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    Uold = U.copy()
    if Linear == 'Yes':
        convX2 = convX**2
        U[1:-1] = U[1:-1] - 0.5*convX*(U[2:] - U[0:-2]) +\
                  0.5*convX2*(U[2:] - 2.0*U[1:-1] + U[0:-2])
    
    elif Linear == 'No':
        E = U*U/2
        convX2 = convX**2
        U[1:-1] = U[1:-1] - 0.5*convX*(E[2:] - E[0:-2]) +\
                  0.25*convX2*((U[2:] + U[1:-1])*(E[2:] - E[1:-1])\
                               - (U[1:-1] + U[0:-2])*(E[1:-1] - E[0:-2]))
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def LaxWendroffMultiStep1D(U, Uhalf, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit multi-step Lax-Wendroff method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    from globalmod import iMax
    Uold = U.copy()
    if Linear == 'Yes':
        for i in range (1,iMax-1):
            Uhalf[i] = 0.5*(U[i+1] + U[i]) - 0.5*convX*(U[i+1] - U[i])
            U[i] = U[i] - convX*(Uhalf[i] - Uhalf[i-1])
    
    elif Linear == 'No':
        print('try this later')
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def MacCormack1D(U, Utemp, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit multi-step Lax-Wendroff method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    from globalmod import iMax
    Uold = U.copy()
    if Linear == 'Yes':
        for i in range (1,iMax-1):
            Utemp[i] = U[i] - convX*(U[i+1] - U[i]) # Predictor step
            U[i] = 0.5*((U[i] + Utemp[i]) - convX*\
                        (Utemp[i] - Utemp[i-1])) # Corrector step
    
    elif Linear == 'No':
        E = U*U/2
        Etemp=Utemp*Utemp/2
        for i in range (1,iMax-1):
            Utemp[i] = U[i] - convX*(E[i+1] - E[i]) # Predictor step
            Etemp[i] = Utemp[i]*Utemp[i]/2
            U[i] = 0.5*((U[i] + Utemp[i]) - convX*\
                        (Etemp[i] - Etemp[i-1])) # Corrector step
    
    else:
        print('Use either Yes or No for argument "Linear"')

    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def ModifiedRungeKutta1D(U, convX, Linear='Yes', TVD='No'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit four-stage Modified Runge-Kutta method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''

    Uold = U.copy()
    
    if Linear == 'Yes':
        # 1st stage
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/8.0
        # -- update BC here
        # 2nd stage
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/6.0
        # -- update BC here
        # 3rd stage
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/4.0
        # -- update BC here
        # 4th stage
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/2.0
        # -- update BC here
    elif Linear == 'No':
        # 1st stage
        E = U*U/2
        U[1:-1] = Uold[1:-1] - convX*(E[2:] - E[0:-2])/8.0
        # -- update BC here
        # 2nd stage
        E = U*U/2
        U[1:-1] = Uold[1:-1] - convX*(E[2:] - E[0:-2])/6.0
        # -- update BC here
        # 3rd stage
        E = U*U/2
        U[1:-1] = Uold[1:-1] - convX*(E[2:] - E[0:-2])/4.0
        # -- update BC here
        # 4th stage
        E = U*U/2
        U[1:-1] = Uold[1:-1] - convX*(E[2:] - E[0:-2])/2.0
        # -- update BC here
        
    else:
        print('Use either Yes or No for argument "Linear"')
        
    if TVD == 'Yes':
        phiPlus = HartenYeeTVD()
        phiMinus = HartenYeeTVD()
        U = U - 0.5*convX*(phiPlus - phiMinus)

    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def EulersBTCS1D(n, U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the implicit Euler's Backward Time Central Space (BTCS) method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    from globalmod import iMax
    import tridiagonal as td
    Uold = U.copy()
    A = U.copy()
    B = U.copy()
    C = U.copy()
    D = U.copy()

    if Linear == 'Yes':
        A[1:-1] = 0.5*convX
        B[1:-1] = -1.0
        C[1:-1] = -0.5*convX
        D[1:-1] = -U[1:-1]

        U = td.TridiagonalSolver(A, B, C, D, U)

    elif Linear == 'No':
        print('try again later')
        
    else:
        print('Use either Yes or No for argument "Linear"')
    
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def CrankNicolson1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the implicit Crank-Nicolson method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    from globalmod import iMax
    import tridiagonal as td
    Uold = U.copy()
    A = U.copy()
    B = U.copy()
    C = U.copy()
    D = U.copy()
        
    if Linear == 'Yes':
        A[1:-1] = 0.25*convX
        B[1:-1] = -1.0
        C[1:-1] = -0.25*convX
        D[1:-1] = -U[1:-1] + 0.25*convX*(U[2:] - U[0:-2])

        U = td.TridiagonalSolver(A, B, C, D, U)

    elif Linear == 'No':
        print('try again later')
        
    else:
        print('Use either Yes or No for argument "Linear"')
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#********************************************************************************
def BeamAndWarming1D(U, convX, Linear='No'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the implicit Beam and Warming method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    from globalmod import iMax
    import tridiagonal as td
    Uold = U.copy()
    A = U.copy()
    B = U.copy()
    C = U.copy()
    D = U.copy()
    
    if Linear == 'Yes':
        print('try again later')

    elif Linear == 'No':
        E = U*U/2
        AA = U
        A[1:-1] = -0.25*convX*AA[0:-2]
        B[1:-1] = 1.0
        C[1:-1] = 0.25*convX*AA[2:]
        D[1:-1] = U[1:-1] - 0.5*convX*(E[2:] - E[0:-2]) +\
                  0.25*convX*(AA[2:]*U[2:] - AA[0:-2]*U[0:-2]) 

        U = td.TridiagonalSolver(A, B, C, D, U)

    else:
        print('Use either Yes or No for argument "Linear"')
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error
    
    return U, Error

#*******************************************************************************
def SecondOrderTVD(U, convX, LimiterFunc, Limiter):
    '''Solve a first-order non-linear 1D hyperbolic partial differential
       equation using the second-order TVD schemes and their various
       Limiter Functions and Limiters.
    '''
    from globalmod import iMax, Eps
    import tvdfunctions as tvd
    from auxilliaryfunc import siFunc
    import sys
    Uold = U.copy()
            
    for i in range (2,iMax-2):

        dUiPlus12 = Uold[i+1] - Uold[i]
        dUiMinus12 = Uold[i] - Uold[i-1]
        dUiPlus12Abs = abs(dUiPlus12)
        dUiMinus12Abs = abs(dUiMinus12)

        dUiPlus32 = Uold[i+2] - Uold[i+1]
        dUiMinus32 = Uold[i-1] - Uold[i-2]
        dUiPlus32Abs = abs(dUiPlus32)
        dUiMinus32Abs = abs(dUiMinus32)

        EiPlus1 = 0.5*Uold[i+1]*Uold[i+1]
        Ei = 0.5*Uold[i]*Uold[i]
        EiMinus1 = 0.5*Uold[i-1]*Uold[i-1]
        EiPlus2 = 0.5*Uold[i+2]*Uold[i+2]
        EiMinus2 = 0.5*Uold[i-2]*Uold[i-2]

        # Equation 6-128
        if dUiPlus12 != 0:
            alphaiPlus12 = (EiPlus1 - Ei)/dUiPlus12
        else:
            alphaiPlus12 = 0.5*(Uold[i+1] + Uold[i])

        if dUiMinus12 != 0:
            alphaiMinus12 = (Ei - EiMinus1)/dUiMinus12
        else:
            alphaiMinus12 = 0.5*(Uold[i] + Uold[i-1])
                
        # .........................................................................
        if LimiterFunc == 'Harten-Yee-Upwind':

            # .....................................................................
            # Calculate some variables required in Harten-Yee Upwind.
            # .....................................................................
            if dUiPlus32 != 0:
                alphaiPlus32 = (EiPlus2 - EiPlus1)/dUiPlus32
            else:
                alphaiPlus32 = 0.5*(Uold[i+2] + Uold[i+1])

            if dUiMinus32 != 0:
                alphaiMinus32 = (EiMinus1 - EiMinus2)/dUiMinus32
            else:
                alphaiMinus32 = 0.5*(Uold[i-1] + Uold[i-2])
                
            # .....................................................................            
            if Limiter == 'G':            
                # Equation 6-130
                Gi = tvd.LimiterGforHYU(convX,alphaiPlus12,alphaiMinus12,\
                                        dUiPlus12,dUiMinus12,Eps)
                GiPlus1 = tvd.LimiterGforHYU(convX,alphaiPlus32,alphaiPlus12,\
                                             dUiPlus32,dUiPlus12,Eps)
                GiMinus1 = tvd.LimiterGforHYU(convX,alphaiMinus12,alphaiMinus32,\
                                              dUiMinus12,dUiMinus32,Eps)
            # .....................................................................
            else:
                print('****** Incorrect TVD Limiter for Harten-Yee Upwind. ******')
                print("Use Limiter = 'G'.")
                sys.exit('Program terminating now.')
                
            # Calculate Equation 6-126
            phiPlus, phiMinus = tvd.HartenYeeUp(dUiPlus12, dUiMinus12, GiPlus1,\
                                                Gi, GiMinus1, alphaiPlus12,\
                                                alphaiMinus12, Eps)
            
        # .........................................................................
        elif LimiterFunc == 'Modified-Harten-Yee-Upwind':
            
            # .....................................................................
            # Calculate some variables required in Modified Harten-Yee Upwind.
            # .....................................................................
            if dUiPlus32 != 0:
                alphaiPlus32 = (EiPlus2 - EiPlus1)/dUiPlus32
            else:
                alphaiPlus32 = 0.5*(Uold[i+2] + Uold[i+1])

            if dUiMinus32 != 0:
                alphaiMinus32 = (EiMinus1 - EiMinus2)/dUiMinus32
            else:
                alphaiMinus32 = 0.5*(Uold[i-1] + Uold[i-2])
                
            # .....................................................................
            if Limiter == 'G1':
                # Equation 6-132
                Gi = tvd.LimiterG1forHYU(dUiPlus12,dUiMinus12)
                GiPlus1 = tvd.LimiterG1forHYU(dUiPlus32,dUiPlus12)
                GiMinus1 = tvd.LimiterG1forHYU(dUiMinus12,dUiMinus32)            
            # .....................................................................
            elif Limiter == 'G2':
                # Equation 6-133
                Gi = tvd.LimiterG2forHYU(dUiPlus12,dUiMinus12)
                GiPlus1 = tvd.LimiterG2forHYU(dUiPlus32,dUiPlus12)
                GiMinus1 = tvd.LimiterG2forHYU(dUiMinus12,dUiMinus32)
            # .....................................................................
            elif Limiter == 'G3':
                # Equation 6-134
                Gi = tvd.LimiterG3forHYU(dUiPlus12,dUiMinus12)
                GiPlus1 = tvd.LimiterG3forHYU(dUiPlus32,dUiPlus12)
                GiMinus1 = tvd.LimiterG3forHYU(dUiMinus12,dUiMinus32)
            # .....................................................................
            elif Limiter == 'G4':
                # Equation 6-135
                Gi = tvd.LimiterG4forHYU(dUiPlus12,dUiMinus12)
                GiPlus1 = tvd.LimiterG4forHYU(dUiPlus32,dUiPlus12)
                GiMinus1 = tvd.LimiterG4forHYU(dUiMinus12,dUiMinus32)
            # .....................................................................
            elif Limiter == 'G5':
                # Equation 6-136
                Gi = tvd.LimiterG5forHYU(dUiPlus12,dUiMinus12)
                GiPlus1 = tvd.LimiterG5forHYU(dUiPlus32,dUiPlus12)
                GiMinus1 = tvd.LimiterG5forHYU(dUiMinus12,dUiMinus32)
            # .....................................................................
            else:
                print('****** Incorrect TVD Limiter for Modified Harten-Yee Upwind. ******')
                print("Options 'G1', 'G2', 'G3', 'G4', 'G5'.")
                sys.exit('Program terminating now.')
                
            # Calculate Equation 6-131                
            phiPlus, phiMinus = tvd.ModHartenYeeUp(dUiPlus12, dUiMinus12, GiPlus1,\
                                                   Gi, GiMinus1, alphaiPlus12,\
                                                   alphaiMinus12, convX, Eps)
            
        # .........................................................................
        elif LimiterFunc == 'Roe-Sweby-Upwind':
            
            # .....................................................................
            # Calculate some variables required in Row-Sweby Upwind.
            # .....................................................................
            zero_filter = 1.e-7 # variable to filter out division by zero
            if alphaiPlus12 != 0:
                sigPlus12 = alphaiPlus12/abs(alphaiPlus12)
            else:
                sigPlus12 = 0.0
                
            if alphaiMinus12 != 0:
                sigMinus12 = alphaiMinus12/abs(alphaiMinus12)
            else:
                sigMinus12 = 0.0
            
            if dUiPlus12 >= zero_filter:
                riPlus = (Uold[i+1+int(sigPlus12)] - Uold[i+int(sigPlus12)])/\
                         dUiPlus12
            else:
                riPlus = 0.0
                
            if dUiMinus12 >= zero_filter:
                riMinus = (Uold[i+int(sigMinus12)] - Uold[i+1+int(sigMinus12)])/\
                          dUiMinus12               
            else:
                riMinus = 0.0
            
            # .....................................................................
            if Limiter == 'G1':
                # Equation 6-138
                Gi = tvd.LimiterG1forRSU(riPlus)
                GiMinus1 = tvd.LimiterG1forRSU(riMinus)         
            # .....................................................................
            elif Limiter == 'G2':
                # Equation 6-139
                Gi = tvd.LimiterG2forRSU(riPlus)
                GiMinus1 = tvd.LimiterG2forRSU(riMinus)      
            # .....................................................................
            elif Limiter == 'G3':
                # Equation 6-140
                Gi = tvd.LimiterG3forRSU(riPlus)
                GiMinus1 = tvd.LimiterG3forRSU(riMinus)
            # .....................................................................
            else:
                print('****** Incorrect TVD Limiter for Roe-Sweby Upwind. ******')
                print("Options: 'G1', 'G2', 'G3'.")
                sys.exit('Program terminating now.')
                
            # Calculate Equation 6-137
            phiPlus, phiMinus = tvd.RoeSwebyUp(dUiPlus12, dUiMinus12, Gi, GiMinus1,\
                                               alphaiPlus12, alphaiMinus12, convX)
        
        # .........................................................................
        elif LimiterFunc == 'Davis-Yee-Symmetric':
            
            # .....................................................................
            if Limiter == 'G1':
                # Equation 6-142
                GiPlus12 = tvd.LimiterG1forDYS(dUiMinus12, dUiPlus12, dUiPlus32)
                GiMinus12 = tvd.LimiterG1forDYS(dUiMinus32, dUiMinus12, dUiPlus12)
            # .....................................................................
            elif Limiter == 'G2':
                # Equation 6-143
                GiPlus12 = tvd.LimiterG2forDYS(dUiMinus12, dUiPlus12, dUiPlus32)
                GiMinus12 = tvd.LimiterG2forDYS(dUiMinus32, dUiMinus12, dUiPlus12)
            # .....................................................................
            elif Limiter == 'G3':
                # Equation 6-144
                GiPlus12 = tvd.LimiterG3forDYS(dUiPlus12, dUiMinus12, dUiPlus32)
                GiMinus12 = tvd.LimiterG3forDYS(dUiMinus12, dUiMinus32, dUiPlus12)
            # .....................................................................
            else:
                print('****** Incorrect TVD Limiter for Davis-Yee Symmetric. ******')
                print("Options: 'G1', 'G2', 'G3'.")
                sys.exit('Program terminating now.')
                
            # Calculate Equation 6-137
            phiPlus, phiMinus = tvd.DavisYeeSym(dUiPlus12, dUiMinus12, GiPlus12,\
                                                GiMinus12, alphaiPlus12, alphaiMinus12,\
                                                convX, Eps)
        
        # .........................................................................
        else:
            print('****** Incorrect Limiter Function Selection ******.')
            print("Options: 'Harten-Yee-Upwind'.")
            print("         'Modified-Harten-Yee-Upwind'.")
            print("         'Roe-Sweby-Upwind'.")
            print("         'Davis-Yee-Symmetric'.")
            sys.exit('Program terminating now.')                   

        # Equation 6-124 and 6-125 in Hoffmann Vol. 1
        hPlus = 0.5*(EiPlus1 + Ei + phiPlus)
        hMinus = 0.5*(Ei + EiMinus1 + phiMinus)

        # Equation 6-123    
        U[i] = Uold[i] - convX*(hPlus - hMinus)
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error

    return U, Error

#*******************************************************************************
