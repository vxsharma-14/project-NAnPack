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
        if n == 1:
            E = U.copy() # initialize E array

        E[:] = 0.5*U[:]**2
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
        convX2 = convX**2
        U[1:-1] = U[1:-1] - 0.5*convX*(E[2:] - E[0:-2]) +\
                  0.25*convX2*((U[2:] + U[1:-1])*(E[2:] - E[1:-1])\
                               - (U[1:-1] + U[0:-2])*(E[1:-1] + E[0:-2]))
    
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
        for i in range (1,iMax-1):
            Utemp[i] = U[i] - convX*(E[i+1] - E[i]) # Predictor step
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
        E[:] = [0.5*i*i for i in U]
        U[1:-1] = Uold[1:-1] - convX*(E[2:] - E[0:-2])/8.0
        # -- update BC here
        # 2nd stage
        E[:] = [0.5*i*i for i in U]
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/6.0
        # -- update BC here
        # 3rd stage
        E[:] = [0.5*i*i for i in U]
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/4.0
        # -- update BC here
        # 4th stage
        E[:] = [0.5*i*i for i in U]
        U[1:-1] = Uold[1:-1] - convX*(U[2:] - U[0:-2])/2.0
        # -- update BC here
        
    else:
        print('Use either Yes or No for argument "Linear"')
        
    if TVD == 'Yes':
        phiPlus = HartenYeeTVD()
        phiMinus = HartenYeeTVD()
        U[:] = U[:] - 0.5*convX*(phiPlus - phiMinus)

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
    if n == 1:
        A = np.zeros(iMax)
        B = np.zeros(iMax)
        C = np.zeros(iMax)
        D = np.zeros(iMax)
    
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
    if n == 1:
        A = np.zeros(iMax)
        B = np.zeros(iMax)
        C = np.zeros(iMax)
        D = np.zeros(iMax)
        E = np.zeros(iMax)
        AA = np.zeros(iMax)
    
    if Linear == 'Yes':
        print('try again later')

    elif Linear == 'No':
        E[:] = [0.5*i*i for i in U]
        AA[:] = U[:]
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

#********************************************************************************    
