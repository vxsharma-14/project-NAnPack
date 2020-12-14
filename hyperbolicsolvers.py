# coding: utf-8
#********************************************************************************
def FirstOrderUpwind1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit first upwind differencing method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
    if Linear == 'Yes':
        if convX >= 0:
            U[1:-1] = U[1:-1] - convX*(U[1:-1] - U[0:-2])
        
        else:
            U[1:-1] = U[1:-1] - convX*(U[2:] - U[1:-1])
    
    elif Linear == 'No':
        U[1:-1] = U[1:-1] - convX*(E[1:-1] - E[0:-2])
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    return U

#********************************************************************************
def Lax1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit Lax method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
    if Linear == 'Yes':
        U[1:-1] = 0.5*(U[2:] + U[0:-2]) - 0.5*convX*(U[2:] - U[0:-2])
    
    elif Linear == 'No':
        U[1:-1] = 0.5*(U[2:] + U[0:-2]) - 0.25*(convX*\
                                                (U[2:]**2 - U[0:-2]**2))
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    return U

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
        
    return U

#********************************************************************************
def LaxWendroffMultiStep1D(U, Uhalf, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit multi-step Lax-Wendroff method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
    if Linear == 'Yes':
        for i in range (1,iMax-1):
            Uhalf[i] = 0.5*(U[i+1] + U[i]) - 0.5*convX*(U[i+1] - U[i])
            U[i] = U[i] - convX*(Uhalf[i] - Uhalf[i-1])
    
    elif Linear == 'No':
        print('try this later')
    
    else:
        print('Use either Yes or No for argument "Linear"')
        
    return U

#********************************************************************************
def MacCormack1D(U, Utemp, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit multi-step Lax-Wendroff method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
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

    return U

#********************************************************************************
def ModifiedRungeKutta1D(U, convX, Linear='Yes'):
    '''Solve a first-order 1D hyperbolic partial differential equation
       using the explicit Modified Runge-Kutta method.
       
       A note of caution
       For linear problems:     convX = conv*dT/dX
       For non-linear problems: convX = dT/dX
    '''
    
    Uold = U.copy()
    
    if Linear == 'Yes':
        print('try again later')
    
    elif Linear == 'No':
        print('try again later')
        
    else:
        print('Use either Yes or No for argument "Linear"')

    return U

#********************************************************************************
    
    
