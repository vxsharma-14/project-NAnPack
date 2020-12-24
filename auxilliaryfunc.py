# coding: utf-8
#*******************************************************************
def siFunc(Param, Ep):
    '''Calculate Equation 6-127 in Hoffmann Vol. 1
    '''
    y = Param
    yAbs = abs(y)
    if yAbs >= Ep:
        si = y
    else:
        si = (y**2 + Ep**2)/(2*Ep)
        
    return si

#*******************************************************************
def FourthOrderDamping(U, DampCoeff,iMax):
    '''Calculate the fourth-order damping term
       Equation 6-77 in Hoffmann Vol. 1
    '''
    D = U.copy() # initialize D
    D[2:-2] = -DampCoeff*(U[0:-4] - 4.0*U[1:-3] + 6.0*U[2:-2]\
                          - 4.0*U[3:-1] + U[4:])
    '''
    for i in range (2,iMax-2):
        D[i] = -DampCoeff*(U[i-2] - 4*U[i-1] + 6*U[i]\
                           - 4*U[i+1] + U[i+2])
    '''
    return D

#*******************************************************************
def SecondOrderDamping(U, DampCoeff):
    '''Calculate the second-order damping term
       Equation 6-80 in Hoffmann Vol. 1
    '''
    D = U.copy() # initialize D
    D[1:-1] = -DampCoeff*(U[0:-2] - 2.0*U[1:-1] + U[2:])
    
    return D

#*******************************************************************
