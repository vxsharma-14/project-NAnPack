# coding: utf-8
import numpy as np

#*********************************************************************************
def BC(U, iMax, jMax, BCType):
    '''Assigns boundary condition at the walls, inlet and outlet'''

    if iMax % 2 == 0:
        iMid = int(iMax/2)
        jMid = int(jMax/2)
    else:
        iMid = int((iMax - 1)/2)
        jMid = int((jMax - 1)/2)

    i1 = iMid
    i2 = iMid+1

    j1 = jMid+5
    j2 = jMid+6
    
    U_0 = 0.0
    U_100 = 100.0

    if BCType == 'Dirichlet': # Dirichlet boundary conditions
        # along j = 1
        U[0:-1, 0] = 40.0

        # along j = jMax
        U[0:-1, -1] = 10.0

        # along i = 1
        U[0, 0:-1] = 0.0

        # along i = iMax
        U[-1, 0:j1] = 0.0

    elif BCType == 'Neumann': # Neumann boundary conditions

        # along j = 1
        U[0:i1, 0] = U_0
        U[i2-1:-1, 0] = U_100

        # along j = jMax
        U[0:-1, -1] = U_0

        # along i = 1
        U[0, 0:-1] = U_0

        # along i = iMax
        U[-1, 0:j1] = U_100
        U[-1, j2-1:-1] = U_0

    elif BCType == 'Mixed': # Mixed type boundary conditions

        # along j = 1
        U[0:i1, 0] = U_0
        U[i2-1:-1, 0] = U_100

        # along j = jMax
        U[0:-1, -1] = U_0

        # along i = 1
        U[0, 0:-1] = U_0

        # along i = iMax
        U[-1, 0:j1] = U_100
        U[-1, j2-1:-1] = U_0

    else:
        print('Incorrect boundary condition type specified')
        print('Try "Dirichlet", "Neumann" or "Mixed" instead')
        print('Program will quit now')
        sys.exit()
    
    return U

#*********************************************************************************
def Initial(iMax, jMax):
    '''Assigns initial condition within the domain'''

    U = np.zeros((iMax,jMax), dtype='float64')

    return U

#*********************************************************************************

#*********************************************************************************
if __name__ == '__main__':
    InputParam(InFileName)
