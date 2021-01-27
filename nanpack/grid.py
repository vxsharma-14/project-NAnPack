#**************************************************************************
def ComputeGridPoints(Dimension, Length, delX, Height=None, delY=None):
    '''Returns the grid points along X and Y direction in the mesh.'''
    iMax = int(Length/delX) + 1
    if Dimension.upper() == '2D':
        jMax = int(Height/delY) + 1
    else:
        jMax = 0
    print('Calculating grid size: Completed.')

    return iMax, jMax

#**************************************************************************
def ComputeGridSteps(Dimension, Length, iMax, Height=None, jMax=None):
    '''Returns the uniform grid steps size along X and Y direction in the
    mesh.
    '''
    delX = Length/(iMax - 1)
    if Dimension.upper() == '2D':
        delY = Height/(jMax - 1)
    else:
        delY = 0.0
    print('Calculating grid step size: Completed.')

    return delX, delY

#**************************************************************************
def RectangularGrid(dX, iMax, dY=None, jMax=None):
    '''Returns X and Y grid point locations in a rectangular, uniform
    mesh in a cartesian coordinate system.
    '''
    import numpy as np

    if isinstance(dY, float) and isinstance(jMax, int):
        X = np.zeros((iMax,jMax), dtype='float')
        Y = np.zeros((iMax,jMax), dtype='float')
        for i in range(0,iMax):
            for j in range(0,jMax):
                X[i][j] = i*dX
                Y[i][j] = j*dY
    else:
        X = np.zeros((iMax), dtype='float')
        for i in range(0,iMax):
            X[i] = -9.0 + i*dX
        Y = 0.0
    print(f'Uniform rectangular grid generation in cartesian\
 coordinate system: Completed.')

    return X, Y

#**************************************************************************
def CurvilinearGrid(dX, iMax, dY=None, jMax=None):
    '''Returns X and Y grid point locations in a rectangular, uniform or a
    non-uniform mesh in a transformed coordinate system (Xi, Eta).
    '''
    print('Calculating X and Y locations of all grid points within\
 the mesh.')
    from .backend import gridmetrics
    from .backend import plotmetrics
    dXi = 1.0
    dEta = 1.0

    X, Y = RectangularGrid(dX, iMax, dY, jMax)
    dim = X.shape

    if len(dim) == 2: # Two dimensional
        Xi = [[i*dXi for j in range(0,jMax)] for i in range(0,iMax)]
        Eta = [[j*dEta for j in range (0,jMax)] for i in range (0,iMax)]
        XiX, XiY, EtaX, EtaY, JJ = gridmetrics.Metrics2D(X, Y)
        print('Grid metrics and Jacobian evaluation: Completed.')
        plotmetrics.PlotMetrics2D(X, Y, XiX, XiY, EtaX, EtaY)

    elif len(dim) == 1:
        Xi = [i*dX for i in range(0,iMax)]
        Eta = 0.0
        Xi, JJ = gridmetrics.Metrics1D(X)
        print('Grid metrics and Jacobian evaluation: Completed.')

    print('Grid transformation to curvilinear coordinate system:\
 Completed.')

    return X, Y

#**************************************************************************
def CalcTimeStep(CFL, diff, conv, dX, dY, Dimension, Model):
    '''Returns the time step size in the solution.
    '''
    #*************** DIFFUSION EQN. ******************
    if Model.upper() == 'DIFFUSION':
        dX2 = dX*dX
        if Dimension.upper() == '1D':
            TimeStep = CFL*dX2/diff
        elif Dimension.upper() == '2D':
            dY2 = dY*dY
            TimeStep = CFL*(1.0/((1/dX2) + (1/dY2)))/diff
    #*************** FIRST-ORDER WAVE EQN. *****************
    elif Model.upper() == 'FO_WAVE':
        if Dimension.upper() == '1D':
            TimeStep = CFL*dX/conv
    #*************** BURGERS EQN. *****************
    elif Model.upper() == 'BURGERS':
        if Dimension.upper() == '1D':
            TimeStep = CFL*dX
    print('Calculating time step size for the simulation: Completed.')

    return TimeStep

#**************************************************************************
def CalcMaxSteps(State, nMax, dT, simTime):
    '''Returns the max iteration/time steps for the solution.
    '''
    if State.upper() == 'TRANSIENT':
        if not simTime > 0.0: # simulation time can't be negative
            error_msg = 'ERROR: INCORRECT DATA FORMAT.\nCheck section\
 [STOP] --> key: [SIM_TIME] in\nfile:'
            raise Exception(f'{error_msg} {InFileName}.')
        try:
            MaxSteps = int(simTime/dT)
        except dT:
            raise Exception('No time step provided.')
    elif State.upper() == 'STEADY':
        MaxSteps = nMax
    print('Calculating maximum iterations/steps for the simulation:\
 Completed.')

    return MaxSteps

