#**************************************************************************
def RectangularGrid(dX, iMax, dY=None, jMax=None):
    '''Compute X and Y grid point locations in a rectangular, uniform
    mesh in a cartesian coordinate system.
    '''
    print('Calculating X and Y locations of all grid points within\
 the mesh.')
    
    if isinstance(dY, float) and isinstance(jMax, int):
        X = [[i*dX for j in range (0,jMax)] for i in range (0,iMax)] 
        Y = [[j*dY for j in range (0,jMax)] for i in range (0,iMax)]
    else:
        X = [i*dX for i in range (0,iMax)]
        Y = 0.0
    print(f'Uniform rectangular grid generation in a cartesian\
 coordinate system: Completed.')

    return X, Y

#**************************************************************************
def ComputeGridPoints(Dimension, Length, delX, Height=None, delY=None):
    '''Compute grid points along X and Y direction in the mesh.
    '''
    iMax = int(Length/delX) + 1
    if Dimension.upper() == '2D':
        jMax = int(Height/delY) + 1
    else:
        jMax = 0
    print('Calculating grid size: Completed.')
    

    return iMax, jMax

#**************************************************************************
def ComputeGridSteps(Dimension, Length, iMax, Height=None, jMax=None):
    '''Compute uniform grid steps size along X and Y direction in the mesh.
    '''
    delX = Length/(iMax - 1)
    if Dimension.upper() == '2D':
        delY = Height/(jMax - 1)
    else:
        delY = 0.0
    print('Calculating grid step size: Completed.')

    return delX, delY

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

