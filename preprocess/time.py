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
