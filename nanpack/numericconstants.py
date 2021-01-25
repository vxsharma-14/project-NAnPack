#**************************************************************************
def CourantNumber(CFL):
    ''' Returns the Courant Number.
    '''
    if Dimension.upper() == '1D':
        convX = CFL # convection coefficient for the x-term.

    elif Dimension.upper() == '2D':
        convX = CFL # convection coefficient for the x-term.
        convY = CFL # convection coefficient for the y-term.
    print('Calculating coefficient in the convective equation: Completed.')

    return convX, convY

#**************************************************************************
def DiffusionNumbers(Dimension, diff, dT, dX, dY=None):
    ''' Returns the diffusion numbers along X and Y directions.

    Call signature:

        DiffusionNumbers(Dimension, diff, dT, dX, dY)

    Parameters
    ----------

    Dimension : string

                Dimension of the simulation domain.
                Available Options = "1D" or "2D"

    Return
    ------

    diffX, diffY : float values

                   Diffusion numbers along X and Y axis, respectively.
    '''
    dX2 = dX*dX
    # diffusion coefficient for the x term.
    diffX = diff*dT/dX2
    diffY = None
    if Dimension.upper() == '2D':
        dY2 = dY*dY
        diffY = diff*dT/dY2 # diffusion coefficient for y term.
    print('Calculating diffusion numbers: Completed.')
        
    return diffX, diffY

#**************************************************************************
