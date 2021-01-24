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
