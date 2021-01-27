#**************************************************************************
def ExternalFlowClustering(Xi,Eta):
    '''Returns the clustered grid near the X-lo or Y-lo wall boundaries.'''
    shapeX = X.shape()
    iMax, jMax = shapeX

    # If clustering is required near Y = 0.0 surface in positive Y
    # direction.
    Beta1 = (Beta+1.0) / (Beta-1.0)

    X[:,:] = Xi[:,:]
    Y[:,:] = Height * ((Beta+1.0) - (Beta-1.0)) * Beta1**(1.0-Eta[:,:])\
             / (Beta1**(1.0-Eta[:,:]) + 1.0)

    return XiX, XiY, EtaX, EtaY, JJ
