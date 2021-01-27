#**************************************************************************
def Metrics2D(X, Y, dXi=1.0, dEta=1.0):
    '''Returns the metrics of the transformation, XiX, XiY, EtaX and EtaY
    and the Jacobian of the transformation, J.
    '''
    shapeX = X.shape
    iMax, jMax = shapeX

    # Initialize
    XXi = X.copy()
    YXi = X.copy()
    XEta = X.copy()
    YEta = X.copy()
    XiX = X.copy()
    XiY = X.copy()
    EtaX = X.copy()
    EtaY = X.copy()
    JJ = X.copy()

    # At interior grid points
    XXi[1:-1,1:-1] = (X[2:,1:-1]-X[0:-2,1:-1]) / 2.0 / dXi
    YXi[1:-1,1:-1] = (Y[2:,1:-1]-Y[0:-2,1:-1]) / 2.0 / dXi

    XEta[1:-1,1:-1] = (X[1:-1,2:]-X[1:-1,0:-2]) / 2.0 / dEta
    YEta[1:-1,1:-1] = (Y[1:-1,2:]-Y[1:-1,0:-2]) / 2.0 / dEta

    # At boundary X = 0.0
    XXi[0,1:-1] = (-3.0*X[0,1:-1] + 4.0*X[1,1:-1] - X[2,1:-1]) / 2.0 / dXi
    YXi[0,1:-1] = (-3.0*Y[0,1:-1] + 4.0*Y[1,1:-1] - Y[2,1:-1]) / 2.0 / dXi

    XEta[0,1:-1] = (X[0,2:]-X[0,0:-2]) / 2.0 / dEta
    YEta[0,1:-1] = (Y[0,2:]-Y[0,0:-2]) / 2.0 / dEta

    # At boundary X = L
    XXi[-1,1:-1] = (3.0*X[-1,1:-1] - 4.0*X[-2,1:-1] + X[-3,1:-1]) / 2.0\
                   / dXi
    YXi[-1,1:-1] = (3.0*Y[-1,1:-1] - 4.0*Y[-2,1:-1] + Y[-3,1:-1]) / 2.0\
                   / dXi

    XEta[-1,1:-1] = (X[-1,2:]-X[-1,0:-2]) / 2.0 / dEta
    YEta[-1,1:-1] = (Y[-1,2:]-Y[-1,0:-2]) / 2.0 / dEta

    # At boundary Y = 0.0
    XXi[1:-1,0] = (X[2:,0]-X[0:-2,0]) / 2.0 / dXi
    YXi[1:-1,0] = (Y[2:,0]-Y[0:-2,0]) / 2.0 / dXi

    XEta[1:-1,0] = (-3.0*X[1:-1,0] + 4.0*X[1:-1,1] - X[1:-1,2]) / 2.0 / dEta
    YEta[1:-1,0] = (-3.0*Y[1:-1,0] + 4.0*Y[1:-1,1] - Y[1:-1,2]) / 2.0 / dEta

    # At boundary Y = H 
    XXi[1:-1,-1] = (X[2:,-1]-X[0:-2,-1]) / 2.0 / dXi
    YXi[1:-1,-1] = (Y[2:,-1]-Y[0:-2,-1]) / 2.0 / dXi

    XEta[1:-1,-1] = (3.0*X[1:-1,-1] - 4.0*X[1:-1,-2] + X[1:-1,-3]) / 2.0\
                    / dEta
    YEta[1:-1,-1] = (3.0*Y[1:-1,-1] - 4.0*Y[1:-1,-2] + Y[1:-1,-3]) / 2.0\
                    / dEta

    # At vertices
    # X=0.0, Y=0.0
    XXi[0,0] = (-3.0*X[0,0] + 4.0*X[1,0] - X[2,0]) / 2.0 / dXi
    YXi[0,0] = (-3.0*Y[0,0] + 4.0*Y[1,0] - Y[2,0]) / 2.0 / dXi

    XEta[0,0] = (-3.0*X[0,0] + 4.0*X[0,1] - X[0,2]) / 2.0 / dEta
    YEta[0,0] = (-3.0*Y[0,0] + 4.0*Y[0,1] - Y[0,2]) / 2.0 / dEta

    # X=L, Y=0.0
    XXi[-1,0] = (3.0*X[-1,0] - 4.0*X[-2,0] + X[-3,0]) / 2.0 / dXi
    YXi[-1,0] = (3.0*Y[-1,0] - 4.0*Y[-2,0] + Y[-3,0]) / 2.0 / dXi

    XEta[-1,0] = (-3.0*X[-1,0] + 4.0*X[-1,1] - X[-1,2]) / 2.0 / dEta
    YEta[-1,0] = (-3.0*Y[-1,0] + 4.0*Y[-1,1] - Y[-1,2]) / 2.0 / dEta

    # X=0.0, Y=H
    XXi[0,-1] = (-3.0*X[0,-1] + 4.0*X[1,-1] - X[2,-1]) / 2.0 / dXi
    YXi[0,-1] = (-3.0*Y[0,-1] + 4.0*Y[1,-1] - Y[2,-1]) / 2.0 / dXi

    XEta[0,-1] = (3.0*X[0,-1] - 4.0*X[0,-2] + X[0,-3]) / 2.0 / dEta
    YEta[0,-1] = (3.0*Y[0,-1] - 4.0*Y[0,-2] + Y[0,-3]) / 2.0 / dEta

    # X=L, Y=H
    XXi[-1,-1] = (3.0*X[-1,-1] - 4.0*X[-2,-1] + X[-3,-1]) / 2.0 / dXi
    YXi[-1,-1] = (3.0*Y[-1,-1] - 4.0*Y[-2,-1] + Y[-3,-1]) / 2.0 / dXi

    XEta[-1,-1] = (3.0*X[-1,-1] - 4.0*X[-1,-2] + X[-1,-3]) / 2.0 / dEta
    YEta[-1,-1] = (3.0*Y[-1,-1] - 4.0*Y[-1,-2] + Y[-1,-3]) / 2.0 / dEta

    # Evaluate metrics and Jacobian
    '''
    JJ[:,:] = 1.0 / (XXi[:,:]*YEta[:,:] - YXi[:,:]*XEta[:,:])
    XiX[:,:] = JJ[:,:] * YEta[:,:]
    XiY[:,:] = -JJ[:,:] * XEta[:,:]
    EtaX[:,:] = -JJ[:,:] * YXi[:,:]
    EtaY[:,:] = JJ[:,:] * XXi[:,:]
    '''
    for i in range(0,iMax):
        for j in range(0,jMax):
            JJ[i][j] = 1.0 / (XXi[i][j]*YEta[i][j] - YXi[i][j]*XEta[i][j])
            XiX[i][j] = JJ[i][j] * YEta[i][j]
            XiY[i][j] = -JJ[i][j] * XEta[i][j]
            EtaX[i][j] = -JJ[i][j] * YXi[i][j]
            EtaY[i][j] = JJ[i][j] * XXi[i][j]


    return XiX, XiY, EtaX, EtaY, JJ
