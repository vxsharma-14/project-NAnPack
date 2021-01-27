class nondimensionalize:
    '''A class to non-dimensionalize various numerical parameters.'''
#**************************************************************************
    def __init__(self):
        '''Initialize class.'''

    def ndXgrid(self, X, L):
        '''Non-dimensionalize x point locations along X axis.'''
        self.xstar[:] = X[:]/L
        return self.xstar

    def ndYgrid(self, Y, L):
        '''Non-dimensionalize y point locations along Y axis.'''
        self.ystar[:] = Y[:]/L
        return self.ystar

    def ndvelU(self, U, L, nu):
        '''Non-dimensionalize velocity U.'''
        shapeU = U.shape
        if len(shapeU) == 1:
            im = shapeU
            self.ustar = [U[i]*L/nu for i in range(im)]
        elif len(shapeU) == 2:
            im, jm = shapeU
            self.ustar = [[U[i][j]*L/nu for j in range(jm)]\
                          for i in range(im)]
        return self.ustar

    def ndtime(self, t, L, nu):
        '''Non-dimensionalize time t.'''
        L2 = L*L
        self.tstar = nu*t/L2
        return self.tstar

#**************************************************************************
class dimensionalize:
    '''A class to dimensionalize various numerical parameters.'''

    def __init__(self):
        '''Initialize class.'''

    def dimXgrid(self, xstar, L):
        '''Dimensionalize x point locations along X axis.'''
        self.x[:] = xstar[:]*L
        return self.x

    def dimYgrid(self, ystar, L):
        '''Dimensionalize y point locations along Y axis.'''
        self.y[:] = ystar[:]*L
        return self.ystar

    def dimvelU(self, ustar, L, nu):
        '''Dimensionalize velocity U.'''
        shapeU = U.shape
        if len(shapeU) == 1:
            im = shapeU
            self.u = [ustar[i]*nu/L for i in range(im)]
        elif len(shapeU) == 2:
            im, jm = shapeU
            self.u = [[ustar[i][j]*nu/L for j in range(jm)]\
                      for i in range(im)]
        return self.u

    def dimtime(self, tstar, L, nu):
        '''Non-dimensionalize time t.'''
        L2 = L*L
        self.t = tstar*L2/nu
        return self.t
