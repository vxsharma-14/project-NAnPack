class FetchOptions:
    '''A class to print/fetch available options for the user to specify in
    the config file or as arguments in the function call.
    '''
#**************************************************************************
    def __init__(self):
        '''Initialize class.'''

    def UnitOptions(self):
        '''Print a list of options avaible for the UNIT key in input
        file.
        '''
        self.options = ['SI', 'BRITISH']
        return self.options

    def StateOptions(self):
        '''Print a list of options avaible for the STATE key in input
        file.
        '''
        self.options = ['STEADY', 'TRANSIENT']
        return self.options

    def ModelOptions(self):
        '''Print a list of options avaible for the MODEL key in input
        file.
        '''
        self.options = ['DIFFUSION', 'POISSONS', 'FO_WAVE', 'BURGERS']
        return self.options

    def SolverOptions(self):
        '''Print a list of options avaible for the SOLVER key in input
        file.
        '''

    def DimensionOptions(self):
        '''Print a list of options avaible for the DIMENSION key in input
        file.
        '''
        self.options = ['1D', '2D']
        return self.options

    def TVDLimiterFunctionOptions(self):
        '''Print a list of options avaible for the LIMITERFUNC argument in
        the call to function SecondOrderTVD().
        '''
        self.options = ["Harten-Yee-Upwind", "Modified-Harten-Yee-Upwind",\
                        "Roe-Sweby-Upwind", "Davis-Yee-Symmetric"]
        return self.options
    
