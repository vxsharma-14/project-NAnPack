"""A module for post-processing operations."""
#   ***********************************************************************
#
#   FILE         preprocess.py
#
#   AUTHOR       Dr. Vishal Sharma
#
#   VERSION      1.0.0-alpha4
#
#   WEBSITE      https://github.com/vxsharma-14/project-NAnPack
#
#   NAnPack Learner's Edition is distributed under the MIT License.
#
#   Copyright (c) 2020 Vishal Sharma
#
#   Permission is hereby granted, free of charge, to any person
#   obtaining a copy of this software and associated documentation
#   files (the "Software"), to deal in the Software without restriction,
#   including without limitation the rights to use, copy, modify, merge,
#   publish, distribute, sublicense, and/or sell copies of the Software,
#   and to permit persons to whom the Software is furnished to do so,
#   subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#
#   You should have received a copy of the MIT License along with
#   NAnPack Learner's Edition.
#
#   ***********************************************************************

import configparser
from .backend.checkconfig import CheckSections
from .backend.exceptions import InputFileError
# from .backend.exceptions import InvalidValueError


class RunConfig:
    """Class to read inputs from configuration file and set-up variables.

    This class is depenedent on the configparser package to read the
    inputs from the file and therefore, for compability,
    always use the most recent version of the config file provided
    with the package.

    Attributes
    ----------
    InFileName: str, Default= "./input/config.ini".

        The string value representing the file to be read for simulation
        inputs. Allowed file extensions: .ini
    """

    def __init__(self, InFileName="./input/config.ini"):
        """Class constructor for the RunConfig class.

        The contructor method calls the members functions to input
        parameters for the simulation set-up from
        the configuration file.

        The input properties for the simulation are stored in the class
        member variables that can be accessed using the
        RunConfig class constructor object.


        Parameters
        ----------
        InFileName: str

            Path to configuration file.
        """
        self.File = InFileName

        print("*******************************************************")
        print("*******************************************************")
        print("Starting configuration.")
        print()
        print("Searching for simulation configuration file in path:")
        print(f'"{InFileName}"')
        self.config = configparser.ConfigParser()
        dataset = self.config.read(InFileName)
        if dataset:
            print("SUCCESS: Configuration file parsing.")
        else:
            raise InputFileError("FileNotFound", InFileName)

        # Upon initialization -
        #       1. check - if all sections exist
        #       2. access numerical setup and and check inputs
        print("Checking whether all sections are included in config file.")
        CheckSections(self.config, self.File)
        # Access all other sections and set variables.
        self.ConfigSolverSetUp()
        self.ConfigGrid()
        self.ConfigInitial()
        self.ConfigBC()
        self.ConfigConstants()
        self.ConfigSimStop()
        self.ConfigOutput()
        self.DisplayConfig()

    def ConfigSolverSetUp(self):
        """Access numerical setup inputs from the SETUP section.

        Call signature :

            self.SolverSetUp()
        """
        from .backend.checkconfig import CheckSetupSection

        print("Checking numerical setup.")
        self.ExpId = self.config['SETUP']['EXPID']
        self.UnitSystem = self.config['SETUP']['UNITS_SYSTEM']
        self.Description = self.config['SETUP']['DESCRIPTION']
        self.State = self.config['SETUP']['STATE']
        self.Model = self.config['SETUP']['MODEL']
        self.Scheme = self.config['SETUP']['SCHEME']
        self.Dimension = self.config['SETUP']['DIMENSION']
        CheckSetupSection(self.config, self.State, self.Model,
                          self.Dimension, self.File)

    def ConfigGrid(self):
        """Access meshing inputs from the DOMAIN and MESH sections.

        Call signature :

            self.GridGen()
        """
        from . import grid
        # ***************** DOMAIN SPEC *****************
        self.Length = float(self.config['DOMAIN']['LENGTH'])
        self.Height = float(self.config['DOMAIN']['HEIGHT'])
        print("Accessing domain geometry configuration: Completed")
        # **************** MESH SPEC *****************
        GridfromFile = self.config['MESH']['GRID_FROM_FILE?']
        if GridfromFile.upper() == 'YES':
            self.GridFName = self.config['MESH']['GRID_FNAME']
            if self.GridFName.lower() == 'none':
                raise InputFileError("FileNotFound", self.GridFName)
            else:
                print("Functionality not available at this time")
                print("Proceeding using other inputs.")
        GridAutoCalc = self.config['MESH']['GRID_AUTO_CALC?']
        if GridAutoCalc.upper() == 'YES':
            self.dX = float(self.config['MESH']['dX'])
            self.dY = float(self.config['MESH']['dY'])
            print("Accessing meshing configuration: Completed.")
            self.iMax, self.jMax = grid.ComputeGridPoints(self.Dimension,
                                                          self.Length,
                                                          self.dX,
                                                          self.Height,
                                                          self.dY)
        elif GridAutoCalc.upper() == 'NO':
            self.iMax = int(self.config['MESH']['iMax'])
            self.jMax = int(self.config['MESH']['jMax'])
            print("Accessing meshing configuration: Completed.")
            self.dX, self.dY = grid.ComputeGridSteps(self.Dimension,
                                                     self.Length,
                                                     self.iMax,
                                                     self.Height,
                                                     self.jMax)

    def ConfigInitial(self):
        """Access intial condition inputs from the IC section.

        Returns the dependent variables with the initial values.

        Call signature :

            RunConfig.Initial()

        Returns
        -------
        U: 1D or 2D array

            Initial conditions for the dependent variable.
        """
        import nanpack.backend.initialize as init
        # *********** INITIAL CONDITIONS *************
        # Firstly, initialize the dependent variable U which is
        # required before assigning any type of intial conditions,
        # whether cold-start or restart.
        self.U = init.InitialCondition(self.Dimension, self.iMax,
                                       self.jMax)
        self.StartOpt = self.config['IC']['START_OPT']
        if self.StartOpt.upper() == 'RESTART':
            print("Functionality does not exists in this version.")
            pass
        if 1 == 0:
            self.RestartFile = self.config['IC']['RESTART_FILE']
            print("Accessing initial condition settings: Completed.")
            if self.RestartFile.lower() == "none":
                raise InputFileError("FileNotFound", self.RestartFile)
                ch = int(input("Do you want to proceed with:\
\n\t1. COLD START conditions, or\n\t2. Use a system default\
 RESTART filename?\nENTER 1 or 2.\n"))
                if ch == 1:
                    # Initialize with zero
                    self.StartOpt = "COLD-START"
                    self.U = init.InitialCondition(self.Dimension,
                                                   self.iMax, self.jMax)
                elif ch == 2:
                    # Use defualt file name
                    # InitFile = './output/restart.dat'
                    print("This functionality is not available at this\
 time.")
                    print("Proceeding to solve using cold-start\
 conditions.")

            print("Accessing initial condition settings: Completed.")

        return self.U

    def ConfigBC(self):
        """Access boundary condition inputs from the BC section.

        Update the boundary conditions of the dependent variable and
        return.

        Call signature :

            RunConfig.BC()

        Returns
        -------
        U: 2D array

            Updated boundary conditions for the dependent variable.
        """
        import nanpack.backend.boundary as bound
        # *********** BOUNDARY CONDITIONS *************
        self.BCfromFile = self.config['BC']['BC_FROM_FILE?']
        self.BCFileName = self.config['BC']['BC_FILE_NAME']
        print("Accessing boundary condition settings: Completed")
        if self.BCfromFile.upper() == 'YES':
            self.BC = bound.ReadBCfromFile(self.BCFileName)
            # Call function to assign 2D BC
            self.U = bound.BC2D(self.U, self.BC, self.dX, self.dY)
            print("Boundary conditions assignment: Completed.")
        elif self.BCfromFile.upper() == 'NO':
            self.U = self.U

        return self.U

    def ConfigConstants(self):
        """Access constant inputs from the CONST section.

        Call signature :

            RunConfig.Constants()
        """
        # *********** CONSTANT COEFFICIENTS ***********
        self.CFL = float(self.config['CONST']['CFL'])
        self.conv = float(self.config['CONST']['CONV'])
        self.diff = float(self.config['CONST']['DIFF'])
        print("Accessing constant data: Completed.")

    def ConfigSimStop(self):
        """Access simulation stop setting inputs from the STOP section.

        Call signature :

            RunConfig.ConfigSimStop()
        """
        import nanpack.grid as grid
        # *********** SIM STOP SETTINGS ***********
        self.totTime = float(self.config['STOP']['SIM_TIME'])
        if self.State.upper() == 'STEADY':
            self.ConvCrit = float(self.config['STOP']['CONV_CRIT'])
        elif self.State.upper() == 'TRANSIENT':
            self.ConvCrit = -0.01
        nMax = int(self.config['STOP']['nMAX'])
        # Execute this block for:
        # diffusion eq., first-order wave eq. and Burgers eq.
        if not self.Model.upper() == 'POISSONS':
            self.dT = grid.CalcTimeStep(self.CFL, self.diff, self.conv,
                                        self.dX, self.dY,
                                        self.Dimension, self.Model)
            self.nMax = grid.CalcMaxSteps(self.State, nMax, self.dT,
                                          self.totTime)
        # Execute this block for Poissons eq.
        elif self.Model.upper() == 'POISSONS':
            self.nMax = nMax

        print("Accessing simulation stop settings: Completed.")

    def ConfigOutput(self):
        """Access output configurations from the OUTPUT section.

        Call signature :

            RunConfig.ConfigOutput()
        """
        # ************* OUTPUT INFORMATION *************
        self.HistFileName = self.config['OUTPUT']['HIST_FILE_NAME']
        self.RestartFile = self.config['OUTPUT']['RESTART_FNAME']
        self.OutFileName = self.config['OUTPUT']['RESULT_FNAME']
        self.nWrite = int(self.config['OUTPUT']['WRITE_EVERY'])
        self.nDisplay = int(self.config['OUTPUT']['DISPLAY_EVERY'])
        self.SaveforAnim = self.config['OUTPUT']['SAVE_FOR_ANIM?']
        if self.SaveforAnim.upper() == 'YES':
            self.nAnime = int(self.config['OUTPUT']['SAVE_EVERY'])
        self.Save1DOut = self.config['OUTPUT']['SAVE_1D_OUTPUT?']
        if self.Save1DOut.upper() == 'YES':
            nodeX = None
            nodeY = None
            try:
                nodeX = self.config['OUTPUT']['X']
                self.nodes = [float(node) for node in nodeX.split(',')]
                self.PrintNodesDir = 'X'
            except:
                nodeY = self.config['OUTPUT']['Y']
                self.nodes = [float(node) for node in nodeY.split(',')]
                self.PrintNodesDir = 'Y'
            self.Out1DFName = self.config['OUTPUT']['SAVE1D_FILENAME']
        print("Accessing settings for storing outputs: Completed.")

    def DisplayConfig(self):
        """Display the configuration to the user for verification.

        Call signature :

            RunConfig.DisplayConfig()
        """
        # ************* PRINT CONFIGURATIONS *************
        print()
        print('**********************************************************')
        print(f'CASE DESCRIPTION                {self.Description}')
        print(f'SOLVER STATE                    {self.State}')
        print(f'MODEL EQUATION                  {self.Model}')
        print(f'DOMAIN DIMENSION                {self.Dimension}')
        print(f'    LENGTH                      {self.Length}')
        if self.Dimension.upper() == '2D':
            print(f'    HEIGHT                      {self.Height}')
        print('GRID STEP SIZE')
        print(f'    dX                          {self.dX:5.3f}')
        if self.Dimension.upper() == '2D':
            print(f'    dY                          {self.dY:5.3f}')
        if not self.Model == 'POISSONS':
            print(f'TIME STEP                       {self.dT:5.3f}')
        print('GRID POINTS')
        print(f'    along X                     {self.iMax}')
        if self.Dimension.upper() == '2D':
            print(f'    along Y                     {self.jMax}')
        if self.Model.upper() == 'DIFFUSION':
            print(f'DIFFUSION CONST.                {self.diff:6.4e}')
            print(f'DIFFUSION NUMBER                {self.CFL}')
        elif self.Model.upper() == 'FO_WAVE':
            print(f'CONVECTION CONST.               {self.conv}')
            print(f'COURANT NUMBER                  {self.CFL}')
        elif self.Model.upper() == 'BURGERS':
            print(f'COURANT NUMBER                  {self.CFL}')
        if self.State.upper() == 'STEADY':
            print(f'CONVERGENCE CRIT                {self.ConvCrit}')
            print('MAXIMUM ITERATIONS')
            print(f'IF CONVERGENCE NOT OBSERVED     {self.nMax}')
        elif self.State.upper() == 'TRANSIENT':
            print(f'TOTAL SIMULATION TIME           {self.totTime}')
            print(f'NUMBER OF TIME STEPS            {self.nMax}')
        if self.BCfromFile.upper() == 'YES':
            print(f'BC COFIGURATION FILE            "{self.BCFileName}"')
        print(f'START CONDITION                 {self.StartOpt}')
        print('**********************************************************')
        print('SUCEESS: Configuration completed.')
        print()


def BC2D(U, BC, dX, dY):
    """Assign boundary conditions. Add more info."""
    from .backend.boundary import BC2D
    u = BC2D(U, BC, dX, dY)

    return u


def CourantNumber(CFL, Dimension):
    """Return the Courant Number.(this function needs modification)."""
    if Dimension.upper() == '1D':
        convX = CFL  # convection coefficient for the x-term.

    elif Dimension.upper() == '2D':
        convX = CFL  # convection coefficient for the x-term.
        convY = CFL  # convection coefficient for the y-term.
    print('Calculating coefficient in the convective equation: Completed.')

    return convX, convY


def DiffusionNumbers(Dimension, diff, dT, dX, dY=None):
    """Return the diffusion numbers along X and Y directions.

    Call signature:

        DiffusionNumbers(Dimension, diff, dT, dX, dY)

    Parameters
    ----------
    Dimension : string

        Dimension of the simulation domain.
        Available Options are "1D" or "2D"

    diff: float

        Constant coefficient in the diffusion equation, such as kinematic
        viscosity.

    dT: float

        Time step size.

    dX: float

        Grid step size in x-direction.

    dY: float

        Grid step size in y-direction. Default=None for 1D applications.

    Returns
    -------
    diffX, diffY : float values

        Diffusion numbers along X and Y axis, respectively.
    """
    dX2 = dX*dX
    # diffusion coefficient for the x term.
    diffX = diff*dT/dX2
    diffY = None
    if Dimension.upper() == "2D":
        dY2 = dY*dY
        diffY = diff*dT/dY2  # diffusion coefficient for y term.
    print("Calculating diffusion numbers: Completed.")

    return diffX, diffY
