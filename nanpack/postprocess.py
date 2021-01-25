'''
+**************************************************************************
+**************************************************************************
+
+   FILE         postprocess.py
+
+   AUTHOR       Vishal Sharma
+
+   VERSION      1.0.0.dev1
+
+   WEBSITE      https://vxsharma-14.github.io/NAnPack/
+
+   NAnPack Learner's Edition is distributed under the MIT License.
+
+   Copyright (c) 2020 Vishal Sharma
+
+   Permission is hereby granted, free of charge, to any person
+   obtaining a copy of this software and associated documentation
+   files (the "Software"), to deal in the Software without restriction,
+   including without limitation the rights to use, copy, modify, merge,
+   publish, distribute, sublicense, and/or sell copies of the Software,
+   and to permit persons to whom the Software is furnished to do so,
+   subject to the following conditions:
+
+   The above copyright notice and this permission notice shall be
+   included in all copies or substantial portions of the Software.
+
+   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
+   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
+   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
+   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
+   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
+   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
+   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
+   SOFTWARE.
+
+   You should have received a copy of the MIT License along with
+   NAnPack Learner's Edition.
+
+**************************************************************************
+**************************************************************************
'''
#**************************************************************************
def FourthOrderDamping(U, DampCoeff):
    '''Returns the fourth-order damping term.

    Calculated using Equation 6-77 in CFD Vol. 1 by Hoffmann.

    Call signature:

        FourthOrderDamping(U, DampCoeff)

    Parameters
    ----------

    U : 1D or 2D array

        The dependent variable from time level (n) within the domain.

    DampCoeff : float

                Damping coefficient. The value must be selected in the
                range of 0 to 0.125 for a stable solution.

    Returns
    -------

    D : 1D or 2D array

        The fourth-order damping term within the entire domain.
    '''
    D = U.copy() # Initialize D
    D[2:-2] = -DampCoeff*(U[0:-4] - 4.0*U[1:-3] + 6.0*U[2:-2]\
                          - 4.0*U[3:-1] + U[4:])

    return D

#**************************************************************************
def SecondOrderDamping(U, DampCoeff):
    '''Returns the second-order damping term.

    Calculated using Equation 6-80 in CFD Vol. 1 by Hoffmann.

    Call signature:

        SecondOrderDamping(U, DampCoeff)

    Parameters
    ----------

    U : 1D or 2D array

        The dependent variable from time level (n) within the domain.

    DampCoeff: float

               Damping coefficient.

    Returns
    -------

    D : 1D or 2D array

        The second-order damping term within the entire domain.
    '''
    D = U.copy() # Initialize D
    D[1:-1] = -DampCoeff*(U[0:-2] - 2.0*U[1:-1] + U[2:])
    
    return D

#**************************************************************************
def LInfNormError(U, Uold):
    '''Calculate L-infinity Norm error.
    '''
    shapeU = U.shape # Obtain shape for Dimension
    if len(shapeU) == 2: # Dimension = 2D
        Abs_err = abs(Uold[1:,1:] - U[1:,1:]).sum(axis=1)
        Error = max(Abs_err)

    elif len(shapeU) == 1: # Dimension = 1D
        Abs_err = abs(Uold[1:] - U[1:]).sum(axis=1)
        Error = max(Abs_err)
    
    return Error

#**************************************************************************
def AbsoluteError(U, Uold):
    '''Calculate absolute error.
    '''
    shapeU = U.shape # Obtain shape for Dimension
    if len(shapeU) == 2: # Dimension = 2D
        Error = abs(Uold[1:-1,1:-1] - U[1:-1,1:-1]).sum()

    elif len(shapeU) == 1: # Dimension = 1D
        Error = abs(Uold[1:-1] - U[1:-1]).sum()

    return Error

#**************************************************************************
def MonitorConvergence(n, nDisplay, Error):
    '''Monitor convergence of data
    '''
    
    if n == 0 or n == 1:
        print(f'{"ITER":>7} {"ERROR":>15}')
        print(f'{"----":>7} {"-----":>15}')
    
    if n % nDisplay == 0:
        print(f'{n:>7} {Error:>15.8f}')

#**************************************************************************
def WriteSolutionToFile(cfg, n, U, OutFileName=None):
    '''Write simulation results of the entire domain to output file.

    Call Signature:

        WriteSolutionToFile(cfg, n, U)

    Parameters
    ----------

    cfg :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    n : int

        Iteration level/time step level.

    U : 1D or 2D array, float

       Solution of the dependent variable to be stored.

    OutFileName : str

                  File name to save output. Default=None
    '''
    if not OutFileName:
        OutFileName = cfg.OutFileName
    if n % cfg.nWrite == 0 or n == cfg.nMax:
        OutFile = open(OutFileName, "w")
        shapeU = U.shape # Obtain shape for Dimension
        if len(shapeU) == 1: # Dimension = 1D
            print(f'IMAX= {cfg.iMax}', file=OutFile)
            print(f'Grid Step= {cfg.dX:>6.4f}', file=OutFile)
            print(f'{"X":^14} {"U":^14}', file=OutFile)
            for i in range (0, cfg.iMax):
                print(f'{cfg.dX*i:>12.8e} {U[i]:>12.8e}', file=OutFile)
        elif len(shapeU) == 2: # Dimension = 2D
            print(f'IMAX= {cfg.iMax}, JMAX= {cfg.jMax}', file=OutFile)
            print(f'Grid Steps= {cfg.dX:>6.4f} {cfg.dY:>6.4f}',\
                  file=OutFile)
            print(f'{"X":^14} {"Y":^14} {"U":^14}', file=OutFile)
            for i in range (0, cfg.iMax):
                for j in range (0, cfg.jMax):
                    print(f'{cfg.dX*i:>12.8e} {cfg.dY*j:>12.8e}\
 {U[i][j]:>12.8e}', file=OutFile)

        OutFile.close()

#**************************************************************************
def WriteConvHistToFile(cfg, n, Error):
    '''Write convergence history log

    Call Signature:

        WriteConvHistToFile(cfg, n, Error)

    Parameters
    ----------

    cfg :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    n : int

        Iteration level (time level).

    Error : float

            Error in the solution at each iteration level.
    '''
    # Write convergence data to file
    if n == 0 or n == 1:
        HistFile = open(cfg.HistFileName, "w")
        print(f'{"ITER":>7} {"ERROR":>15}', file=HistFile)
        print(f'{"----":>7} {"-----":>15}', file=HistFile)

    if n % cfg.nDisplay == 0 or n == cfg.nMax:
        HistFile = open(cfg.HistFileName, "a")
        print(f'{n:>7} {Error:>15.8f}', file=HistFile)

        HistFile.close()
        
    # Write convergence status to output file at the end if STEADY STATE
    if cfg.State.upper() == 'STEADY':
        if n == cfg.nMax and Error > cfg.ConvCrit:
            PrintMsg = f"\nSTATUS: NOT CONVERGED\nMAX. ITERATIONS= {cfg.nMax}\
\nNUMERICAL METHOD= {cfg.Scheme}"
            print(PrintMsg)
            print()
            HistFile = open(cfg.HistFileName, "a")
            print(PrintMsg, file=HistFile)
            print('Writing convergence log file: Completed')
            print('Files saved:')
            print(f'"{cfg.HistFileName}".')

            HistFile.close()

        elif Error < cfg.ConvCrit:
            PrintMsg = f'\nSTATUS: CONVERGED SOLUTION OBTAINED AT\nITERATIONS=\
 {n},\nMAX. ERROR= {cfg.ConvCrit}\nNUMERICAL METHOD= {cfg.Scheme}'
            print(PrintMsg)
            print()

            HistFile = open(cfg.HistFileName, "a")
            print(PrintMsg, file=HistFile)
            print('Writing convergence log file: Completed')
            print('Files saved:')
            print(f'"{cfg.HistFileName}".')

            HistFile.close()

    elif cfg.State.upper() == 'TRANSIENT' and n == cfg.nMax:
        PrintMsg = f'\nSTATUS: SOLUTION OBTAINED AT\nTIME LEVEL= {cfg.totTime}\
 s.\nTIME STEPS= {cfg.nMax}\nNUMERICAL METHOD= {cfg.Scheme}'
        print(PrintMsg)
        print()

        HistFile = open(cfg.HistFileName, "a")
        print(PrintMsg, file=HistFile)
        print('Writing convergence log file: Completed')
        print('Files saved:')
        print(f'"{cfg.HistFileName}".')

        HistFile.close()

#**************************************************************************
def WriteSolutionIn1DFormat(cfg, U, Out1DFName=None):
    '''Write 2D output data in 1D format at locations of X or Y.
    
    Call Signature:

        WriteSolutionIn1DFormat(cfg, U)

    Parameters
    ----------

    cfg :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    U : numpy 2D array, float

       Solution of the dependent variable to be stored.

    Out1DFName : str

                 File name to save output. Default=None
    '''
    if not Out1DFName:
        Out1DFName = cfg.Out1DFName
    OutFile = open(Out1DFName, "w") # Open file to write
    line = [] # initialize a var to store linewise data
    # Format line # 1 in output file which print location points
    locations = ["{0:^12.2f}".format(l) for l in cfg.nodes]
    # Format line # 2 in output file which prints dashes
    dashes = ["{0:^12}".format('------') for l in cfg.nodes]

    if cfg.PrintNodesDir.upper() == 'X':
        iX = [int(i/cfg.dX) for i in cfg.nodes]
        print('Writing data in 1D format at X locations parallel to Y-axis.')
        print(f'{"Y":<5} {"X=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}',*dashes, sep=' ', file=OutFile)
        for j in range(0,cfg.jMax):
            del line[:] # Delete elements in the list for next Y location
            for i in iX:
                # Read output at different X locations at a constant Y and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float and format the output
                lines = ["{0:12.8f}".format(float(l)) for l in (line)]
            print(f'{cfg.dY*j:<8.4f}', *lines, sep=' ', file=OutFile)

    elif cfg.Direction.upper() == 'Y':
        jY = [int(j/cfg.dY) for j in cfg.nodes]
        print('Writing data in 1D format at Y locations parallel to X-axis.')
        print(f'{"X":<5} {"Y=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}',*dashes, sep=' ', file=OutFile)
        for i in range(0,cfg.iMax):
            del line[:] # Delete elements in the list for next X location
            for j in jY:
                # Read output at different Y locations at a constant X and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float and format the output
                lines = ["{0:12.8f}".format(float(l)) for l in (line)]
            print(f'{cfg.dX*i:<8.4f}', *lines, sep=' ', file=OutFile)

    print('Writing output in 1D format: Completed')
    print('Files saved:')
    print(f'"{Out1DFName}".')
    
    OutFile.close()
#**************************************************************************
