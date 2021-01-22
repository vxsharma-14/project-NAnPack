'''
+**************************************************************************
+**************************************************************************
+
+   FILE         writetofiles.py
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
def WriteSolutionToFile(init, n, U, OutFileName=None):
    '''Write simulation results of the entire domain to output file.

    Call Signature:

        WriteSolutionToFile(init, n, U)

    Parameters
    ----------

    init :

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
        OutFileName = init.OutFileName
    if n % init.nWrite == 0 or n == init.nMax:
        OutFile = open(OutFileName, "w")
        shapeU = U.shape # Obtain shape for Dimension
        if len(shapeU) == 1: # Dimension = 1D
            print(f'IMAX= {init.iMax}', file=OutFile)
            print(f'Grid Step= {init.dX:>6.4f}', file=OutFile)
            print(f'{"X":^14} {"U":^14}', file=OutFile)
            for i in range (0, init.iMax):
                print(f'{init.dX*i:>12.8e} {U[i]:>12.8e}', file=OutFile)
        elif len(shapeU) == 2: # Dimension = 2D
            print(f'IMAX= {init.iMax}, JMAX= {init.jMax}', file=OutFile)
            print(f'Grid Steps= {init.dX:>6.4f} {init.dY:>6.4f}',\
                  file=OutFile)
            print(f'{"X":^14} {"Y":^14} {"U":^14}', file=OutFile)
            for i in range (0, init.iMax):
                for j in range (0, init.jMax):
                    print(f'{init.dX*i:>12.8e} {init.dY*j:>12.8e}\
 {U[i][j]:>12.8e}', file=OutFile)

        #print('Writing 2D output: Completed')
        #print('Files saved:')
        #print(f'"{OutFileName}".')

        OutFile.close()

#**************************************************************************
def WriteConvHistToFile(init, n, Error):
    '''Write convergence history log

    Call Signature:

        WriteConvHistToFile(n, nMax, nDisplay, HistFName, Error, ConvCrit,\
                            State, Scheme)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    n : int

        Iteration level (time level).

    Error : float

            Error in the solution at each iteration level.
    '''
    # Write convergence data to file
    if n == 0 or n == 1:
        HistFile = open(init.HistFileName, "w")
        print(f'{"ITER":>7} {"ERROR":>15}', file=HistFile)
        print(f'{"----":>7} {"-----":>15}', file=HistFile)

    if n % init.nDisplay == 0 or n == init.nMax:
        HistFile = open(init.HistFileName, "a")
        print(f'{n:>7} {Error:>15.8f}', file=HistFile)

        HistFile.close()
        
    # Write convergence status to output file at the end if STEADY STATE
    if init.State.upper() == 'STEADY':
        if n == init.nMax and Error > init.ConvCrit:
            PrintMsg = f"\nSTATUS: NOT CONVERGED\nMAX. ITERATIONS= {init.nMax}\
\nNUMERICAL METHOD= {init.Scheme}"
            print(PrintMsg)
            print()
            HistFile = open(init.HistFileName, "a")
            print(PrintMsg, file=HistFile)
            print('Writing convergence log file: Completed')
            print('Files saved:')
            print(f'"{init.HistFileName}".')

            HistFile.close()

        elif Error < init.ConvCrit:
            PrintMsg = f'\nSTATUS: CONVERGED SOLUTION OBTAINED AT\nITERATIONS=\
 {n},\nMAX. ERROR= {init.ConvCrit}\nNUMERICAL METHOD= {init.Scheme}'
            print(PrintMsg)
            print()

            HistFile = open(init.HistFileName, "a")
            print(PrintMsg, file=HistFile)
            print('Writing convergence log file: Completed')
            print('Files saved:')
            print(f'"{init.HistFileName}".')

            HistFile.close()

    elif init.State.upper() == 'TRANSIENT' and n == init.nMax:
        PrintMsg = f'\nSTATUS: SOLUTION OBTAINED AT\nTIME LEVEL= {init.totTime}\
 s.\nTIME STEPS= {init.nMax}\nNUMERICAL METHOD= {init.Scheme}'
        print(PrintMsg)
        print()

        HistFile = open(init.HistFileName, "a")
        print(PrintMsg, file=HistFile)
        print('Writing convergence log file: Completed')
        print('Files saved:')
        print(f'"{init.HistFileName}".')

        HistFile.close()

#**************************************************************************
def WriteSolutionIn1DFormat(init, U, Out1DFName=None):
    '''Write 2D output data in 1D format at locations of X or Y.
    
    Call Signature:

        WriteSolutionIn1DFormat(init, U)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    U : numpy 2D array, float

       Solution of the dependent variable to be stored.

    Out1DFName : str

                 File name to save output. Default=None
    '''
    if not Out1DFName:
        Out1DFName = init.Out1DFName
    OutFile = open(Out1DFName, "w") # Open file to write
    line = [] # initialize a var to store linewise data
    # Format line # 1 in output file which print location points
    locations = ["{0:^12.2f}".format(l) for l in init.nodes]
    # Format line # 2 in output file which prints dashes
    dashes = ["{0:^12}".format('------') for l in init.nodes]

    if init.PrintNodesDir.upper() == 'X':
        iX = [int(i/init.dX) for i in init.nodes]
        print('Writing data in 1D format at X locations parallel to Y-axis.')
        print(f'{"Y":<5} {"X=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}',*dashes, sep=' ', file=OutFile)
        for j in range(0,init.jMax):
            del line[:] # Delete elements in the list for next Y location
            for i in iX:
                # Read output at different X locations at a constant Y and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float and format the output
                lines = ["{0:12.8f}".format(float(l)) for l in (line)]
            print(f'{init.dY*j:<8.4f}', *lines, sep=' ', file=OutFile)

    elif init.Direction.upper() == 'Y':
        jY = [int(j/init.dY) for j in init.nodes]
        print('Writing data in 1D format at Y locations parallel to X-axis.')
        print(f'{"X":<5} {"Y=":^2}', *locations, sep=' ', file=OutFile)
        print(f'{"---":<8}',*dashes, sep=' ', file=OutFile)
        for i in range(0,init.iMax):
            del line[:] # Delete elements in the list for next X location
            for j in jY:
                # Read output at different Y locations at a constant X and
                # store in one line
                line.append(U[i][j])
                # Convert "line" list elements to float and format the output
                lines = ["{0:12.8f}".format(float(l)) for l in (line)]
            print(f'{init.dX*i:<8.4f}', *lines, sep=' ', file=OutFile)

    print('Writing output in 1D format: Completed')
    print('Files saved:')
    print(f'"{Out1DFName}".')
    
    OutFile.close()
#**************************************************************************
