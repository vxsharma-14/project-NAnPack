# coding: utf-8
import numpy as np

#*********************************************************************************
def BC(U, iMax, jMax, BCType):
    '''Assigns boundary condition at the walls, inlet and outlet'''

    if iMax % 2 == 0:
        iMid = int(iMax/2)
        jMid = int(jMax/2)
    else:
        iMid = int((iMax - 1)/2)
        jMid = int((jMax - 1)/2)

    i1 = iMid
    i2 = iMid+1

    j1 = jMid+5
    j2 = jMid+6
    
    U_0 = 0.0
    U_100 = 100.0

    if BCType == 'Dirichlet': # Dirichlet boundary conditions
        # along j = 1
        U[0:-1, 0] = U_100

        # along j = jMax
        U[0:-1, -1] = U_0

        # along i = 1
        U[0, 0:-1] = U_0

        # along i = iMax
        U[-1, 0:j1] = U_0

    elif BCType == 'Neumann': # Neumann boundary conditions

        # along j = 1
        U[0:i1, 0] = U_0
        U[i2-1:-1, 0] = U_100

        # along j = jMax
        U[0:-1, -1] = U_0

        # along i = 1
        U[0, 0:-1] = U_0

        # along i = iMax
        U[-1, 0:j1] = U_100
        U[-1, j2-1:-1] = U_0

    elif BCType == 'Mixed': # Mixed type boundary conditions

        # along j = 1
        U[0:i1, 0] = U_0
        U[i2-1:-1, 0] = U_100

        # along j = jMax
        U[0:-1, -1] = U_0

        # along i = 1
        U[0, 0:-1] = U_0

        # along i = iMax
        U[-1, 0:j1] = U_100
        U[-1, j2-1:-1] = U_0

    else:
        print('Incorrect boundary condition type specified')
        print('Try "Dirichlet", "Neumann" or "Mixed" instead')
        print('Program will quit now')
        sys.exit()
    
    return U

#*********************************************************************************
def Initial(iMax, jMax):
    '''Assigns initial condition within the domain'''

    U = np.zeros((iMax,jMax), dtype='float64')

    return U

#*********************************************************************************
def InputParam(InFileName='./Input/InputParameters.dat'):
    '''Read input parameters for the simulation set-up from a saved file

    Call signature;

        InputParam(InFileName)

    Parameters
    __________

    InFileName: str, default: 'InputParameters.dat'
                The string value representing the file to be read for
                simulation inputs. Use ".dat" or ".txt" file extensions
    '''

    import sys
    
    print(f'Reading Input Parameters from file {InFileName}.')
    InFile = open(InFileName, 'r')
    #InFile.seek(0)
    
    Line1 = (InFile.readline()).split()
    if Line1[0] != "EXPERIME":
        print(f'{Line1[0]}: Record Missing or Misplaced in Line 1')
        print(f'Check File {InFileName}')
        sys.exit(0)
        
    Line2 = (InFile.readline()).split()
    if Line2[0] != "GRIDPOINT":
        print(f'{Line2[0]}: Record Missing or Misplaced in Line 2')
        print(f'Check File {InFileName}')
        sys.exit(0)
    
    Line3 = (InFile.readline()).split()
    if Line3[0] != "DIMENSION":
        print(f'{Line3[0]}: Record Missing or Misplaced in Line 3')
        print(f'Check File {InFileName}')
        sys.exit(0)
        
    Line4 = (InFile.readline()).split()
    if Line4[0] != "COURANT":
        print(f'{Line4[0]}: Record Missing or Misplaced in Line 4')
        print(f'Check File {InFileName}')
        sys.exit(0)
        
    Line5 = (InFile.readline()).split()
    if Line5[0] != "CONVERG":
        print(f'{Line5[0]}: Record Missing or Misplaced in Line 5')
        print(f'Check File {InFileName}')
        sys.exit(0)
    
    Line6 = (InFile.readline()).split()
    if Line6[0] != "FILENAM":
        print(f'{Line6[0]}: Record Missing or Misplaced in Line 6')
        print(f'Check File {InFileName}')
        sys.exit(0)
    
    Line7 = (InFile.readline()).split()
    if Line7[0] != "FRAMESCAP":
        print(f'{Line7[0]}: Record Missing or Misplaced in Line 7')
        print(f'Check File {InFileName}')
        sys.exit(0)
    
    ExpNumber = Line1[1]
    iMax, jMax = [int(Line2[1]), int(Line2[2])]
    L, H = [float(Line3[1]), float(Line3[2])]
    CFL = float(Line4[1])
    ConvCriteria, nMax = [float(Line5[1]), int(Line5[2])]
    nWrite, OutFileName, nDisplay, HistFileName = [int(Line6[1]), Line6[2]\
                                                  ,int(Line6[3]), Line6[4]]
    FrameOpt, FrameWrite = [int(Line7[1]), (Line7[2])]
    print('Simulation set up done')

    OutFileName = './Output/' + OutFileName
    HistFileName = './Output/' + HistFileName
    
    # Return the data using dataclass method
    return (ExpNumber, iMax, jMax, L, H, ConvCriteria, nMax,\
            nWrite, OutFileName, nDisplay, HistFileName, FrameOpt,\
            FrameWrite)

#*********************************************************************************
if __name__ == '__main__':
    InputParam(InFileName)
