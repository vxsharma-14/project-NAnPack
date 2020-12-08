# coding: utf-8
import sys

#*******************************************************************************
def InputParam(InFileName='./Input/SimulationSetUp.dat'):
    '''Read input parameters for the simulation set-up from a saved file

    Call signature;

        InputParam(InFileName)

    Parameters
    __________

    InFileName: str, default: './Input/SimulationSetUp.dat'
                The string value representing the file to be read for
                simulation inputs. Use ".dat" or ".txt" file extensions
    ''' 
    global iMax, jMax
    global Length, Height

    print(f'Reading Input Parameters from file {InFileName}.')
    InFile = open(InFileName, 'r')
    #InFile.seek(0)
    
    Line1 = (InFile.readline()).split()
    if Line1[0] != "EXPERIME":
        print(f'{Line1[0]}: Record Missing or Misplaced in Line 1')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
        
    Line2 = (InFile.readline()).split()
    if Line2[0] != "GRIDPOINT":
        print(f'{Line2[0]}: Record Missing or Misplaced in Line 2')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
    
    Line3 = (InFile.readline()).split()
    if Line3[0] != "DIMENSION":
        print(f'{Line3[0]}: Record Missing or Misplaced in Line 3')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
        
    Line4 = (InFile.readline()).split()
    if Line4[0] != "CONVERG":
        print(f'{Line4[0]}: Record Missing or Misplaced in Line 5')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
    
    Line5 = (InFile.readline()).split()
    if Line5[0] != "FILENAM":
        print(f'{Line5[0]}: Record Missing or Misplaced in Line 6')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
    
    Line6 = (InFile.readline()).split()
    if Line6[0] != "FRAMESCAP":
        print(f'{Line6[0]}: Record Missing or Misplaced in Line 7')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
    
    ExpNumber = Line1[1]
    iMax, jMax = [int(Line2[1]), int(Line2[2])]
    Length, Height = [float(Line3[1]), float(Line3[2])]
    ConvCriteria, nMax = [float(Line4[1]), int(Line4[2])]
    nWrite, OutFileName, nDisplay, HistFileName = [int(Line5[1]), Line5[2]\
                                                  ,int(Line5[3]), Line5[4]]
    FrameOpt, FrameWrite = [int(Line6[1]), int(Line6[2])]
    print('Simulation set up done.')

    OutFileName = './Output/' + OutFileName
    HistFileName = './Output/' + HistFileName
    
    # Return the data using dataclass method
    return (ExpNumber, ConvCriteria, nMax, nWrite, OutFileName,\
            nDisplay, HistFileName, FrameOpt, FrameWrite)

#*********************************************************************************
def Coefficients(InFileName='./Input/Coefficients.dat'):
    '''Read the coefficients - diffusivity that occur in the parabolic
       PDEs and the CFL number provided as user input
       
    Call signature;

        ReadCoefficients(InFileName)

    Parameters
    __________

    InFileName: str, default: './Input/Coefficients.dat'
                The string value representing the file to be read for
                simulation inputs. Use ".dat" or ".txt" file extensions
    '''
    
    print(f'Reading Coefficients from file {InFileName}.')
    InFile = open(InFileName, 'r')
    
    Line1 = InFile.readline()
        
    Line2 = (InFile.readline()).split()
    if Line2[0] != "COURANT":
        print(f'{Line2[0]}: Record Missing or Misplaced in Line 2')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
    
    Line3 = (InFile.readline()).split()
    if Line3[0] != "DIFFUSIV":
        print(f'{Line3[0]}: Record Missing or Misplaced in Line 3')
        print(f'Check File {InFileName}')
        sys.exit(STOP)
        
    print('Coefficient values assigned.')
    
    CFL = float(Line2[1])
    Diff = float(Line3[1])
    
    # Return the data using dataclass method
    return (CFL, Diff)

#*******************************************************************************
def Read1DNodesFile(InFileName='./Input/NodesFor1DOutput.dat'):
    '''Read the user inputs for storing 1D data to a file
       
    Call signature;

        Read1DNodesFile(InFileName)

    Parameters
    __________

    InFileName: str, default: './Input/NodesFor1DOutput.dat'
                The string value representing the file to be read.
                Use ".dat" or ".txt" file extensions
    '''

    print(f'Reading nodes data from file {InFileName}.')
    InFile = open(InFileName, 'r')
    
    Line1 = (InFile.readline()).split()
    if Line1[0] != "FILENAM1D":
        print(f'{Line1[0]}: Record Missing or Misplaced in Line 1')
        print(f'Check File {InFileName}')
        sys.exit('STOP')
        
    Line2 = (InFile.readline()).split()
    if Line2[0] != "NODES":
        print(f'{Line2[0]}: Record Missing or Misplaced in Line 2')
        print(f'Check File {InFileName}')
        sys.exit('STOP')

    print('Nodes data successfully read.')
    
    OutFile = Line1[1]
    X = [int(Line2[1]), int(Line2[2]), int(Line2[3]),\
         int(Line2[4]), int(Line2[5])]

    return OutFile, X

    
#*******************************************************************************
#if __name__ == '__main__':
#    InputParam()
