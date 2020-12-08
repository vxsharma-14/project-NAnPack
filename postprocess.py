# coding: utf-8
import numpy as np

#***********************************************************************************
def WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U):
    '''Write simulation results of the entire domain to output file.

    Call Signature:

        WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    Parameters
    __________

    OutFileName: string

                 File path and name of the output file.

    iMax, jMax: int, int

                Grid points in x and y directions, respectively.

    dX, dY: float, float

            Grid step sizes in x and  y direction, respectively.

    U: numpy 2D array, float

       Solution of the dependent variable to be stored.
    '''
    
    OutFile = open(OutFileName, "w")
    
    for i in range (0, iMax):
        for j in range (0, jMax):
            print(f'{dX*i:>12.8f} {dY*j:>12.8f} {U[i][j]:>12.8f}', file=OutFile)
    
    OutFile.close()

#***********************************************************************************
def WriteConvHistToFile(HistFileName, n, Error, WriteFlag='NO', PrintMsg='NONE'):
    '''Write convergence history log

    Call Signature:

        WriteConvHistToFile(HistFileName, n, Error, WriteFlag, PrintMsg)

    Parameters
    __________

    HistFileName: string

                  File path and name of the output file.

    n: int

       Iteration level (time level).

    Error: float

           Error in the solution at each iteration level.

    WriteFlag: boolean (YES/NO), default: 'NO'

               Internal flag to control print message at the end.

    PrintMsg: string, default: 'NONE'

              Internal message stored at the end of the convergence file.
    '''
    
    if WriteFlag == 'NO':
        if n == 1:
            HistFile = open(HistFileName, "w")
            #print('', file=HistFile)
            print(f'{"ITER":>7} {"ERROR":>15}', file=HistFile)
            print(f'{"----":>7} {"-----":>15}', file=HistFile)
        else:
            HistFile = open(HistFileName, "a")

        print(f'{n:>7}{Error:>15.8f}', file=HistFile)

    elif WriteFlag == 'YES':
        HistFile = open(HistFileName, "a")
        #HistFile.seek(0,0)
        print(PrintMsg, file=HistFile)

    HistFile.close()

#***********************************************************************************
def Write1DSolutionToFile(Out1DFileName, iX, U, Direction):
    '''Write output data in 1D format at nodes of X or Y

    Call Signature:

        WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    Parameters
    __________

    Out1DFileName: string

                   File path and name of the 1D output file.

    U: numpy 2D array, float

       Solution of the dependent variable to be stored.

    Direction: string, either 'X' or 'Y'

               Define directions in which 1D results will be stored.
               If Direction = 'X', results will be stored along X-axis at
                                   various Y-nodes (max 5 nodes).
               If Direction = 'Y', results will be stored along Y-axis at
                                   various X-nodes (max 5 nodes).
    '''

    from readuserinputs import iMax, jMax, Length, Height

    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)
    
    OutFile = open(Out1DFileName, "w")
    
    if Direction == 'Y':

        X = [dX*i for i in iX]
        print(f'{"Y":^10} {"X=":>6}{X[0]:<6} {"X=":>6}{X[1]:<6} {"X=":>6}{X[2]:<6}\
{"X=":>6}{X[3]:<6} {"X=":>6}{X[4]:<6}',file=OutFile)
        for j in range(0,jMax):
            print(f'{dY*j:<10.6f} {U[iX[0]][j]:>12.8f} {U[iX[1]][j]:>12.8f}\
{U[iX[2]][j]:>12.8f} {U[iX[3]][j]:>12.8f} {U[iX[4]][j]:>12.8f}', file=OutFile)
            
    elif Direction == 'X':

        Y = [dY*i for i in iX]
        print(f'{"X":^10} {"Y=":>6}{Y[0]:<6} {"Y=":>6}{Y[1]:<6} {"Y=":>6}{Y[2]:<6}\
{"Y=":>6}{Y[3]:<6} {"Y=":>6}{Y[4]:<6}',file=OutFile)
        for i in range(0,iMax):
            print(f'{dX*i:<10.6f} {U[i][iX[0]]:>12.8f} {U[i][iX[1]]:>12.8f}\
{U[i][iX[2]]:>12.8f} {U[i][iX[3]]:>12.8f} {U[i][iX[4]]:>12.8f}', file=OutFile)
            
    else:
        print("Check Direction input again, use either 'X' or 'Y' (case-sensitive)")

    OutFile.close()
