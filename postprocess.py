# coding: utf-8
import numpy as np

#***********************************************************************************
def WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U):
    '''Function to write simulation output to files'''
    
    OutFile = open(OutFileName, "w")
    
    for i in range (0, iMax):
        for j in range (0, jMax):
            print(f'{dX*i:>12.8f} {dY*j:>12.8f} {U[i][j]:>12.8f}', file=OutFile)
    
    OutFile.close()

#***********************************************************************************
def WriteConvHistToFile(HistFileName, n, Error, WriteFlag='NO', PrintMsg='NONE'):
    '''Function to write convergence history log'''
    
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
