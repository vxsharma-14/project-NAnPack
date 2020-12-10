# coding: utf-8
import readuserinputs as ri
import postprocess as post
import benchmark as bm
from importlib import reload
import globalmod
import ellipticsolvers as ep
from colorama import Fore

globalmod.iMax = 21
globalmod.jMax = 41
globalmod.Length = 1.0
globalmod.Height = 2.0

print('Calculate analytical solution for validating numerical solvers.')
Ta = bm.HeatConduction(Length=1.0, Height=2.0, iMax=21, jMax=41, T1=100.0, T2=0.0, T3=0.0, T4=0.0)
OutF, arr = ri.Read1DNodesFile()
print(arr)
post.Write1DSolutionToFile(OutF, arr, Ta, 'Y')
print('Analytical calculation complete.')

print('Start numerical solutions.')

OutFile, X = ri.Read1DNodesFile()
filename = ['Inp_PGS.dat', 'Inp_LGSi.dat', 'Inp_LGSj.dat', 'Inp_PSOR.dat', 'Inp_LSORi.dat', 'Inp_LSORj.dat']

print("1. Point Gauss-Seidel Method")
print("2. Line Gauss-Seidel Method along Y")
print("3. Line Gauss-Seidel Method along X")
print("4. Point SOR Method")
print("5. Line SOR Method along Y")
print("6. Line SOR Method along X")
print("7. ADI SOR Method")
i = int(input("Enter choice"))
flag = 'YES'

Inset = ri.InputParam(InFileName = './Input/'+filename[i])

while flag == 'YES':

    if i == 1:
        print(Fore.YELLOW + '')
        print('************************************')
        print('     POINT GAUSS-SEIDEL METHOD')
        print('************************************')
        U = ep.PointGaussSeidel(InputSettings = Inset)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 2:
        print(Fore.RED + '')
        print('************************************')
        print(' LINE GAUSS-SEIDEL METHOD ALONG Y')
        print('************************************')
        U = ep.LineGaussSeidel_i(InputSettings = Inset)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 3:
        print(Fore.GREEN + '')
        print('************************************')
        print(' LINE GAUSS-SEIDEL METHOD ALONG X')
        print('************************************')
        U = ep.LineGaussSeidel_j(InputSettings = Inset)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 4:
        print(Fore.BLUE + '')
        print('************************************')
        print('          POINT SOR METHOD')
        print('************************************')
        U = ep.PSOR(InputSettings = Inset, RelaxParam = 1.78)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 5:
        print(Fore.CYAN + '')
        print('************************************')
        print('      LINE SOR METHOD ALONG Y')
        print('************************************')
        U = ep.LSOR_i(InputSettings = Inset, RelaxParam = 1.265)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 6:
        print(Fore.MAGENTA + '')
        print('************************************')
        print('      LINE SOR METHOD ALONG X')
        print('************************************')
        U = ep.LSOR_j(InputSettings = Inset, RelaxParam = 1.265)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')

    elif i == 7:
        print('************************************')
        print('      LINE SOR METHOD ALONG X')
        print('************************************')
        U = ep.ADISOR(InputSettings = Inset, RelaxParam = 1.265)
        post.Write1DSolutionToFile(OutFile, X, U, 'Y')
        
    print('Do you want to solve using other scheme? [Y/N]')
    ch = input('Enter Y/N')
    
    while ch != 'Y' and ch != 'N':
        print("Incorrect character, try again")
        print('Do you want to solve using other scheme? [Y/N]')
        ch = input('Enter Y/N')

    if ch == 'N':
        flag = 'NO'
        print('Program Terminating')
        
    elif ch == 'Y':
        print("1. Point Gauss-Seidel Method")
        print("2. Line Gauss-Seidel Method along Y")
        print("3. Line Gauss-Seidel Method along X")
        print("4. Point SOR Method")
        print("5. Line SOR Method along Y")
        print("6. Line SOR Method along X")
        print("7. ADI SOR Method")
        i = int(input("Enter choice"))
