# coding: utf-8
import preprocess as pre
import ellipticsolvers as ep
from colorama import Fore

filename = ['Inp_PGS.dat', 'Inp_LGSi.dat', 'Inp_LGSj.dat', 'Inp_PSOR.dat', 'Inp_LSORi.dat', 'Inp_LSORj.dat']
for i in range (0,6):
    Inset = pre.InputParam(InFileName = './Input/'+filename[i])
    
    if i == 0:
        print(Fore.YELLOW + '')
        print('************************************')
        print('     POINT GAUSS-SEIDEL METHOD')
        print('************************************')
        ep.PointGaussSeidel(InputSettings = Inset)
    elif i == 1:
        print(Fore.RED + '')
        print('************************************')
        print(' LINE GAUSS-SEIDEL METHOD ALONG Y')
        print('************************************')
        ep.LineGaussSeidel_i(InputSettings = Inset)
    elif i == 2:
        print(Fore.GREEN + '')
        print('************************************')
        print(' LINE GAUSS-SEIDEL METHOD ALONG X')
        print('************************************')
        ep.LineGaussSeidel_j(InputSettings = Inset)
    elif i == 3:
        print(Fore.BLUE + '')
        print('************************************')
        print('          POINT SOR METHOD')
        print('************************************')
        ep.PSOR(InputSettings = Inset, RelaxParam = 1.78)
    elif i == 4:
        print(Fore.CYAN + '')
        print('************************************')
        print('      LINE SOR METHOD ALONG Y')
        print('************************************')
        ep.LSOR_i(InputSettings = Inset, RelaxParam = 1.265)
    elif i == 5:
        print(Fore.MAGENTA + '')
        print('************************************')
        print('      LINE SOR METHOD ALONG X')
        print('************************************')
        ep.LSOR_j(InputSettings = Inset, RelaxParam = 1.265)
