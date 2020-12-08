# coding: utf-8

import math
import numpy as np
import preprocess as pre
import postprocess as post
import readuserinputs as ri

#*******************************************************************************
def PointGaussSeidel(InputSettings, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Point-Gauss
       Seidel method

    Call signature:

        PointGaussSeidel(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list

                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings

    BCType: {'Dirichlet', 'Neumann', 'Mixed'}
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''

    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height

    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)

    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)

    Beta = dX/dY
    B2 = Beta**2.0
    A = 0.5/(1 + B2)
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax:
        Error = 0.0
        n = n + 1
        
        for i in range (1, iMax-1):
            for j in range (1, jMax-1):
                Uold = U[i][j]
                U[i][j] = A*(U[i+1][j] + U[i-1][j] + B2*(U[i][j+1] + U[i][j-1])) # Point Gauss_Seidel method
                Error = Error + abs(Uold - U[i][j]) # Calculate error

        # Write convergence history log to a file
        post.WriteConvHistToFile(HistFileName, n, Error)
   

        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')
        
        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    flag = 'YES' # Flag to write end of file convergence report
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Point Gauss Seidel method.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def LineGaussSeidel_i(InputSettings, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Line-Gauss
       Seidel method along constant i direction (parallel to y-axis)    

    Call signature:

        LineGaussSeidel_i(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    H = np.zeros(jMax)
    G = np.zeros(jMax)
    
    Beta = dX/dY
    A = Beta**2.0
    B = -2.0*(1.0 + A)
    C = Beta**2.0
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        for i in range (1, iMax-1):
            H[0] = 0.0
            G[0] = U[i][0]
            for j in range (1, jMax-1):
                H[j] = C/(B - A*H[j-1])
                D = -(U[i+1][j] + U[i-1][j])
                G[j] = (D - A*G[j-1])/(B - A*H[j-1])
            for j in range (jMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[j]*U[i][j+1] + G[j] # Line Gauss_Seidel method along constant i
                Error = Error + abs(Uold - U[i][j]) # Calculate error

        # Write convergence history log to a file
        post.WriteConvHistToFile(HistFileName, n, Error)
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    flag = 'YES' # Flag to write end of file convergence report
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Line Gauss Seidel method along\nconstant y-direction.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def LineGaussSeidel_j(InputSettings, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Line-Gauss
       Seidel method along constant j direction (parallel to x-axis)
    
    Call signature:

        LineGaussSeidel_j(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height

    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    H = np.zeros(iMax)
    G = np.zeros(iMax)
    
    Beta = dX/dY
    A = Beta**2.0
    B = -2.0*(1.0 + A)
    C = Beta**2.0
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        for j in range (1, jMax-1):
            H[0] = 0.0
            G[0] = U[0][j]
            for i in range (1, iMax-1):
                H[i] = C/(B - A*H[i-1])
                D = -(U[i][j+1] + U[i][j-1])
                G[i] = (D - A*G[i-1])/(B - A*H[i-1])
            for i in range (iMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[i]*U[i+1][j] + G[i] # Line Gauss_Seidel method along constant j 
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        # Write convergence history log to a file
        post.WriteConvHistToFile(HistFileName, n, Error)
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    flag = 'YES' # Flag to write end of file convergence report
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Line Gauss Seidel method along\nconstant x-direction.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def PSOR(InputSettings, RelaxParam, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Point
       Successive Over-Relaxation (PSOR) method
    
    Call signature:

        PSOR(InputSettings, RelaxParam, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    RelaxParam: float
    
                Relaxation Parameter is used for faster convergence of PSOR method.
                Specify RelaxParam values between 0 and 2 to obtain convergence.
                If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
                If RelaxParam = 1: Point Gauss-Seidel method is recovered.
                An optimum value is determined by performing numerical experimentations.
                RelaxParam = 1.78 was found to be an optimum value for PSOR method for
                the problem with a rectangular domain having uniform grid step with the
                Dirichlet BC imposed (see Hoffmann Vol. 1, pg 164, 170).
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    
    Beta = dX/dY
    B2 = Beta**2.0
    A = 0.5/(1.0 + B2)
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        for i in range (1, iMax-1):
            for j in range (1, jMax-1):
                Uold = U[i][j]
                U[i][j] = (1 - RelaxParam)*U[i][j] + RelaxParam*\
                          A*(U[i+1][j] + U[i-1][j] +\
                          B2*(U[i][j+1] + U[i][j-1])) # Point Successive Over Relaxation method
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Point Successive Over-Relaxation method\nwith Relaxation Parameter = {RelaxParam}.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def LSOR_i(InputSettings, RelaxParam, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Line
       Successive-Over Relaxation (LSOR) method along constant i direction
       (parallel to y-axis)
    
    Call signature:

        LSOR_i(InputSettings, RelaxParam, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    RelaxParam: float
    
                Relaxation Parameter is used for faster convergence of LSOR method.
                Specify RelaxParam values between 0 and 2 to obtain convergence.
                If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
                If RelaxParam = 1: Line Gauss-Seidel method is recovered.
                An optimum value is determined by performing numerical experimentations.
                RelaxParam = 1.265 was found to be an optimum value for LSOR_i method 
                for the problem with a rectangular domain having uniform grid step with
                the Dirichlet BC imposed (see Hoffmann Vol. 1, pg 165, 171).
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    H = np.zeros(jMax)
    G = np.zeros(jMax)
    
    Beta = dX/dY
    B2 = Beta**2.0
    A = RelaxParam*B2
    B = -2.0*(1+B2)
    C = RelaxParam*B2
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        for i in range (1, iMax-1):
            H[0] = 0.0
            G[0] = U[i][0]
            for j in range (1, jMax-1):
                H[j] = C/(B - A*H[j-1])
                D = (1 - RelaxParam)*B*U[i][j] -\
                    RelaxParam*(U[i+1][j] + U[i-1][j])
                G[j] = (D - A*G[j-1])/(B - A*H[j-1])
            for j in range (jMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[j]*U[i][j+1] + G[j] # Line Successive Over Relaxation method along constant i
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Line Successive Over-Relaxation method along\nconstant y-direction with\
 Relaxation Parameter = {RelaxParam}.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def LSOR_j(InputSettings, RelaxParam, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Line
       Successive Over-Relaxation (LSOR) method along constant j direction
       (parallel to x-axis)
    
    Call signature:

        LSOR_j(InputSettings, RelaxParam, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    RelaxParam: float
    
                Relaxation Parameter is used for faster convergence of LSOR method.
                Specify RelaxParam values between 0 and 2 to obtain convergence.
                If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
                If RelaxParam = 1: Line Gauss-Seidel method is recovered.
                An optimum value is determined by performing numerical experimentations.
                RelaxParam = 1.265 was found to be an optimum value for LSOR_j method 
                for the problem with a rectangular domain having uniform grid step with
                the Dirichlet BC imposed (see Hoffmann Vol. 1, pg 165, 171).
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    H = np.zeros(iMax)
    G = np.zeros(iMax)
    
    Beta = dX/dY
    B2 = Beta**2.0
    A = RelaxParam
    B = -2.0*(1+B2)
    C = RelaxParam
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        for j in range (1, jMax-1):
            H[0] = 0.0
            G[0] = U[0][j]
            for i in range (1, iMax-1):
                H[i] = C/(B - A*H[i-1])
                D = (1 - RelaxParam)*B*U[i][j] -\
                    RelaxParam*B2*(U[i][j+1] + U[i][j-1])
                G[i] = (D - A*G[i-1])/(B - A*H[i-1])
            for i in range (iMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[i]*U[i+1][j] + G[i] # Line Successive Over Relaxation method along constant j
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Line Successive Over-Relaxation method along\nconstant x-direction with\
 Relaxation Parameter = {RelaxParam}.'
        print(msg)
        print()
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def ADI(InputSettings, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Alternating
       Direction Implicit (ADI) method.
    
    Call signature:

        ADI(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    Uhalf = pre.Initial (iMax, jMax) # Uhalf is U at time level (n + 1/2)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    Uhalf = pre.BC(Uhalf, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    ijMax = max(iMax, jMax)
    H = np.zeros(ijMax)
    G = np.zeros(ijMax)
    
    Beta = dX/dY
    B2 = Beta**2.0
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        #***********************************************************************
        # This block of codes solves for U at time level n + 1/2 (i.e. = Uhalf) 
        # along constant j line 
        # Eq. 5.22 using Tridiagonal system Appendix B in Hoffmann CFD Vol.1
        #***********************************************************************
        A = 1.0
        B = -2.0*(1 + B2)
        C = 1.0
        for j in range (1, jMax-1):
            H[0] = 0.0
            G[0] = U[0][j]
            for i in range (1, iMax-1):
                H[i] = C/(B - A*H[i-1])
                D = -B2*(U[i][j+1] + Uhalf[i][j-1])
                G[i] = (D - A*G[i-1])/(B - A*H[i-1])
            for i in range (iMax-2,0,-1):
                Uhalf[i][j] = -H[i]*Uhalf[i+1][j] + G[i] # Alternating Direction Implicit method in x-direction
                
        Uhalf = pre.BC(Uhalf, iMax, jMax, BCType) # Update boundary conditions for Uhalf
        
        #***********************************************************************
        # This block of codes solves for U at time level n + 1  
        # along constant i line
        # Eq. 5.23 using Tridiagonal system Appendix B in Hoffmann CFD Vol.1
        #***********************************************************************
        A = B2
        B = -2.0*(1 + B2)
        C = B2
        for i in range (1, iMax-1):
            H[0] = 0.0
            G[0] = U[i][0]
            for j in range (1, jMax-1):
                H[j] = C/(B - A*H[j-1])
                D = -(Uhalf[i+1][j] + U[i-1][j])
                G[j] = (D - A*G[j-1])/(B - A*H[j-1])
            for j in range (jMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[j]*U[i][j+1] + G[j] # Alternating Direction Implicit method in y-direction
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Alternating Direction Implicit method'
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
def ADISOR(InputSettings, RelaxParam, BCType='Dirichlet'):
    '''Solve a 2D elliptic partial differential equation using the Alternating
       Direction Implicit Successive Over-Relaxation method.
    
    Call signature:

        ADI(InputSettings, RelaxParam, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax,, nWrite, OutFileName,\
                   nDisplay, HistFileName, FrameOpt, FrameWrite = InputSettings
                   
    RelaxParam: float
    
                Relaxation Parameter is used for faster convergence of ADISOR method.
                Specify RelaxParam values between 0 and 2 to obtain convergence.
                If 0 < RelaxParam < 1: it is called UNDER-RELAXATION.
                If RelaxParam = 1: Alternating Direction Implicit method is recovered.
                An optimum value is determined by performing numerical experimentations.
                RelaxParam = 1.27 was found to be an optimum value for ADISOR method 
                for the problem with a rectangular domain having uniform grid step with
                the Dirichlet BC imposed (see Hoffmann Vol. 1, pg 172, 183).
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
               HistFileName, FrameOpt, FrameWrite = InputSettings

    iMax = ri.iMax
    jMax = ri.jMax
    Length = ri.Length
    Height = ri.Height
    
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    Uhalf = pre.Initial (iMax, jMax) # Uhalf is U at time level (n + 1/2)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    Uhalf = pre.BC(Uhalf, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    ijMax = max(iMax, jMax)
    H = np.zeros(ijMax)
    G = np.zeros(ijMax)
    
    Beta = dX/dY
    B2 = Beta**2.0
    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1

        #***********************************************************************
        # This block of codes solves for U at time level n + 1/2 (i.e. = Uhalf) 
        # along constant j line 
        # Eq. 5.24 using Tridiagonal system Appendix B in Hoffmann CFD Vol.1
        #***********************************************************************
        A = RelaxParam
        B = -2.0*(1 + B2)
        C = RelaxParam
        for j in range (1, jMax-1):
            H[0] = 0.0
            G[0] = U[0][j]
            for i in range (1, iMax-1):
                H[i] = C/(B - A*H[i-1])
                D = (1 - RelaxParam)*B*U[i][j] - RelaxParam*B2*(U[i][j+1] + Uhalf[i][j-1])
                G[i] = (D - A*G[i-1])/(B - A*H[i-1])
            for i in range (iMax-2,0,-1):
                Uhalf[i][j] = -H[i]*Uhalf[i+1][j] + G[i] # Alternating Direction Implicit SOR method in x-direction
                
        Uhalf = pre.BC(Uhalf, iMax, jMax, BCType) # Update boundary conditions for Uhalf
        
        #***********************************************************************
        # This block of codes solves for U at time level n + 1  
        # along constant i line
        # Eq. 5.25 using Tridiagonal system Appendix B in Hoffmann CFD Vol.1
        #***********************************************************************
        A = RelaxParam*B2
        B = -2.0*(1 + B2)
        C = RelaxParam*B2
        for i in range (1, iMax-1):
            H[0] = 0.0
            G[0] = U[i][0]
            for j in range (1, jMax-1):
                H[j] = C/(B - A*H[j-1])
                D = (1 - RelaxParam)*B*Uhalf[i][j] - RelaxParam*(Uhalf[i+1][j] + U[i-1][j])
                G[j] = (D - A*G[j-1])/(B - A*H[j-1]) 
            for j in range (jMax-2,0,-1):
                Uold = U[i][j]
                U[i][j] = -H[j]*U[i][j+1] + G[j] # Alternating Direction Implicit SOR method in y-direction
                Error = Error + abs(Uold - U[i][j]) # Calculate error
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations"
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Alternating Direction Implicit SOR method with\nRelaxation Parameter = {RelaxParam} '
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)

    return U

#*******************************************************************************
