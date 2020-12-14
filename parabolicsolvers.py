# coding: utf-8
import numpy as np
import preprocess as pre
import postprocess as post
import readuserinputs as ri

#*******************************************************************************
def FTCS(InputSettings, StateOpt=1, totTime=360.0, BCType='Dirichlet'):
    '''Solve a 2D linear parabolic partial differential equation using the
       explicit Forward Time/Central Space method.
       
       Equation solved: UNSTEADY HEAT CONDUCTION EQUATION
    
    Call signature:

        FTCS(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
                   HistFileName, FrameOpt, FrameWrite = InputSettings
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    from globalmod import iMax, jMax, Length, Height

    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
    HistFileName, FrameOpt, FrameWrite = InputSettings
    
    # Read CFL number and diffusion coefficient from file.
    CFL, diff = ri.Coefficients()

    # Calculate step sizes.
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0.
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions.
    U = pre.BC(U, iMax, jMax, BCType)
    
    # Calculate time step from CFL number, alpha and grid step.
    dT = CFL*(1.0/((1/dX**2) + (1/dY**2)))/diff
    
    diffX = diff*dT/dX**2 # diffusion number for x term.
    diffY = diff*dT/dY**2 # diffusion number for y term.

    # Calculate number of time steps, transient solution is desired.
    if StateOpt == 2:
        nMax = int(totTime/dT) + 1

    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        Uold = U.copy()
        
        U[1:-1,1:-1] = U[1:-1,1:-1]\
                       + diffX*(U[2:,1:-1] - 2*U[1:-1,1:-1] + U[0:-2,1:-1])\
                       + diffY*(U[1:-1,2:] - 2*U[1:-1,1:-1] + U[1:-1,0:-2]) # FTCS method
        Error = abs(Uold[1:-1,1:-1] - U[1:-1,1:-1]).sum() # absolute error
        
        # L-inf norm error calculation
        
        #Abs_err = abs(Uold[1:,1:] - U[1:,1:]).sum(axis=1)
        #Error = max(Abs_err)
   
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations."
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using Forward Time/Central Spacing method.'
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
    
    return U

#*******************************************************************************
def DuFortFrankel(InputSettings, StateOpt=1, totTime=360.0, BCType='Dirichlet'):
    '''Solve a 2D linear parabolic partial differential equation using the
       explicit DuFort-Frankel method.
       
       Equation solved: UNSTEADY HEAT CONDUCTION EQUATION
    
    Call signature:

        DuFortFrankel(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
                   HistFileName, FrameOpt, FrameWrite = InputSettings
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    from globalmod import iMax, jMax, Length, Height
    
    ExpNumber, ConvCriteria, nMax, nWrite, OutFileName, nDisplay,\
    HistFileName, FrameOpt, FrameWrite = InputSettings
    
    # Read CFL number and diffusion coefficient from file.
    CFL, diff = ri.Coefficients()

    # Calculate step sizes.
    dX = Length/(iMax - 1)
    dY = Height/(jMax - 1)

    # Initialize fields everywhere at t = 0.
    U = pre.Initial(iMax, jMax)
    
    # Assign boundary conditions.
    U = pre.BC(U, iMax, jMax, BCType)
    Uold = U.copy()
    
    # Calculate time step from CFL number, alpha and grid step.
    dT = CFL*(1.0/((1/dX**2) + (1/dY**2)))/diff
    
    diffX = diff*dT/dX**2 # diffusion number for x term.
    diffY = diff*dT/dY**2 # diffusion number for y term.
    
    # Calculate number of time steps, transient solution is desired.
    if StateOpt == 2:
        nMax = int(totTime/dT) + 1

    n = 0

    Error = 1.0
    print(f'{"ITER":>7} {"ERROR":>15}')
    print(f'{"----":>7} {"-----":>15}')
    print()
    while Error > ConvCriteria and n < nMax: # Start iteration
        Error = 0.0
        n = n + 1
        
        Uold2 = Uold.copy()
        Uold  = U.copy()
        
        if n == 1: # Determine U using FTCS method at n = 1 using n = 0
            U[1:-1,1:-1] = U[1:-1,1:-1]\
                           + diffX*(U[2:,1:-1] - 2*U[1:-1,1:-1] + U[0:-2,1:-1])\
                           + diffY*(U[1:-1,2:] - 2*U[1:-1,1:-1] + U[1:-1,0:-2])
        
        U[1:-1,1:-1] = ((1.0 - 2*diffX - 2*diffY)*Uold2[1:-1,1:-1]\
                       + 2*diffX*(U[2:,1:-1] + U[0:-2,1:-1])\
                       + 2*diffY*(U[1:-1,2:] + U[1:-1,0:-2]))\
                       / (1 + 2*diffX + 2*diffY) # DuFort-Frankel method
        
        Error = abs(Uold[1:-1,1:-1] - U[1:-1,1:-1]).sum() # absolute error
        
        
        post.WriteConvHistToFile(HistFileName, n, Error) # Write convergence history log to a file
        
        U = pre.BC(U, iMax, jMax, BCType) # Update boundary conditions
        
        if n % nDisplay == 0: # Display convergence monitors
            print(f'{n:>7} {Error:>15.8f}')

        # Write output to files
        if n % nWrite == 0:
            post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
            
    flag = 'YES' # Flag to write end of file convergence report 
    if n == nMax:
        msg = f"Solution didn't converge in {nMax} iterations."
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)
    else:
        msg = f'Solution Converged in {n} Iterations with a Maximum Error of\n{ConvCriteria}\
 using explicit DuFort-Frankel method.'
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
    
    return U

#*******************************************************************************
def ADI(InputSettings, BCType='Dirichlet'):
    '''Solve a 2D linear parabolic partial differential equation using the
       implicit Alternating Direction Implicit method.
       
       Equation solved: UNSTEADY HEAT CONDUCTION EQUATION
    
    Call signature:

        ADI(InputSettings, BCType)

    Parameters
    __________

    InputSettings: list
                   
                   The inputs are packed in this list. This is designed to be
                   obtained from function "preprocess.InputSettings(InFileName)"
                   packed into a tuple.

                   This tuple is unpacked as
                   ExpNumber, iMax, jMax, L, H, ConvCriteria, nMax,\
                   nWrite, OutFileName, nDisplay, HistFileName, FrameOpt,\
                   FrameWrite = InputSettings
    
    BCType: {'Dirichlet', 'Neumann', 'Mixed'}, default: 'Dirichlet'
            
            Boundary conditions type for the simulation. Defaults to 'Dirichlet'.
            Possible values:
            'Dirichlet': Dirichlet Boundary Conditions
            'Neumann': Neumann Boundary Conditions
            'Mixed': Mixed-type Boundary Conditions
    '''
    
    ExpNumber, iMax, jMax, L, H, ConvCriteria, nMax,\
            nWrite, OutFileName, nDisplay, HistFileName, FrameOpt,\
            FrameWrite = InputSettings
    
    # Read CFL number and diffusion coefficient from file
    CFL, diff = ri.Coefficients()

    dX = L/(iMax - 1)
    dY = H/(jMax - 1)

    # Initialize fields everywhere at t = 0
    U = pre.Initial(iMax, jMax)
    Uhalf = pre.Initial(iMax, jMax) # Uhalf is U at time level (n + 1/2)
    
    # Assign boundary conditions 
    U = pre.BC(U, iMax, jMax, BCType)
    Uold = U.copy()
    Uhalf = pre.BC(Uhalf, iMax, jMax, BCType)
    
    # Initialize coefficient H and G arrays for solution of
    # tridiagonal system of equations
    ijMax = max(iMax, jMax)
    H = np.zeros(ijMax)
    G = np.zeros(ijMax)
    
    # Calculate time step from CFL number, alpha and grid step.
    dT = CFL*((1/dX**2) + 1/(dY**2))/diff
    
    diffX = diff*dT/dX**2 # diffusion number for x term
    diffY = diff*dT/dY**2 # diffusion number for y term
    
    d1 = 0.5*diffX
    d2 = 0.5*diffY

    # Calculate number of time steps, transient solution is desired.
    if StateOpt == 2:
        nMax = int(totTime/dT) + 1

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
        A = -d1
        B = 1 + 2*d1
        C = -d1
        for j in range (1, jMax-1):
            H[0] = 0.0
            G[0] = U[0][j]
            for i in range (1, iMax-1):
                H[i] = C/(B - A*H[i-1])
                D = d2*U[i][j+1] + (1 - 2*d2)*U[i][j] + d2*U[i][j-1]
                G[i] = (D - A*G[i-1])/(B - A*H[i-1])
            for i in range (iMax-2,0,-1):
                Uhalf[i][j] = -H[i]*Uhalf[i+1][j] + G[i] # Alternating Direction Implicit SOR method in x-direction
                
        Uhalf = pre.BC(Uhalf, iMax, jMax, BCType) # Update boundary conditions for Uhalf
        
        #***********************************************************************
        # This block of codes solves for U at time level n + 1  
        # along constant i line
        # Eq. 5.25 using Tridiagonal system Appendix B in Hoffmann CFD Vol.1
        #***********************************************************************
        A = -d2
        B = 1 + 2*d2
        C = -d2
        for i in range (1, iMax-1):
            H[0] = 0.0
            G[0] = U[i][0]
            for j in range (1, jMax-1):
                H[j] = C/(B - A*H[j-1])
                D = d1*Uhalf[i+1][j] + (1 - 2*d1)*Uhalf[i][j] + d1*Uhalf[i-1][j]
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
 using Alternating Direction Implicit method'
        print(msg)
        post.WriteConvHistToFile(HistFileName, n, Error, WriteFlag=flag, PrintMsg=msg)

    # Write output to files
    post.WriteSolutionToFile(OutFileName, iMax, jMax, dX, dY, U)
