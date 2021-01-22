'''
+**************************************************************************
+**************************************************************************
+
+   FILE         hyperbolicsolvers.py
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
def ExplicitFirstUpwind(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burger equation
    using the explicit first upwind differencing method.

    The wave equation and the Burger equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burger equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        FirstOrderUpwind(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import numpy as np
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        positive_a = 1 + np.sign(init.conv) # for positive a
        negative_a = 1 - np.sign(init.conv) # for negative a
        U[1:-1] = Uo[1:-1]\
                  - 0.5*Courant*positive_a*(Uo[1:-1] - Uo[0:-2])\
                  - 0.5*Courant*negative_a*(Uo[2:] - Uo[1:-1])

    elif init.Model.upper() == 'BURGERS':
        E = Uo*Uo/2
        U[1:-1] = Uo[1:-1] - Courant*(E[1:-1] - E[0:-2])

    return U

#**************************************************************************
def Lax(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burger equation
    using the explicit Lax method.
       
    The wave equation and the Burger equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burger equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        Lax(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        U[1:-1] = 0.5*(Uo[2:] + Uo[0:-2])\
                  - 0.5*Courant*(Uo[2:] - Uo[0:-2])
    
    elif init.Model.upper() == 'BURGERS':
        E = Uo*Uo/2
        U[1:-1] = 0.5*(Uo[2:] + Uo[0:-2])\
                  - 0.25*Courant*(E[2:] - E[0:-2])

    return U

#**************************************************************************
def MidpointLeapfrog(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the explicit Midpoint Leapfrog method.
       
    The wave equation and the Burgers equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        MidpointLeapfrog(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")
    
    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")
        
    #return U

#**************************************************************************
def LaxWendroff(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the explicit Lax-Wendroff method.

    The wave equation and the Burgers equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        LaxWendroff(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        Courant2 = Courant*Courant
        U[1:-1] = Uo[1:-1] - 0.5*Courant*(Uo[2:] - Uo[0:-2]) +\
                  0.5*Courant2*(Uo[2:] - 2.0*Uo[1:-1] + Uo[0:-2])
    
    elif init.Model.upper() == 'BURGERS':
        E = Uo*Uo/2
        Courant2 = Courant*Courant
        U[1:-1] = Uo[1:-1] - 0.5*Courant*(E[2:] - E[0:-2]) +\
                  0.25*Courant2*((Uo[2:] + Uo[1:-1])*(E[2:] - E[1:-1])\
                               - (Uo[1:-1] + Uo[0:-2])*(E[1:-1] - E[0:-2]))

    return U

#**************************************************************************
def LaxWendroffMultiStep(init, Uo, Courant):
    '''Solve a first-order 1D wave equation using the explicit multi-step
    Lax-Wendroff method.
  
    The wave equation is the hyperbolic partial differential equation.
    The first-order wave equation is a linear equation which is expressed
    as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed
                
    Call signature:

        LaxWendroffMultiStep(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    Uhalf = Uo.copy()
    if init.Model.upper() == 'FO_WAVE':
        for i in range (1,init.iMax-1):
            Uhalf[i] = 0.5*(Uo[i+1] + Uo[i]) - 0.5*Courant*(Uo[i+1] - Uo[i])
            U[i] = Uo[i] - Courant*(Uhalf[i] - Uhalf[i-1])
    
    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U

#**************************************************************************
def MacCormack(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the explicit MacCormack method.
   
    The wave equation and the Burgers equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        MacCormack(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    Utemp = Uo.copy()
    if init.Model.upper() == 'FO_WAVE':
        for i in range (1,init.iMax-1):
            Utemp[i] = Uo[i] - Courant*(Uo[i+1] - Uo[i]) # Predictor step
            U[i] = 0.5*((Uo[i] + Utemp[i]) -\
                        Courant*(Utemp[i] - Utemp[i-1])) # Corrector step

    elif init.Model.upper() == 'BURGERS':
        E = Uo*Uo/2
        Etemp=Utemp*Utemp/2
        for i in range (1,init.iMax-1):
            Utemp[i] = Uo[i] - Courant*(E[i+1] - E[i]) # Predictor step
            Etemp[i] = Utemp[i]*Utemp[i]/2
            U[i] = 0.5*((Uo[i] + Utemp[i]) -\
                        Courant*(Etemp[i] - Etemp[i-1])) # Corrector step

    return U

#**************************************************************************
def FourthOrderRungeKutta(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the explicit four-stage Runge-Kutta method.

    The wave equation and the Burgers equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        ModifiedRungeKutta(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    if init.Model.upper() == 'FO_WAVE':
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    U = Uo.copy()# Initialize U
    if init.Model.upper() == 'FO_WAVE':
        U1 = Uo.copy()
        U2 = Uo.copy()
        U3 = Uo.copy()
        # 1st stage
        U1[1:-1] = Uo[1:-1] - 0.5*Courant*(Uo[2:] - Uo[0:-2])/2.0
        # -- update BC here
        # 2nd stage
        U2[1:-1] = Uo[1:-1] - 0.5*Courant*(U1[2:] - U1[0:-2])/2.0
        # -- update BC here
        # 3rd stage
        U3[1:-1] = Uo[1:-1] - Courant*(U2[2:] - U2[0:-2])/2.0
        # -- update BC here
        # 4th stage
        U[1:-1] = Uo[1:-1] - 0.5*Courant*\
                  (((1.0/6)*(Uo[2:] - Uo[0:-2])) +\
                   ((1.0/3)*(U1[2:] - U1[0:-2])) +\
                   ((1.0/3)*(U2[2:] - U2[0:-2])) +\
                   ((1.0/6)*(U3[2:] - U3[0:-2])))
        # -- update BC here
    
    elif init.Model.upper() == 'BURGERS':
        # 1st stage
        U1 = Uo.copy()
        U2 = Uo.copy()
        U3 = Uo.copy()
        E1 = Uo*Uo/2
        U1[1:-1] = Uo[1:-1] - 0.5*Courant*(E1[2:] - E1[0:-2])/2.0
        # -- update BC here
        # 2nd stage
        E2 = U1*U1/2
        U2[1:-1] = Uo[1:-1] - 0.5*Courant*(E2[2:] - E2[0:-2])/2.0
        # -- update BC here
        # 3rd stage
        E3 = U2*U2/2
        U3[1:-1] = Uo[1:-1] - Courant*(E3[2:] - E3[0:-2])/2.0
        # -- update BC here
        # 4th stage
        E4 = U3*U3/2
        U[1:-1] = Uo[1:-1] - 0.5*Courant*\
                  (((1.0/6)*(E1[2:] - E1[0:-2])) +\
                   ((1.0/3)*(E2[2:] - E2[0:-2])) +\
                   ((1.0/3)*(E3[2:] - E3[0:-2])) +\
                   ((1.0/6)*(E4[2:] - E4[0:-2])))
        # -- update BC here

    '''if TVD == 'Yes':
        phiPlus = HartenYeeTVD()
        phiMinus = HartenYeeTVD()
        U = U - 0.5*convX*(phiPlus - phiMinus)'''
    
    return U

#**************************************************************************
def ModifiedRungeKutta(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the explicit four-stage Modified Runge-Kutta method.

    The wave equation and the Burgers equation are the hyperbolic partial
    differential equation. The first-order wave equation is a linear
    equation which is expressed as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                E = u^2/2
                
    Call signature:

        ModifiedRungeKutta(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        # 1st stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:] - U[0:-2])/8.0
        # -- update BC here
        # 2nd stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:] - U[0:-2])/6.0
        # -- update BC here
        # 3rd stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:] - U[0:-2])/4.0
        # -- update BC here
        # 4th stage
        U[1:-1] = Uo[1:-1] - Courant*(U[2:] - U[0:-2])/2.0
        # -- update BC here
    
    elif init.Model.upper() == 'BURGERS':
        # 1st stage
        E = Uo*Uo/2
        U[1:-1] = Uo[1:-1] - Courant*(E[2:] - E[0:-2])/8.0
        # -- update BC here
        # 2nd stage
        E = U*U/2
        U[1:-1] = Uo[1:-1] - Courant*(E[2:] - E[0:-2])/6.0
        # -- update BC here
        # 3rd stage
        E = U*U/2
        U[1:-1] = Uo[1:-1] - Courant*(E[2:] - E[0:-2])/4.0
        # -- update BC here
        # 4th stage
        E = U*U/2
        U[1:-1] = Uo[1:-1] - Courant*(E[2:] - E[0:-2])/2.0
        # -- update BC here
    
    return U

#**************************************************************************
def EulersBTCS(init, Uo, Courant):
    '''Solve a first-order 1D wave equation using the implicit Euler's
    Backward Time Central Space (BTCS) method.
    
    The wave equation is the hyperbolic partial differential equation.
    The first-order wave equation is a linear equation which is expressed
    as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    Call signature:

        EulersBTCS(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import fluid.tridiagonal as trid

    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        cc = 0.5*Courant
        A = [cc for i in range(init.iMax)]
        B = [-1 for i in range(init.iMax)]
        C = [-cc for i in range(init.iMax)]
        D = -Uo
        UU = Uo.copy()

        U = trid.TridiagonalSolver(init.iMax,A, B, C, D, UU)

    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U

#**************************************************************************
def CrankNicolson(init, Uo, Courant):
    '''Solve a first-order 1D wave equation using the implicit
    Crank-Nicolson method.

    The wave equation is the hyperbolic partial differential equation.
    The first-order wave equation is a linear equation which is expressed
    as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    Call signature:

        CrankNicolson(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import fluid.tridiagonal as trid

    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        cc = 0.25*Courant
        A = [cc for i in range(init.iMax)]
        B = [-1.0 for i in range(init.iMax)]
        C = [-cc for i in range(init.iMax)]
        D = [0 for i in range (init.iMax)]
        D[1:-1] = -Uo[1:-1] + 0.25*Courant*(Uo[2:] - Uo[0:-2])
        UU = Uo.copy()

        U = trid.TridiagonalSolver(init.iMax, A, B, C, D, UU)

    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U

#**************************************************************************
def BeamAndWarming(init, Uo, Courant):
    '''Solve a first-order 1D wave equation or inviscid Burgers equation
    using the implicit Beam and Warming method.

    The Burgers equation is the hyperbolic partial differential equation.
    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                u: measurable quanity
                E = u^2/2
                
    Call signature:

        BeamAndWarming(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import fluid.tridiagonal as trid

    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    if init.Model.upper() == 'FO_WAVE':
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    elif init.Model.upper() == 'BURGERS':
        E = Uo*Uo/2
        cc = 0.25*Courant
        A = [0.0 for i in range(init.iMax)]
        B = [1.0 for i in range(init.iMax)]
        C = [0.0 for i in range(init.iMax)]
        D = [0.0 for i in range(init.iMax)]
        A[1:-1] = -cc*Uo[0:-2]
        C[1:-1] = cc*Uo[2:]
        D[1:-1] = Uo[1:-1] - 0.5*Courant*(E[2:] - E[0:-2]) +\
                  0.25*Courant*(Uo[2:]*Uo[2:] - Uo[0:-2]*Uo[0:-2])
        UU = Uo.copy()

        U = trid.TridiagonalSolver(init.iMax, A, B, C, D, UU)

    return U

#**************************************************************************
def FirstOrderTVD(init, Uo, Courant):
    '''Solve a first-order inviscid Burgers equation using the second-
    order TVD schemes and their various Limiter Functions and Limiters.

    The Burgers equation is the hyperbolic partial differential equation.
    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                u: measurable quanity
                E = u^2/2
                
    Call signature:

        FirstOrderTVD(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import fluid.secondaryfunctions as sec

    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    if init.Model.upper() == 'FO_WAVE':
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    U = Uo.copy() # Initialize U
    #Utemp = Uo.copy()
    E = Uo*Uo/2

    for i in range (1,init.iMax-1):
        dUiPlus12 = sec.CalcUi(Uo[i+1], Uo[i])
        dUiMinus12 = sec.CalcUi(Uo[i], Uo[i-1])
        #Ei = sec.CalcE(Uo[i])
        #EiPlus1 = sec.CalcE(Uo[i+1])
        #EiMinus1 = sec.CalcE(Uo[i-1])
        # -- Calculate alpha(i+1/2) using equation 6-98
        if dUiPlus12 == 0:
            alphaiPlus12 = Uo[i]
        else:
            alphaiPlus12 = (E[i+1] - E[i])/dUiPlus12
        # -- Calculate alpha(i-1/2) using equation 6-100
        if dUiMinus12 == 0:
            alphaiMinus12 = Uo[i]
        else:
            alphaiMinus12 = (E[i] - E[i-1])/dUiMinus12
        # Equation 6-119 and 6-120 in CFD Vol. 1 by Hoffmann
        phiPlus = abs(alphaiPlus12)*dUiPlus12
        phiMinus = abs(alphaiMinus12)*dUiMinus12
        # Equation 6-117 and 6-118 in CFD Vol. 1 by Hoffmann
        Utemp = Uo[i] - 0.5*Courant*(E[i+1] - E[i-1])
        U[i] = Utemp + 0.5*Courant*(phiPlus - phiMinus)

    return U

#**************************************************************************
def SecondOrderTVD(init, Uo, Courant, LimiterFunc, Limiter, Eps=0.1):
    '''Solve a first-order inviscid Burgers equation using the second-
    order TVD schemes and their various Limiter Functions and Limiters.

    The Burgers equation is the hyperbolic partial differential equation.
    The first-order inviscid Burgers equation is a non-linear equation
    which is expressed as:

                            du/dt = -u(du/dx)   or,

                            du/dt = -dE/dx
    where,
                u: measurable quanity
                E = u^2/2
                
    Call signature:

        SecondOrderTVD(init, Uo, Courant)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
    import fluid.tvdfunctions as tvd
    import preprocess.fetchoptions as opt

    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    if init.Model.upper() == 'FO_WAVE':
        raise Exception("This formulation is not available for WAVE\
 equation in this version.")

    U = Uo.copy() # Initialize U
    E = Uo*Uo/2

    fetch = opt.FetchOptions()
    limfunc_options = fetch.TVDLimiterFunctionOptions()

    if not LimiterFunc in limfunc_options:
        raise Exception(f"Invalid flux limiter function selection in the call\
 to function\nSecondOrderTVD().\Valid options for LimiterFunc are:\
 {limfunc_options}.")
            
    for i in range (2,init.iMax-2):
        phiPlus, phiMinus = tvd.CalculateTVD(i, Uo, E, Eps, Courant,\
                                             Limiter, LimiterFunc)

        # Equation 6-124 and 6-125 in Hoffmann Vol. 1
        hPlus = 0.5*(E[i+1] + E[i] + phiPlus)
        hMinus = 0.5*(E[i] + E[i-1] + phiMinus)

        # Equation 6-123
        U[i] = Uo[i] - Courant*(hPlus - hMinus)

    return U

#**************************************************************************
#def FluxCorrectedTransport(init, Uo, Courant, Damp1, Damp2):
    '''Solve a first-order 1D wave equation using the Flux Corrected
    Transport scheme for the Lax-Wendroff method.
   
    The wave equation is the hyperbolic partial differential equation.
    The first-order wave equation is a linear equation which is expressed
    as:

                            du/dt = -a(du/dx)       for a>0
    where,
                u: measurable quanity
                a: constant speed

    Call signature:

        FluxCorrectedTransport(init, Uo, Courant, Damp1, Damp2)

    Parameters
    ----------

    init :

           Class object of RunConfig class which was created at the
           beginning of the simulation.

    Uo : 1D array

         The dependent variable from time level (n) within the domain.

    Courant : float

              Courant number (entered as user input in file).

    Damp1 : float

            Damping term which is added to the predictor step.

    Damp2 : float

            Antii-diffusive term which is added to the corrector step
            to remove excessive damping.

    Returns
    -------

    U : 1D array

        The dependent variable calculated at time level (n+1) within the
        entire domain.
    '''
'''    import fluid.secondaryfunctions as sf
    shapeU = Uo.shape # Obtain Dimension
    if len(shapeU) == 2:
        raise Exception("This formulation is only available for 1D first\
 order wave equation or the inviscid Burgers equation in this version.")

    U = Uo.copy() # Initialize U
    Utemp = Uo.copy()
    if init.Model.upper() == 'FO_WAVE':
        Courant2 = Courant*Courant
        for i in range (1,init.iMax-1):
            # Predictor Step
            Utemp[i] = Uo[i]\
                       - 0.5*Courant*(Uo[i+1] - Uo[i-1])\
                       + (Damp1 + 0.5*Courant2)*\
                           (Uo[i+1] - 2.0*Uo[i] + Uo[i-1])
        for i in range (2,init.iMax-2):
            # Corrector step
            U[i] = Utemp[i]\
                   - Damp2*(Utemp[i+1] - 2.0*Utemp[i] + Utemp[i-1])

    elif init.Model.upper() == 'BURGERS':
        raise Exception("This formulation is not available for BURGERS\
 equation in this version.")

    return U'''
#**************************************************************************
