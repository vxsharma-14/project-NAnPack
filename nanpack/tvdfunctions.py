''''
+**************************************************************************
+**************************************************************************
+
+   FILE         tvdfunctions.py
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
def CalculateTVD(i, Uo, E, Eps, Courant, Limiter, LimFunction):
    ''' A dictionary to call a function based on the required TVD scheme
    which is given as input to this function.
    '''
    tvd = {
        "Harten-Yee-Upwind" : HartenYeeUpwind,
        "Modified-Harten-Yee-Upwind" : ModifiedHartenYeeUpwind,
        "Roe-Sweby-Upwind" : RoeSwebyUpwind,
        "Davis-Yee-Symmetric" : DavisYeeSymmetric
        }

    SelectedFunction = tvd.get(LimFunction)
    return SelectedFunction(i,Uo,E,Eps,Courant,Limiter)

#**************************************************************************
def HartenYeeUpwind(i, Uo, E, Eps, Courant, Limiter):
    '''Returns the value of flux limiter function phi at (i+1/2) and (i-1/2)
    using the Harten-Yee Upwind TVD Limiters.

    Call signature:

        HartenYeeUpwind(i, Uo, E, Eps, Courant, Limiter)

    Parameters
    ----------

    i : int

        Integer value from the for loop to advance in space.

    Uo : 1D or 2D array

         The dependent variable from time level (n) within the domain.

    E : 1D or 2D array

        The flux vector for the non-linear term in the inviscid Burgers
        equation, which is E = U^2/2

    Eps : float

          A positive constant value within the range 0.0 and 0.125.

    Courant : float

              Courant number (entered as user input in file).

    Limiter : str

               The limiter for the TVD function.
               Options available for Harten-Yee Upwind TVD Limiters = "G".

    Returns
    -------

    phiPlus : float

              Flux limiter function at i+1/2 location.

    phiMinus : float

               Flux limiter function at i-1/2 location.
    '''
    import nanpack.secondaryfunctions as sf
    import nanpack.limiters as tvd

    dUiPlus12 = sf.CalcUi(Uo[i+1], Uo[i])
    dUiPlus32 = sf.CalcUi(Uo[i+2], Uo[i+1])

    dUiMinus12 = sf.CalcUi(Uo[i], Uo[i-1])
    dUiMinus32 = sf.CalcUi(Uo[i-1], Uo[i-2])

    alphaiPlus12 = sf.CalcAlpha(E[i+1], E[i], dUiPlus12,\
                                Uo[i+1], Uo[i])
    alphaiMinus12 = sf.CalcAlpha(E[i], E[i-1], dUiMinus12,\
                                 Uo[i], Uo[i-1])
    # .............................................................
    # Calculate some variables required in Harten-Yee Upwind.
    # .............................................................
    alphaiPlus32 = sf.CalcAlpha(E[i+2], E[i+1], dUiPlus32,\
                                Uo[i+2], Uo[i+1])
    alphaiMinus32 = sf.CalcAlpha(E[i-1], E[i-2], dUiMinus32,\
                                 Uo[i-1], Uo[i-2])

    if not Limiter == 'G':
        raise Exception('Invalid TVD limiter for Harten-Yee\
 Upwind TVD.\nUse Limiter = "G".')

    # .............................................................
    # Equation 6-130
    Gi = tvd.LimiterGforHYU(alphaiPlus12, alphaiMinus12, dUiPlus12,\
                            dUiMinus12, Courant, Eps)
    GiPlus1 = tvd.LimiterGforHYU(alphaiPlus32, alphaiPlus12, dUiPlus32,\
                                 dUiPlus12, Courant, Eps)
    GiMinus1 = tvd.LimiterGforHYU(alphaiMinus12, alphaiMinus32, dUiMinus12,\
                                  dUiMinus32, Courant, Eps)
    # Equation 6-129
    if dUiPlus12 != 0:
        betaiPlus12 = (GiPlus1 - Gi)/dUiPlus12
    else:
        betaiPlus12 = 0.0

    if dUiMinus12 != 0:
        betaiMinus12 = (Gi - GiMinus1)/dUiMinus12
    else:
        betaiMinus12 = 0.0

    # Calculate function si(alpha + beta) in Equation 6-126
    abPlus = alphaiPlus12 + betaiPlus12
    abMinus = alphaiMinus12 + betaiMinus12
    siPlus = sf.EntropyCorrectionFunction(abPlus, Eps)
    siMinus = sf.EntropyCorrectionFunction(abMinus, Eps)

    # Calculate the flux limiter function, Equation 6-126
    phiPlus = (GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = (Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus

#**************************************************************************
def ModifiedHartenYeeUpwind(i, Uo, E, Eps, Courant, Limiter):
    '''Returns the value of flux limiter function phi at (i+1/2) and (i-1/2)
    using the Modified Harten-Yee Upwind TVD Limiters.

    Call signature:

        ModifiedHartenYeeUpwind(i, Uo, E, Eps, Courant, Limiter)

    Parameters
    ----------

    i : int

        Integer value from the for loop to advance in space.

    Uo : 1D or 2D array

         The dependent variable from time level (n) within the domain.

    E : 1D or 2D array

        The flux vector for the non-linear term in the inviscid Burgers
        equation, which is E = U^2/2

    Eps : float

          A positive constant value within the range 0.0 and 0.125.

    Courant : float

              Courant number (entered as user input in file).

    Limiter : str

               The limiter for the TVD function.
               Options available for Modified Harten-Yee Upwind TVD Limiters
               are "G1", "G2", "G3", "G4", "G5".

    Returns
    -------

    phiPlus : float

              Flux limiter function at i+1/2 location.

    phiMinus : float

               Flux limiter function at i-1/2 location.
    '''
    import nanpack.secondaryfunctions as sf
    import nanpack.limiters as tvd

    dUiPlus12 = sf.CalcUi(Uo[i+1], Uo[i])
    dUiPlus32 = sf.CalcUi(Uo[i+2], Uo[i+1])

    dUiMinus12 = sf.CalcUi(Uo[i], Uo[i-1])
    dUiMinus32 = sf.CalcUi(Uo[i-1], Uo[i-2])

    # Equation 6-128
    alphaiPlus12 = sf.CalcAlpha(E[i+1], E[i], dUiPlus12,\
                                Uo[i+1], Uo[i])
    alphaiMinus12 = sf.CalcAlpha(E[i], E[i-1], dUiMinus12,\
                                 Uo[i], Uo[i-1])

    # .................................................................
    # Calculate some variables required in Modified Harten-Yee Upwind.
    # .................................................................
    alphaiPlus32 = sf.CalcAlpha(E[i+2], E[i+1], dUiPlus32,\
                                Uo[i+2], Uo[i+1])
    alphaiMinus32 = sf.CalcAlpha(E[i-1], E[i-2], dUiMinus32,\
                                 Uo[i-1], Uo[i-2])

    # .............................................................
    if not Limiter in ["G1", "G2", "G3", "G4", "G5"]:
        raise Exception('Invalid TVD limiter for the Modified Harten\
-Yee Upwind TVD.\nValid options are "G1"/"G2"/"G3"/"G4"/"G5".')

    # .............................................................
    Gi = tvd.LimiterforHYU(dUiPlus12, dUiMinus12, Limiter)
    GiPlus1 = tvd.LimiterforHYU(dUiPlus32, dUiPlus12, Limiter)
    GiMinus1 = tvd.LimiterforHYU(dUiMinus12, dUiMinus32, Limiter)

    # Calculate si(alpha) and sigma(si(alpha))
    siAlphaP = sf.EntropyCorrectionFunction(alphaiPlus12, Eps)
    siAlphaM = sf.EntropyCorrectionFunction(alphaiMinus12, Eps)
    sigmaP = 0.5*siAlphaP + Courant*alphaiPlus12*alphaiPlus12
    sigmaM = 0.5*siAlphaM + Courant*alphaiMinus12*alphaiMinus12

    # Calculate betaiPlus12
    if dUiPlus12 != 0:
        betaiPlus12 = sigmaP*(GiPlus1 - Gi)/dUiPlus12
    else:
        betaiPlus12 = 0.0

    # Calculate betaiMinus12
    if dUiMinus12 != 0:
        betaiMinus12 = sigmaM*(Gi - GiMinus1)/dUiMinus12
    else:
        betaiMinus12 = 0.0

    # Calculate function si(alpha + beta) in Equation 6-131
    abPlus = alphaiPlus12 + betaiPlus12
    abMinus = alphaiMinus12 + betaiMinus12
    siPlus = sf.EntropyCorrectionFunction(abPlus, Eps)
    siMinus = sf.EntropyCorrectionFunction(abMinus, Eps)

    # Calculate the flux limiter function, Equation 6-131
    phiPlus = sigmaP*(GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = sigmaM*(Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus

#**************************************************************************
def RoeSwebyUpwind(i, Uo, E, Eps, Courant, Limiter):
    '''Returns the value of flux limiter function phi at (i+1/2) and (i-1/2)
    using the Roe-Sweby Upwind TVD Limiters.

    Call signature:

        RoeSwebyUpwind(i, Uo, E, Eps, Courant, Limiter)

    Parameters
    ----------

    i : int

        Integer value from the for loop to advance in space.

    Uo : 1D or 2D array

         The dependent variable from time level (n) within the domain.

    E : 1D or 2D array

        The flux vector for the non-linear term in the inviscid Burgers
        equation, which is E = U^2/2

    Eps : float

          A positive constant value within the range 0.0 and 0.125.

    Courant : float

              Courant number (entered as user input in file).

    Limiter : str

               The limiter for the TVD function.
               Options available for Roe-Sweby Upwind TVD limiters are
               "G1", "G2", "G3".

    Returns
    -------

    phiPlus : float

              Flux limiter function at i+1/2 location.

    phiMinus : float

               Flux limiter function at i-1/2 location.
    '''
    import nanpack.secondaryfunctions as sf
    import nanpack.limiters as tvd

    dUiPlus12 = sf.CalcUi(Uo[i+1], Uo[i])

    dUiMinus12 = sf.CalcUi(Uo[i], Uo[i-1])

    # Equation 6-128
    alphaiPlus12 = sf.CalcAlpha(E[i+1], E[i], dUiPlus12,\
                                Uo[i+1], Uo[i])
    alphaiMinus12 = sf.CalcAlpha(E[i], E[i-1], dUiMinus12,\
                                 Uo[i], Uo[i-1])

    # .............................................................
    # Calculate some variables required in Row-Sweby Upwind.
    # .............................................................
    zero_filter = 1.e-7 # variable to filter out division by zero
    if alphaiPlus12 != 0:
        sigPlus12 = alphaiPlus12/abs(alphaiPlus12)
    else:
        sigPlus12 = 0.0

    if alphaiMinus12 != 0:
        sigMinus12 = alphaiMinus12/abs(alphaiMinus12)
    else:
        sigMinus12 = 0.0

    if dUiPlus12 >= zero_filter:
        riPlus = (Uold[i+1+int(sigPlus12)] -\
                  Uold[i+int(sigPlus12)])/dUiPlus12
    else:
        riPlus = 0.0

    if dUiMinus12 >= zero_filter:
        riMinus = (Uold[i+int(sigMinus12)] -\
                   Uold[i+1+int(sigMinus12)])/dUiMinus12
    else:
        riMinus = 0.0

    # .............................................................
    if not Limiter in ["G1", "G2", "G3"]:
        raise Exception('Invalid TVD limiter for the Roe-\
Sweby Upwind TVD.\nValid options are "G1"/"G2"/"G3".')

    # .............................................................
    Gi = tvd.LimiterforRSU(riPlus, Limiter)
    GiMinus1 = tvd.LimiterforRSU(riMinus, Limiter)

    # Calculate the flux limiter function, Equation 6-137
    phiPlus = ((Gi/2.0)*\
               (abs(alphaiPlus12) + Courant*alphaiPlus12**2)\
               - abs(alphaiPlus12))*dUiPlus12
    phiMinus = ((GiMinus1/2.0)*\
                (abs(alphaiMinus12) + Courant*alphaiMinus12**2)\
                - abs(alphaiMinus12))*dUiMinus12

    return phiPlus, phiMinus

#**************************************************************************
def DavisYeeSymmetric(i, Uo, E, Eps, Courant, Limiter):
    '''Returns the value of flux limiter function phi at (i+1/2) and (i-1/2)
    using the Davis-Yee Symmetric Limiter function.

    Call signature:

        DavisYeeSymmetric(i, Uo, E, Eps, Courant, Limiter)

    Parameters
    ----------

    i : int

        Integer value from the for loop to advance in space.

    Uo : 1D or 2D array

         The dependent variable from time level (n) within the domain.

    E : 1D or 2D array

        The flux vector for the non-linear term in the inviscid Burgers
        equation, which is E = U^2/2

    Eps : float

          A positive constant value within the range 0.0 and 0.125.

    Courant : float

              Courant number (entered as user input in file).

    Limiter : str

               The limiter for the TVD function.
               Options available for Davis-Yee Symmetric TVD limiters
               are "G1", "G2", "G3".

    Returns
    -------

    phiPlus : float

              Flux limiter function at i+1/2 location.

    phiMinus : float

               Flux limiter function at i-1/2 location.
    '''
    import nanpack.secondaryfunctions as sf
    import nanpack.limiters as tvd

    dUiPlus12 = sf.CalcUi(Uo[i+1], Uo[i])
    dUiPlus32 = sf.CalcUi(Uo[i+2], Uo[i+1])

    dUiMinus12 = sf.CalcUi(Uo[i], Uo[i-1])
    dUiMinus32 = sf.CalcUi(Uo[i-1], Uo[i-2])

    alphaiPlus12 = sf.CalcAlpha(E[i+1], E[i], dUiPlus12,\
                                Uo[i+1], Uo[i])
    alphaiMinus12 = sf.CalcAlpha(E[i], E[i-1], dUiMinus12,\
                                 Uo[i], Uo[i-1])

    # .............................................................
    if not Limiter in ["G1", "G2", "G3"]:
        raise Exception('Invalid TVD limiter for the Davis-\
Yee Symmetric TVD.\nValid options are "G1"/"G2"/"G3".')

    # .............................................................
    GiPlus12 = tvd.LimiterforDYS(dUiMinus12, dUiPlus12, dUiPlus32, Limiter)
    GiMinus12 = tvd.LimiterforDYS(dUiMinus32, dUiMinus12, dUiPlus12, Limiter)

    # Calculate function si(alpha) in Equation 6-141
    siPlus = sf.EntropyCorrectionFunction(alphaiPlus12, Eps)
    siMinus = sf.EntropyCorrectionFunction(alphaiMinus12, Eps)
    # Calculate the flux limiter function, Equation 6-141
    phiPlus = -((Courant*alphaiPlus12**2*GiPlus12) +\
                (siPlus*(dUiPlus12 - GiPlus12)))
    phiMinus = -((Courant*alphaiMinus12**2*GiMinus12) +\
                (siMinus*(dUiMinus12 - GiMinus12)))

    return phiPlus, phiMinus

#**************************************************************************







