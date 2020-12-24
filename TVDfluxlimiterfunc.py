# coding: utf-8
#************************************************************************
def LimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi, GiMinus1,\
                     alphaiPlus12, alphaiMinus12, Eps):
    '''Calculate the Harten-Yee Upwind TVD limiter function
       in Equation 6-126 in Hoffmann Vol. 1.
    '''
    # Equation 6-129
    if dUiPlus12 != 0:
        beta1Plus12 = (GiPlus1 - Gi)/dUiPlus12
    else:
        betaiPlus12 = 0.0

    if dUiMinus12 != 0:
        betaiMinus12 = (Gi - GiMinus1)/dUiMinus12
    else:
        betaiMinus12 = 0.0

    # Calculate function si(alpha + beta) in Equation 6-126
    abPlus = alphaiPlus12 + betaiPlus12
    abMinus = alphaiMinus12 + betaiMinus12
    siPlus = siFunc(abPlus,Eps)
    siMinus = siFunc(abMinus,Eps)

    # Ccalculate the flux limiter function, Equation 6-126
    phiPlus = (GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = (Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus

#************************************************************************
def ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi, GiMinus1,\
                         alphaiPlus12, alphaiMinus12, convX, Eps):
    '''Calculate the Modified Harten-Yee Upwind TVD limiter function
       in Equation 6-131 in Hoffmann Vol. 1.
    '''
    from auxilliaryfunc import siFunc
    
    # Calculate si(alpha) and sigma(si(alpha))
    siAlphaP = siFunc(alphaiPlus12,Eps)
    siAlphaM = siFunc(alphaiMinus12,Eps)
    sigmaP = 0.5*siAlphaP + convX*alphaiPlus12**2
    sigmaM = 0.5*siAlphaM + convX*alphaiMinus12**2
    
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
    siPlus = siFunc(abPlus,Eps)
    siMinus = siFunc(abMinus,Eps)

    # Ccalculate the flux limiter function, Equation 6-131
    phiPlus = sigmaP*(GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = sigmaM*(Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus
