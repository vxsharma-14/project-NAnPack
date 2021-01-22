# coding: utf-8
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#               FLUX LIMITER FUNCTIONS FOR TVD               +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. HARTEN YEE UPWIND
# 2. MODIFIED HARTEN YEE UPWIND
# 3. ROE-SWEBY UPWIND
# 4. DAVIE-YEE SYMMETRIC

#*************************************************************
def HartenYeeUp(dUiPlus12, dUiMinus12, GiPlus1, Gi, GiMinus1,\
                alphaiPlus12, alphaiMinus12, Eps):
    '''Calculate the Harten-Yee Upwind TVD limiter function
       in Equation 6-126 in Hoffmann Vol. 1.
    '''
    from auxilliaryfunc import siFunc
    
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
    siPlus = siFunc(abPlus,Eps)
    siMinus = siFunc(abMinus,Eps)

    # Calculate the flux limiter function, Equation 6-126
    phiPlus = (GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = (Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus

#*************************************************************
def ModHartenYeeUp(dUiPlus12, dUiMinus12, GiPlus1, Gi, GiMinus1,\
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

    # Calculate the flux limiter function, Equation 6-131
    phiPlus = sigmaP*(GiPlus1 + Gi) - siPlus*dUiPlus12
    phiMinus = sigmaM*(Gi + GiMinus1) - siMinus*dUiMinus12

    return phiPlus, phiMinus

#*************************************************************
def RoeSwebyUp(dUiPlus12, dUiMinus12, Gi, GiMinus1,\
               alphaiPlus12, alphaiMinus12, convX):
    '''Calculate the Roe-Sweby Upwind TVD limiter function
       in Equation 6-137 in Hoffmann Vol. 1.
    '''
    # Calculate the flux limiter function, Equation 6-137
    phiPlus = ((Gi/2.0)*\
               (abs(alphaiPlus12) + convX*alphaiPlus12**2)\
               - abs(alphaiPlus12))*dUiPlus12
    phiMinus = ((GiMinus1/2.0)*\
                (abs(alphaiMinus12) + convX*alphaiMinus12**2)\
                - abs(alphaiMinus12))*dUiMinus12

    return phiPlus, phiMinus

#*************************************************************
def DavisYeeSym(dUiPlus12, dUiMinus12, GiPlus12, GiMinus12,\
                alphaiPlus12, alphaiMinus12, convX, Eps):
    '''Calculate the Davis-Yee Symmetric TVD limiter function
       in Equation 6-141 in Hoffmann Vol. 1.
    '''
    from auxilliaryfunc import siFunc
    
    # Calculate function si(alpha) in Equation 6-141
    siPlus = siFunc(alphaiPlus12,Eps)
    siMinus = siFunc(alphaiMinus12,Eps)
    # Calculate the flux limiter function, Equation 6-141
    phiPlus = -((convX*alphaiPlus12**2*GiPlus12) +\
                (siPlus*(dUiPlus12 - GiPlus12)))
    phiMinus = -((convX*alphaiMinus12**2*GiMinus12) +\
                (siMinus*(dUiMinus12 - GiMinus12)))

    return phiPlus, phiMinus

#*************************************************************

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                   LIMITER FUNCTIONS FOR TVD                +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. LIMITER (G) IN HARTEN YEE UPWIND
# 2. LIMITERS (G1...G5) FOR MODIFIED HARTEN YEE UPWIND
# 3. LIMITERS (G1...G3) FOR ROE-SWEBY UPWIND
# 4. LIMITERS FOR DAVIS-YEE SYMMETRIC

#*************************************************************
def LimiterGforHYU(convX,alpha1,alpha2,dU1,dU2,Ep):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-130 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-126.
    '''
    from auxilliaryfunc import siFunc
    
    # Calculate si(alpha) in sigma and S
    siAlpha1 = siFunc(alpha1,Ep)
    siAlpha2 = siFunc(alpha2,Ep)
    # Calculate sigma
    sigma1 = 0.5*(siAlpha1 - convX*alpha1**2)
    sigma2 = 0.5*(siAlpha2 - convX*alpha2**2)
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    # Equation 6-130
    term1 = sigma1*abs(dU1)
    term2 = S*sigma2*dU2
    G = S*max(0.0, min(term1, term2))  
    
    return G

#*************************************************************
def LimiterG1forHYU(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-132 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    if dU2 != 0:
        S = dU2/abs(dU2)
    else:
        S = 0.0
    # Equation 6-132
    term1 = abs(dU2)
    term2 = S*dU1
    G = S*max(0.0, min(term1, term2))  
    
    return G

#*************************************************************
def LimiterG2forHYU(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-133 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    term1 = dU1*dU2
    term2 = abs(term1)
    denom = dU1 + dU2

    if denom != 0:
        # Equation 6-133
        G = (term1 + term2)/denom
    else:
        G = 0.0  
    
    return G

#*************************************************************
def LimiterG3forHYU(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-134 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    omeg = 1.e-7 # use between 1.e-7 and 1.e-5
    term1 = dU2*(dU1*dU1 + omeg)
    term2 = dU1*(dU2*dU2 + omeg)
    denom = dU1*dU1 + dU2*dU2 + 2.0*omeg

    # Equation 6-134
    G = (term1 + term2)/denom  
    
    return G

#*************************************************************
def LimiterG4forHYU(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-135 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    if dU2 != 0:
        S = dU2/abs(dU2)
    else:
        S = 0.0

    term1 = abs(2.0*dU2)
    term2 = S*2.0*dU1
    term3 = S*0.5*(dU1 + dU2)

    # Equation 6-135
    G = S*max(0.0, min(term1, term2, term3))

    return G

#*************************************************************
def LimiterG5forHYU(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-136 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    
    term1 = 2.0*abs(dU1)
    term2 = S*dU2
    term3 = abs(dU1)
    term4 = 2.0*S*dU2

    # Equation 6-136
    G = S*max(0.0, min(term1, term2), min(term3, term4))

    return G

#*************************************************************
def LimiterG1forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-138 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-138
    G = max(0.0, min(1.0, r))
    
    return G

#*************************************************************
def LimiterG2forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-139 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-139
    G = (r + abs(r))/(1.0 + r)
    
    return G

#*************************************************************
def LimiterG3forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-140 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-140
    G = max(0.0, min(2.0*r, 1.0), min(r, 2.0))
    
    return G

#*************************************************************
def LimiterG1forDYS(dU1, dU2, dU3):
    '''Calculate the Davis-Yee Symmetric TVD limiter in Equation
       6-142 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-141.
    '''
    if 2.0*dU1 != 0:
        S = 2.0*dU1/abs(2.0*dU1)
    else:
        S = 0.0
    term1 = abs(2.0*dU1)
    term2 = S*2.0*dU2
    term3 = S*2.0*dU3
    term4 = S*0.5*(dU1 + dU3)
    # Equation 6-142
    G = S*max(0.0, min(term1, term2, term3, term4))
    
    return G

#*************************************************************
def LimiterG2forDYS(dU1, dU2, dU3):
    '''Calculate the Davis-Yee Symmetric TVD limiter in Equation
       6-143 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-141.
    '''
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    term1 = abs(dU1)
    term2 = S*dU2
    term3 = S*dU3
    # Equation 6-143
    G = S*max(0.0, min(term1, term2, term3))
    
    return G

#*************************************************************
def LimiterG3forDYS(dU1, dU2, dU3):
    '''Calculate the Davis-Yee Symmetric TVD limiter in Equation
       6-144 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-141.
    '''
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    term11 = abs(dU1) # term 1 in 1st minmod
    term12 = S*dU2 # term 2 in 1st minmod
    term21 = abs(dU1) # term 1 in 2nd minmod
    term22 = S*dU3 # term 2 in 2nd minmod
    term3 = S*dU2
    # Equation 6-144
    G = S*max(0.0, min(term11, term12)) +\
        S*max(0.0, min(term21, term22)) - term3
    
    return G

#*************************************************************
