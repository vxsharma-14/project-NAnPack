#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                   DEFINE DICTIONARIES FOR LIMITERS                +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. DICTIONARY FOR LIMITERS IN HARTEN-YEE UPWIND AND MODIFIED TVD
# 2. DICTIONARY FOR LIMITERS IN ROE-SWEBY UPWIND TVD
# 3. DICTIONARY FOR LIMITERS IN DAVIS-YEE SYMMETRIC TVD

#**************************************************************************
def LimiterforHYU(dU1, dU2, Limiter):
    ''' A dictionary to call a function based on the required limiter for
    the Modified Harten-Yee Upwind TVD scheme.
    '''
    lim = {
        "G1" : LimiterG1forHYU,
        "G2" : LimiterG2forHYU,
        "G3" : LimiterG3forHYU,
        "G4" : LimiterG4forHYU,
        "G5" : LimiterG5forHYU
        }

    SelectedFunction = lim.get(Limiter)
    return SelectedFunction(dU1, dU2)

#**************************************************************************
def LimiterforRSU(r, Limiter):
    ''' A dictionary to call a function based on the required limiter for
    the Roe-Sweby Upwind TVD scheme.
    '''
    lim = {
        "G1" : LimiterG1forRSU,
        "G2" : LimiterG2forRSU,
        "G3" : LimiterG3forRSU
        }

    SelectedFunction = lim.get(Limiter)
    return SelectedFunction(r)

#**************************************************************************
def LimiterforDYS(dU1, dU2, dU3, Limiter):
    ''' A dictionary to call a function based on the required limiter for
    the Davis-Yee Symmetric TVD scheme.
    '''
    lim = {
        "G1" : LimiterG1forDYS,
        "G2" : LimiterG2forDYS,
        "G3" : LimiterG3forDYS
        }

    SelectedFunction = lim.get(Limiter)
    return SelectedFunction(dU1, dU2, dU3)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                   LIMITER FUNCTIONS FOR TVD                +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1. LIMITER (G) IN HARTEN YEE UPWIND
# 2. LIMITERS (G1...G5) FOR MODIFIED HARTEN YEE UPWIND
# 3. LIMITERS (G1...G3) FOR ROE-SWEBY UPWIND
# 4. LIMITERS FOR DAVIS-YEE SYMMETRIC

#**************************************************************************
def LimiterGforHYU(alpha1, alpha2, dU1, dU2, Courant, Ep):
    '''Returns the Harten-Yee Upwind TVD limiter in Equation
       6-130 in CFD Vol. 1 by Hoffmann for the flux limiter
       function given by Equation 6-126.
    '''
    import fluid.secondaryfunctions as sf
    
    # Calculate si(alpha) in sigma and S
    siAlpha1 = sf.EntropyCorrectionFunction(alpha1, Ep)
    siAlpha2 = sf.EntropyCorrectionFunction(alpha2, Ep)
    # Calculate sigma
    sigma1 = 0.5*(siAlpha1 - Courant*alpha1*alpha1)
    sigma2 = 0.5*(siAlpha2 - Courant*alpha2*alpha2)
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    # Equation 6-130
    term1 = sigma1*abs(dU1)
    term2 = S*sigma2*dU2
    G = S*max(0.0, min(term1, term2))  
    
    return G

#**************************************************************************
def LimiterG1forHYU(dU1, dU2):
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

#**************************************************************************
def LimiterG2forHYU(dU1, dU2):
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

#**************************************************************************
def LimiterG3forHYU(dU1, dU2):
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

#**************************************************************************
def LimiterG4forHYU(dU1, dU2):
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

#**************************************************************************
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

#**************************************************************************
def LimiterG1forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-138 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-138
    G = max(0.0, min(1.0, r))
    
    return G

#**************************************************************************
def LimiterG2forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-139 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-139
    G = (r + abs(r))/(1.0 + r)
    
    return G

#**************************************************************************
def LimiterG3forRSU(r):
    '''Calculate the Roe-Sweby Upwind TVD limiter in Equation
       6-140 in Hoffmann Vol. 1 for the flux limiter function
       given by Equation 6-137.
    '''
    # Equation 6-140
    G = max(0.0, min(2.0*r, 1.0), min(r, 2.0))
    
    return G

#**************************************************************************
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

#**************************************************************************
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

#**************************************************************************
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

#**************************************************************************

