# coding: utf-8
#*************************************************************
def LimiterGforHYUpwind(convX,alpha1,alpha2,dU1,dU2,Ep):
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
    G = S*max(0, min(term1, term2))  
    
    return G

#*************************************************************
def LimiterG1forHYUpwind(dU1,dU2):
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
    G = S*max(0, min(term1, term2))  
    
    return G

#*************************************************************
def LimiterG2forHYUpwind(dU1,dU2):
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
def LimiterG3forHYUpwind(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-134 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    omeg = 1.e-7 # use between 1.e-7 and 1.e-5
    term1 = dU2*(dU1*dU1 + omeg)
    term2 = dU1*(dU2*dU2 + omeg)
    denom = dU1*dU1 + dU2*dU2 + 2*omeg

    # Equation 6-134
    G = (term1 + term2)/denom  
    
    return G

#*************************************************************
def LimiterG4forHYUpwind(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-135 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    if dU2 != 0:
        S = dU2/abs(dU2)
    else:
        S = 0.0

    term1 = abs(2*dU2)
    term2 = S*2*dU1
    term3 = S*0.5*(dU1 + dU2)

    # Equation 6-135
    G = S*max(0, min(term1, term2, term3))

    return G

#*************************************************************
def LimiterG5forHYUpwind(dU1,dU2):
    '''Calculate the Harten-Yee Upwind TVD limiter in Equation
       6-136 in Hoffmann Vol. 1 for the modified flux limiter
       function given by Equation 6-131.
    '''
    if dU1 != 0:
        S = dU1/abs(dU1)
    else:
        S = 0.0
    
    term1 = 2*abs(dU1)
    term2 = S*dU2
    term3 = abs(dU1)
    term4 = 2*S*dU2

    # Equation 6-136
    G = S*max(0, min(term1, term2), min(term3, term4))

    return G
#*************************************************************

