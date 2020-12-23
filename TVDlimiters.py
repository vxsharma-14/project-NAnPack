# coding: utf-8
def LimiterGforHYUpwind(convX,alpha1,alpha2,dU1,dU2,Ep):
    '''Calculate the limiter in Harten-Yee TVD using
       Equation 6-130 in Hoffmann Vol. 1.
    '''
    # Calculate si(alpha) in sigma and S
    siAlpha1 = siFunc(alpha1,Ep)
    siAlpha2 = siFunc(alpha2,Ep)
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
