# coding: utf-8
def HartenYeeTVD(U, convX):
    '''Calculate second-order TVD Limiter based on Harten-
       Yee Upwind TVD Limiters.
    '''
    Error = 0.0
    Uold = U
    
    for i in range (2,gm.iMax-2):
    
        dUiPlus12 = Uold[i+1] - Uold[i]
        dUiMinus12 = Uold[i] - Uold[i-1]
        dUiPlus12Abs = abs(dUiPlus12)
        dUiMinus12Abs = abs(dUiMinus12)
        
        dUiPlus32 = Uold[i+2] - Uold[i+1]
        dUiMinus32 = Uold[i-1] - Uold[i-2]
        dUiPlus32Abs = abs(dUiPlus32)
        dUiMinus32Abs = abs(dUiMinus32)
        
        EiPlus1 = 0.5*Uold[i+1]*Uold[i+1]
        Ei = 0.5*Uold[i]*Uold[i]
        EiMinus1 = 0.5*Uold[i-1]*Uold[i-1]
        EiPlus2 = 0.5*Uold[i+2]*Uold[i+2]
        EiMinus2 = 0.5*Uold[i-2]*Uold[i-2]
        

        # Equation 6-128
        if dUiPlus12 != 0:
            alphaiPlus12 = (EiPlus1 - Ei)/dUiPlus12
        else:
            alphaiPlus12 = 0.5*(Uold[i+1] + Uold[i])

        if dUiMinus12 != 0:
            alphaiMinus12 = (Ei - EiMinus1)/dUiMinus12
        else:
            alphaiMinus12 = 0.5*(Uold[i] + Uold[i-1])
            
        if dUiPlus32 != 0:
            alphaiPlus32 = (EiPlus2 - EiPlus1)/dUiPlus32
        else:
            alphaiPlus32 = 0.5*(Uold[i+2] + Uold[i+1])

        if dUiMinus32 != 0:
            alphaiMinus32 = (EiMinus1 - EiMinus2)/dUiMinus32
        else:
            alphaiMinus32 = 0.5*(Uold[i-1] + Uold[i-2])
            
        Gi = LimiterGforHYUpwind(convX,alphaiPlus12,alphaiMinus12,\
                                 dUiPlus12,dUiMinus12,Eps)
        GiPlus1 = LimiterGforHYUpwind(convX,alphaiPlus32,alphaiPlus12,\
                                      dUiPlus32,dUiPlus12,Eps)
        GiMinus1 = LimiterGforHYUpwind(convX,alphaiMinus12,alphaiMinus32,\
                                       dUiMinus12,dUiMinus32,Eps)
        
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

        # Equation 6-124 and 6-125 in Hoffmann Vol. 1
        hPlus = 0.5*(EiPlus1 + Ei + phiPlus)
        hMinus = 0.5*(Ei + EiMinus1 + phiMinus)

        # Equation 6-123    
        U[i] = Uold[i] - convX*(hPlus - hMinus)
        Error = Error + abs(Uold[i] - U[i]) # Calculate error

    return U, Error
