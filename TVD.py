#*******************************************************************************
def HartenYeeTVD(U, convX, Limiter, LimiterFunc):
    '''Calculate second-order TVD Limiter based on Harten-
       Yee Upwind TVD Limiters.
    '''
    from globalmod import iMax, Eps
    import TVDlimiters as tvd
    import TVDfluxlimiterfunc as flux
    from auxilliaryfunc import siFunc
    import sys
    Uold = U.copy()
            
    for i in range (2,iMax-2):

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

        if Limiter == 'G' and LimiterFunc == 'HY':
            
            # Equation 6-130
            Gi = tvd.LimiterGforHYUpwind(convX,alphaiPlus12,alphaiMinus12,\
                                         dUiPlus12,dUiMinus12,Eps)
            GiPlus1 = tvd.LimiterGforHYUpwind(convX,alphaiPlus32,alphaiPlus12,\
                                              dUiPlus32,dUiPlus12,Eps)
            GiMinus1 = tvd.LimiterGforHYUpwind(convX,alphaiMinus12,alphaiMinus32,\
                                               dUiMinus12,dUiMinus32,Eps)
            
            phiPlus, phiMinus = flux.LimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                       GiMinus1, alphaiPlus12, alphaiMinus12, Eps)
            
        elif Limiter == 'G1' and LimiterFunc == 'HYModified':
            
            # Equation 6-132
            Gi = tvd.LimiterG1forHYUpwind(dUiPlus12,dUiMinus12)
            GiPlus1 = tvd.LimiterG1forHYUpwind(dUiPlus32,dUiPlus12)
            GiMinus1 = tvd.LimiterG1forHYUpwind(dUiMinus12,dUiMinus32)
            
            phiPlus, phiMinus = flux.ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                          GiMinus1, alphaiPlus12, alphaiMinus12, convX, Eps)
            
        elif Limiter == 'G2' and LimiterFunc == 'HYModified':
            
            # Equation 6-133
            Gi = tvd.LimiterG2forHYUpwind(dUiPlus12,dUiMinus12)
            GiPlus1 = tvd.LimiterG2forHYUpwind(dUiPlus32,dUiPlus12)
            GiMinus1 = tvd.LimiterG2forHYUpwind(dUiMinus12,dUiMinus32)
            
            phiPlus, phiMinus = flux.ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                          GiMinus1, alphaiPlus12, alphaiMinus12, convX, Eps)
            
        elif Limiter == 'G3' and LimiterFunc == 'HYModified':
            
            # Equation 6-134
            Gi = tvd.LimiterG3forHYUpwind(dUiPlus12,dUiMinus12)
            GiPlus1 = tvd.LimiterG3forHYUpwind(dUiPlus32,dUiPlus12)
            GiMinus1 = tvd.LimiterG3forHYUpwind(dUiMinus12,dUiMinus32)
            
            phiPlus, phiMinus = flux.ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                          GiMinus1, alphaiPlus12, alphaiMinus12, convX, Eps)
            
        elif Limiter == 'G4' and LimiterFunc == 'HYModified':
            
            # Equation 6-135
            Gi = tvd.LimiterG4forHYUpwind(dUiPlus12,dUiMinus12)
            GiPlus1 = tvd.LimiterG4forHYUpwind(dUiPlus32,dUiPlus12)
            GiMinus1 = tvd.LimiterG4forHYUpwind(dUiMinus12,dUiMinus32)
            
            phiPlus, phiMinus = flux.ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                          GiMinus1, alphaiPlus12, alphaiMinus12, convX, Eps)
            
        elif Limiter == 'G5' and LimiterFunc == 'HYModified':
            
            # Equation 6-136
            Gi = tvd.LimiterG5forHYUpwind(dUiPlus12,dUiMinus12)
            GiPlus1 = tvd.LimiterG5forHYUpwind(dUiPlus32,dUiPlus12)
            GiMinus1 = tvd.LimiterG5forHYUpwind(dUiMinus12,dUiMinus32)
            
            phiPlus, phiMinus = flux.ModLimiterFuncforHYU(dUiPlus12, dUiMinus12, GiPlus1, Gi,\
                                                          GiMinus1, alphaiPlus12, alphaiMinus12, convX, Eps)
            
        else:
            print('****** Incorrect TVD Limiter and Limiter Function Selection **********.')
            sys.exit('Program terminating now.')                   

        # Equation 6-124 and 6-125 in Hoffmann Vol. 1
        hPlus = 0.5*(EiPlus1 + Ei + phiPlus)
        hMinus = 0.5*(Ei + EiMinus1 + phiMinus)

        # Equation 6-123    
        U[i] = Uold[i] - convX*(hPlus - hMinus)
        
    Error = abs(Uold[1:-1] - U[1:-1]).sum() # absolute error

    return U, Error

#*******************************************************************************
