# coding: utf-8
def siFunc(Param, Ep):
    '''Calculate Equation 6-127 in Hoffmann Vol. 1
    '''
    y = Param
    yAbs = abs(y)
    if yAbs >= Eps:
        si = y
    else:
        si = (y**2 + Eps**2)/(2*Eps)
        
    return si
