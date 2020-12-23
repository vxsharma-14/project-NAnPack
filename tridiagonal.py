# coding: utf-8
def TridiagonalSolver(A, B, C, D, U):
    '''Solve a tridiagonal matrix for a system of linear
       equations in implicit methods for solution of one-
       dimensional equations.
    '''
    from globalmod import iMax
    H = U.copy()
    G = U.copy()
    H[0] = 0.0
    G[0] = U[0]
    for i in range (1,iMax-1):
        H[i] = C[i]/(B[i] - A[i]*H[i-1])
        G[i] = (D[i] - A[i]*G[i-1])/(B[i] - A[i]*H[i-1])
        
    for i in range (iMax-2,0,-1):
        U[i] = -H[i]*U[i+1] + G[i]
        
    return U
