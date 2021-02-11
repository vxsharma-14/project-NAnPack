"""Not a public module in version 1.0.0a4."""
#   ***********************************************************************
#
#   FILE         tridiagonal.py
#
#   AUTHOR       Dr. Vishal Sharma
#
#   VERSION      1.0.0-alpha4
#
#   WEBSITE      https://github.com/vxsharma-14/project-NAnPack
#
#   NAnPack Learner's Edition is distributed under the MIT License.
#
#   Copyright (c) 2020 Vishal Sharma
#
#   Permission is hereby granted, free of charge, to any person
#   obtaining a copy of this software and associated documentation
#   files (the "Software"), to deal in the Software without restriction,
#   including without limitation the rights to use, copy, modify, merge,
#   publish, distribute, sublicense, and/or sell copies of the Software,
#   and to permit persons to whom the Software is furnished to do so,
#   subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#
#   You should have received a copy of the MIT License along with
#   NAnPack Learner's Edition.
#
#   ***********************************************************************


def TridiagonalSolver(tMax, A, B, C, D, UU):
    """Solve a tridiagonal matrix for a system of linear equations.

    The function is used by implicit methods for the obtaining the
    numerical of solution of PDEs. The approach to solve
    such a system is given in Appendix B of CFD Vol. 1 by Klaus Hoffmann.

    The formulation for the implicit methods is of the form:

    A(i,n)U(i-1,n+1) + B(i,n)U(i,n+1) + C(i,n)U(i+1,n+1) = D(i,n)

    Call signature:

        TridiagonalSolver(U, diffX, diffY)

    Parameters
    ----------
    tMax : int

        Grid points in a given axis direction.

    A : 1D array

        Coefficient of u(i-1, n+1) in the implicit formulation.

    B : 1D array

        Coefficient of u(i, n+1) in the implicit formulation.

    C : 1D array

        Coefficient of u(i+1, n+1) in the implicit formulation.

    D : 1D array

        Right hand side equations in the implicit formulation.

    UU : 1D array

        The dependent variable at time level (n) within the domain
        along a given axis as required by the implicit method.

    Returns
    -------
    UU : 1D array

        The dependent variable within the domain calculated using
        equation (B-5) in CFD Vol.1 by Klaus Hoffmann.
    """
    H = [0 for t in range(tMax)]  # initialize H
    G = [0 for t in range(tMax)]  # initialize G

    H[0] = 0.0
    G[0] = UU[0]
    for t in range(1, tMax-1):
        # Equation B-8 and B-9 in CFD Vol. 1 by Klaus Hoffmann
        H[t] = C[t]/(B[t] - A[t]*H[t-1])
        G[t] = (D[t] - A[t]*G[t-1])/(B[t] - A[t]*H[t-1])

    for t in range(tMax-2, 0, -1):
        # Equation B-5 in CFD Vol. 1 by Klaus Hoffmann
        UU[t] = -H[t]*UU[t+1] + G[t]

    return UU
