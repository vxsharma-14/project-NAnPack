'''
+**************************************************************************
+**************************************************************************
+
+   FILE         residual.py
+
+   AUTHOR       Vishal Sharma
+
+   VERSION      1.0.0.dev1
+
+   WEBSITE      https://vxsharma-14.github.io/NAnPack/
+
+   NAnPack Learner's Edition is distributed under the MIT License.
+
+   Copyright (c) 2020 Vishal Sharma
+
+   Permission is hereby granted, free of charge, to any person
+   obtaining a copy of this software and associated documentation
+   files (the "Software"), to deal in the Software without restriction,
+   including without limitation the rights to use, copy, modify, merge,
+   publish, distribute, sublicense, and/or sell copies of the Software,
+   and to permit persons to whom the Software is furnished to do so,
+   subject to the following conditions:
+
+   The above copyright notice and this permission notice shall be
+   included in all copies or substantial portions of the Software.
+
+   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
+   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
+   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
+   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
+   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
+   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
+   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
+   SOFTWARE.
+
+   You should have received a copy of the MIT License along with
+   NAnPack Learner's Edition.
+
+**************************************************************************
+**************************************************************************
'''
#**************************************************************************
def LInfNormError(U, Uold):
    '''Calculate L-infinity Norm error.
    '''
    shapeU = U.shape # Obtain shape for Dimension
    if len(shapeU) == 2: # Dimension = 2D
        Abs_err = abs(Uold[1:,1:] - U[1:,1:]).sum(axis=1)
        Error = max(Abs_err)

    elif len(shapeU) == 1: # Dimension = 1D
        Abs_err = abs(Uold[1:] - U[1:]).sum(axis=1)
        Error = max(Abs_err)
    
    return Error

#**************************************************************************
def AbsoluteError(U, Uold):
    '''Calculate absolute error.
    '''
    shapeU = U.shape # Obtain shape for Dimension
    if len(shapeU) == 2: # Dimension = 2D
        Error = abs(Uold[1:-1,1:-1] - U[1:-1,1:-1]).sum()

    elif len(shapeU) == 1: # Dimension = 1D
        Error = abs(Uold[1:-1] - U[1:-1]).sum()

    return Error

#**************************************************************************
def MonitorConvergence(n, nDisplay, Error):
    '''Monitor convergence of data
    '''
    
    if n == 0 or n == 1:
        print(f'{"ITER":>7} {"ERROR":>15}')
        print(f'{"----":>7} {"-----":>15}')
    
    if n % nDisplay == 0:
        print(f'{n:>7} {Error:>15.8f}')
    
