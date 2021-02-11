"""Not a public module in version 1.0.0a4."""
#   ***********************************************************************
#
#   FILE         dimensionalize.py
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


class nondimensionalize:
    """A class to non-dimensionalize various numerical parameters.

    This is not complete or tested for accuracy.
    """

    def __init__(self):
        """Class construcotr."""

    def ndXgrid(self, X, L):
        """Non-dimensionalize x point locations along X axis."""
        self.xstar[:] = X[:]/L
        return self.xstar

    def ndYgrid(self, Y, L):
        """Non-dimensionalize y point locations along Y axis."""
        self.ystar[:] = Y[:]/L
        return self.ystar

    def ndvelU(self, U, L, nu):
        """Non-dimensionalize velocity U."""
        shapeU = U.shape
        if len(shapeU) == 1:
            im = shapeU
            self.ustar = [U[i]*L/nu for i in range(im)]
        elif len(shapeU) == 2:
            im, jm = shapeU
            self.ustar = [[U[i][j]*L/nu for j in range(jm)]
                          for i in range(im)]
        return self.ustar

    def ndtime(self, t, L, nu):
        """Non-dimensionalize time t."""
        L2 = L*L
        self.tstar = nu*t/L2
        return self.tstar


class dimensionalize:
    """A class to dimensionalize various numerical parameters.

    This is not complete or tested for accuracy.
    """

    def __init__(self):
        """Class constructor."""

    def dimXgrid(self, xstar, L):
        """Dimensionalize x point locations along X axis."""
        self.x[:] = xstar[:]*L
        return self.x

    def dimYgrid(self, ystar, L):
        """Dimensionalize y point locations along Y axis."""
        self.y[:] = ystar[:]*L
        return self.ystar

    def dimvelU(self, ustar, L, nu):
        """Dimensionalize velocity U."""
        shapeU = ustar.shape
        if len(shapeU) == 1:
            im = shapeU
            self.u = [ustar[i]*nu/L for i in range(im)]
        elif len(shapeU) == 2:
            im, jm = shapeU
            self.u = [[ustar[i][j]*nu/L for j in range(jm)]
                      for i in range(im)]
        return self.u

    def dimtime(self, tstar, L, nu):
        """Non-dimensionalize time t."""
        L2 = L*L
        self.t = tstar*L2/nu
        return self.t
