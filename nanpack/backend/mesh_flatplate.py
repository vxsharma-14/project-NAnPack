#   ***********************************************************************
#
#   FILE         mesh_flatplate.py
#
#   AUTHOR       Dr. Vishal Sharma
#
#   VERSION      1.0.0-alpha5
#
#   WEBSITE      https://github.com/vxsharma-14/project-NAnPack
#
#   NAnPack Learner's Edition is distributed under the MIT License.
#
#   Copyright (c) 2022 Vishal Sharma
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

from .util_clustering import calculate_beta1, calculate_eta1


def flateplate(A, B, x, y, dX, dY, iM, jM, Xi, Eta, beta):
    beta1 = calculate_beta1(beta)
    for j in range(0, jM):
        # Grid points along wall X = 0.0
        eta1 = calculate_eta1(Eta[0, j], Eta[0, -1])
        x[0][j] = 0.0
        y[0][j] = B * (((beta+1.0) - (beta-1.0)*(beta1**eta1))
                       / ((beta1**eta1)+1.0))
        # Grid points along wall X = A
        eta1 = calculate_eta1(Eta[-1, j], Eta[-1, -1])
        x[-1][j] = A
        y[-1][j] = B * (((beta+1) - (beta-1)*(beta1**eta1))
                        / ((beta1**eta1)+1.0))

        for i in range(1, iM-1):
            x[i][j] = x[i-1][j] + dX
            y[i][j] = y[0][j]

    return x, y
