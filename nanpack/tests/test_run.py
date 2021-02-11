#   ***********************************************************************
#
#   FILE         test_run.py
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

# System imports

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt

# nanpack imports
import nanpack.ellipticsolvers as ep
import nanpack.postprocess as pp
import nanpack.benchmark as bm
from nanpack.backend.initialize import InitialCondition
from nanpack.grid import RectangularGrid

class Properties:
    """Use this class to store input parameters."""
    def __init__(self):
        """Class constructor method."""
        self.display()
    def display(self):
        self.nDisplay = 2
        

def test_install():
    """Test nanpack installation by running example solution of Poisson's
    equation using ADISOR method.
    """
    # Define grid and domain
    iM = 21
    jM = 41
    dX = 0.05
    dY = 0.05
    Beta = dX/dY
    X, Y = RectangularGrid(dX, iM, dY, jM)
    cfg = Properties()
    print(cfg.nDisplay)

    # Initialize U
    T = InitialCondition('2D', iM, jM)

    # Boundary conditions
    T = BC(T)

    Error = 1.0
    ConvCrit = 0.01  # Convergence criteria
    nM = 100
    # Start iterations
    n = 0
    while n < nM and Error > ConvCrit:
        Error = 0.0
        n = n + 1
        Told = T.copy()
        T = ep.ADISOR(Told, Beta)
        Error = pp.AbsoluteError(T, Told)
        T = BC(T)
        pp.MonitorConvergence(cfg, n, Error)

    # Calculate Analytical solution
    Tana = bm.HeatConduction(X, Y, 100.0, 0.0, 0.0, 0.0, 20)
    print("Test run execution SUCCESS.")
    plot1D(Y, dX, T, Tana)


def BC(T):
    """Assign the boundary conditions for the test example solved in this
    module.
    """
    # Along the boundary edges
    T[:, 0] = 100.0
    T[:, -1] = 0.0
    T[0, :] = 0.0
    T[-1, 0] = 0.0

    # Along the boundary vertices
    T[0, 0] = 100.0
    T[-1, 0] = 100.0
    T[0, -1] = 0.0
    T[-1, -1] = 0.0

    return T


def plot1D(Y, dx, T, Tana):
    """Plot along 1D axis."""

    iM, jM = Y.shape
    Tn = np.zeros(jM, dtype="float")
    Ta = np.zeros(jM, dtype="float")
    y = np.zeros(jM, dtype="float")

    plt.rc("font", family="serif", size=8)
    fig, ax = plt.subplots(dpi=150)

    marker = ["-k", "o:r"]

    xloc = [0.2, 0.4, 0.6, 0.8]
    i = [int(x/dx) for x in xloc]

    plt_ttl = [f"X = {x}" for x in xloc]

    for p in range(1, 5):
        plt.subplot(2, 2, p)
        for j in range(0, jM):
            Tn[j] = T[i[p-1]][j]
            Ta[j] = Tana[i[p-1]][j]
            y[j] = Y[i[p-1]][j]
        if p == 1:
            plt.plot(Tn, y, marker[0], linewidth=0.5, label="Numerical",
                     markersize=2, markevery=5)
            plt.plot(Tn, y, marker[1], linewidth=0.5, label="Analytical",
                     markersize=2, markevery=5)
        if p == 2:
            plt.plot(Tn, y, marker[0], linewidth=0.5, label="Numerical",
                     markersize=2, markevery=5)
            plt.plot(Tn, y, marker[1], linewidth=0.5, label="Analytical",
                     markersize=2, markevery=5)

        if p == 3:
            plt.plot(Tn, y, marker[0], linewidth=0.5, label="Numerical",
                     markersize=2, markevery=5)
            plt.plot(Tn, y, marker[1], linewidth=0.5, label="Analytical",
                     markersize=2, markevery=5)

        if p == 4:
            plt.plot(Tn, y, marker[0], linewidth=0.5, label="Numerical",
                     markersize=2, markevery=5)
            plt.plot(Tn, y, marker[1], linewidth=0.5, label="Analytical",
                     markersize=2, markevery=5)

    # Format and customize plot
        plt.grid(which="major", axis="both", color="lightgrey",
                 linestyle=":", linewidth=0.5)
        plt.tight_layout()
        plt.xlabel('Temperature (K)')
        plt.ylabel('Domain height (m)')
        plt.title(f"{plt_ttl[p-1]} m", fontsize=8)
        plt.legend(fontsize=6)
    plt.subplots_adjust(left=0.02, right=0.95, hspace=0.45, wspace=0.3)
    plt.show()

if __name__ == '__main__':
    """Define test_install as the main function in test_run.py module."""
    test_install()
