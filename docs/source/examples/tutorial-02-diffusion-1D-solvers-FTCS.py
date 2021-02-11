#   ***********************************************************************
#
#   FILE         tutorial-02-diffusion-1D-solvers-FTCS.py
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
"""Script to solve 1D diffusion equation using a parabolic solver.

Example Description:
    This is a classical application in fluid mechanics where the fluid is
    bounded between two parallel plates. The upper plate remains stationary
    and
    the lower plate is suddenly accelerated in x-direction at velocity u.
    It is required to find the velocity profile between the plates for the
    given initial and boundary conditions.

    U at (t=0.0, 0.0<y<H) = 0.0 m/s # Initial condition
    U at (t=0.0, y=0.0) = 40.0 m/s  # Initial condition
    U at (t>0.0, y=0.0) = 40.0 m/s  # Boundary condition at lower plate
    U at (t>0.0, y=H) = 0.0         # Boundary condition at upper plate

    Viscosity of fluid= 2.17e-4 m2/s
    H = 0.04 m
    dY = 0.001 m
"""
# Import modules
import matplotlib.pyplot as plt
import nanpack.preprocess as pre
from nanpack.grid import RectangularGrid
from nanpack.parabolicsolvers import FTCS
import nanpack.postprocess as post
from nanpack.benchmark import ParallelPlateFlow


def diffusion1D():
    """Compute the numerical solution."""
    config_file = "path/to/project/input/config.ini"
    cfg = pre.RunConfig(config_file)
    # Define initial conditions
    cfg.U[0] = 40.0
    cfg.U[1:] = 0.0

    U = BC(cfg.U)

    X, _ = RectangularGrid(cfg.dX, cfg.iMax)
    diffX, _ = pre.DiffusionNumbers(cfg.Dimension, cfg.diff, cfg.dT,
                                    cfg.dX)

    # Start iterations
    Error = 1.0
    n = 0

    while n <= cfg.nMax and Error > cfg.ConvCrit:
        Error = 0.0
        n += 1
        Uold = U.copy()
        U = FTCS(Uold, diffX)
        Error = post.AbsoluteError(U, Uold)
        # Update BC
        U = BC(U)
        post.MonitorConvergence(cfg, n, Error)
        # Write output to file
        post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax,
                                 cfg.OutFileName, cfg.dX)
        # Write convergence history log to a file
        post.WriteConvHistToFile(cfg, n, Error)

    # Write output to file
    post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax,
                             cfg.OutFileName, cfg.dX)
    # Obtain analytical solution
    Uana = ParallelPlateFlow(40.0, X, cfg.diff, cfg.totTime, 20)
    # Write convergence history log to a file
    post.WriteConvHistToFile(cfg, n, Error)
    plt.rc("font", family="serif", size=8)
    fig, ax = plt.subplots(dpi=150)
    plt.plot(U, X, ">-.b", linewidth=0.5, label="FTCS",
             markersize=5, markevery=5)
    plt.plot(Uana, X, "o:r", linewidth=0.5, label="Analytical",
             markersize=5, markevery=5)
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Plate distance (m)')
    plt.legend()
    plt.title(f"Velocity profile\nat t={cfg.totTime} sec", fontsize=8)
    plt.show()


def BC(U):
    """Return the dependent variable with the updated boundary values."""
    U[0] = 40.0
    U[-1] = 0.0

    return U


if __name__ == "__main__":
    diffusion1D()
