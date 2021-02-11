#   ***********************************************************************
#
#   FILE         tutorial-03-diffusion-1D-solvers-all.py
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
"""Script to solve 1D diffusion equation using all parabolic solvers.

Code Description:
    We will repeat the example in "tutorial-02-diffusion-solvers-FTCS"
    using all available solvers. To perform this efficiently, we must
    automate the script which will be the objective of this tutorial.

    The available numerical methods for the 1D diffusion equation are:
        Explicit Forward Time Central Spacing (FTCS) method
        Explicit DuFort-Frankel method
        Implicit Laasonenn method
        Implicit Cran-Nicolson method

    I will post separate tutorials explaining these methods individually.
    You may also find a brief description on one of the blogs on my
    LinkedIn profile.
"""
# Import modules
import nanpack.preprocess as pre
from nanpack.grid import RectangularGrid
import nanpack.parabolicsolvers as pb
import nanpack.postprocess as post
from nanpack.benchmark import ParallelPlateFlow


def diffusion1D():
    """Compute the numerical solution using parabolic solvers."""
    config_file = "path/to/project/input/config.ini"
    cfg = pre.RunConfig(config_file)

    X, _ = RectangularGrid(cfg.dX, cfg.iMax)
    diffX, _ = pre.DiffusionNumbers(cfg.Dimension, cfg.diff, cfg.dT,
                                    cfg.dX)

    # Create a list of functions that we will be calling in this script
    func = [pb.FTCS, pb.DuFortFrankel, pb.Laasonen, pb.CrankNicolson]
    files = ["FTCS", "DuFortFrankel", "Laasonen", "CrankNicolson"]

    # Start a loop for 4 solver functions

    for f in range(len(func)):
        # Define initial conditions
        cfg.U[0] = 40.0
        cfg.U[1:] = 0.0
        # Define boundary conditions
        U = BC(cfg.U)
        # Start iterations
        Error = 1.0
        n = 0

        while n <= cfg.nMax and Error > cfg.ConvCrit:
            Error = 0.0
            n += 1
            if f == 1:  # DuFort Frankel requires an extra step

                if n == 1:
                    Uold = U.copy()
                    U = pb.FTCS(Uold, diffX)
                Uold2 = Uold.copy()
                Uold = U.copy()
                U = func[f](Uold, Uold2, diffX)
            else:
                Uold = U.copy()
                U = func[f](Uold, diffX)
            Error = post.AbsoluteError(U, Uold)
            # Update BC
            U = BC(U)
            # Write output to file
            fname = f"path/to/project/output/\
{files[f]}1D.dat"
            convfname = f"path/to/project/output/\
hist{files[f]}1D.dat"
            post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax, fname,
                                     cfg.dX)
            # Write convergence history log to a file
            post.WriteConvHistToFile(cfg, n, Error, convfname)

        # Write output to file
        post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax, fname, cfg.dX)
        # Write convergence history log to a file
        post.WriteConvHistToFile(cfg, n, Error, convfname)
        print()
    # Obtain analytical solution
    Uana = ParallelPlateFlow(40.0, X, cfg.diff, cfg.totTime, 20)
    ana_file = "path/to/project/output/analytical1D.dat"
    post.WriteSolutionToFile(Uana, 10, cfg.nWrite, cfg.nMax,
                             ana_file, cfg.dX)
    fn = []
    files.append("analytical")
    for f in range(len(files)):
        fn.append(f"path/to/project/output/\
{files[f]}1D.dat")
    post.Plot1DResults(dataFiles=fn,
                       uAxis="X",
                       Markers="Default",
                       Legend=files,
                       xLabel="Velocity (m/s)",
                       yLabel="Plate distance (m)",
                       Title=f"Comparison of Numerical Methods\nat t=\
{cfg.totTime} s.")


def BC(U):
    """Return the dependent variable with the updated boundary values."""
    U[0] = 40.0
    U[-1] = 0.0

    return U


if __name__ == "__main__":
    diffusion1D()
