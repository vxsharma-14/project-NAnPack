#   ***********************************************************************
#
#   FILE         tutorial-05-poissons-solvers.py
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

"""Script to solve Poissons equation using an elliptic solver.

Example Description:
    The steady-state temperature distribution is required to be obtained
    within a rectangular bar with the given boundary conditions.

    T at (t=0.0, x, y)     = 0.0   R  # Initial condition
    T at (t>0.0, x, y=0.0) = 100.0 R  # Boundary condition at lower wall
    T at (t>0.0, x, y=H)   = 0.0   R  # Boundary condition at upper wall
    T at (t>0.0, x=0.0, y) = 0.0   R  # Boundary condition at left wall
    T at (t>0.0, x=L, y)   = 0.0   R  # Boundary condition at right wall

    Length = 1 ft
    Heigth = 2 ft
    dX  = 0.05 ft
    dY  = 0.05 ft
    Convergence Criteria = 0.01
"""
# Import modules
import numpy as np
import matplotlib.pyplot as plt
from nanpack.grid import RectangularGrid
import nanpack.ellipticsolvers as ep
import nanpack.postprocess as post
import nanpack.preprocess as pre
from nanpack.benchmark import HeatConduction

configfile = "D:/MyProjects/projectroot/nanpack/input/config-steady\
-state-temp.ini"
# Create class instance
cfg = pre.RunConfig(configfile)
U = cfg.U

X, Y = RectangularGrid(cfg.dX, cfg.iMax, cfg.dY, cfg.jMax)
Beta = cfg.dX/cfg.dY

# start iterations
Error = 1.0
n = 0
while n < cfg.nMax and Error > cfg.ConvCrit:
    Error = 0.0
    n = n + 1
    Uold = U.copy()
    U = ep.PointGaussSeidel(Uold, Beta)
    Error = post.AbsoluteError(U, Uold)
    # Update BC
    U = pre.BC2D(U, cfg.BC, cfg.dX, cfg.dY)
    post.MonitorConvergence(cfg, n, Error)
    # Write output to file
    post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax,
                             cfg.OutFileName, cfg.dX, cfg.dY)
    # Write convergence history log to a file
    post.WriteConvHistToFile(cfg, n, Error)

# Write output to file
post.WriteSolutionToFile(U, n, cfg.nWrite, cfg.nMax,
                         cfg.OutFileName, cfg.dX, cfg.dY)
# Write convergence history log to a file
post.WriteConvHistToFile(cfg, n, Error)
post.WriteSolutionIn1DFormat(cfg, U)
# Calculate Analytical solution
Tana = HeatConduction(X, Y, 100.0, 0.0, 0.0, 0.0, 20)
file_name = "D:/MyProjects/projectroot/nanpack/output/analytic1Dx.dat"
post.WriteSolutionIn1DFormat(cfg, Tana, file_name)

# Plot from files
files = [
    "D:/MyProjects/projectroot/nanpack/output/pgs1Dx.dat",
    "D:/MyProjects/projectroot/nanpack/output/analytic1Dx.dat"
    ]
x = np.loadtxt(files[0], unpack=True, skiprows=2, usecols=(0))
tn = np.loadtxt(files[0], unpack=True, skiprows=2, usecols=(2))
ta = np.loadtxt(files[1], unpack=True, skiprows=2, usecols=(2))
plt.rc("font", family="serif", size=8)
fig, ax = plt.subplots(dpi=150)
plt.plot(tn, x, "-.b", linewidth=0.5, label="Numerical", markersize=2,
         markevery=5)
plt.plot(ta, x, ">r", linewidth=0.5, label="Analytical", markersize=2,
         markevery=5)
plt.xlabel('Temperature (deg R)')
plt.ylabel('Domain height (m)')
plt.title("Temperature distribution", fontsize=8)
plt.show()
