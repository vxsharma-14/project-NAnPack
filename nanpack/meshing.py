"""A module consisting of various meshing functions."""
#   ***********************************************************************
#
#   FILE         mesh.py
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

from .backend.exceptions import InvalidValueError
import numpy as np


def CalcGridPoints(Dimension, Length, delX, Height=None, delY=None):
    """Return the grid points along X and Y direction in the mesh.

    Call signature:
        CalcGridPoints(Dimension, Length, delX, Height=None, delY=None)

    Parameters
    ----------
    Dimension: str
        Dimension of the domain. Allowed inputs are "1D" or "2D".
    Length: float
        Length of the domain.
    delX: float
        Grid step size along X-axis.
    Height: float
        Height of the domain. Value required for 2D applications.
    delY: float
        Grid step size along Y-axis. Value required for 2D applications.

    Returns
    -------
    iMax: int
        Number of grid points along X-axis within the domain.
    jMax: int
        Number of grid points along Y-axis within the domain. Returns 0 for
        1D applications.
    """
    iMax = int(Length/delX) + 1
    if Dimension.upper() == "2D":
        jMax = int(Height/delY) + 1
    else:
        jMax = 0
    print("Calculating grid size: Completed.")

    return iMax, jMax


def CalcGridStepsize(Dimension, Length, iMax, Height=None, jMax=None):
    """Return the uniform grid steps size along X and Y axis.

    Call signature:
        CalcGridStepsize(Dimension, Length, iMax, Height=None, jMax=None)

    Parameters
    ----------
    Dimension: str
        Dimension of the domain. Allowed inputs are "1D" or "2D".
    Length: float
        Length of the domain.
    iMax: int
        Number of grid points along X-axis within the domain.
    Height: float
        Height of the domain. Value required for 2D applications.
    jMax: int
        Number of grid points along Y-axis within the domain. Value
        required for 2D applications.

    Returns
    -------
    delX: float
        Grid step size along X-axis.
    delY: float
        Grid step size along Y-axis. Returns 0.0 for 1D applications.
    """
    delX = Length/(iMax - 1)
    if Dimension.upper() == "2D":
        delY = Height/(jMax - 1)
    else:
        delY = 0.0
    print("Calculating grid step size: Completed.")

    return delX, delY


def RectangularMesh(dX, iMax, dY=None, jMax=None):
    """Return a rectangular uniform rectangular mesh.

    X and/or Y grid point locations are computed in a cartesian coordinate
    system using the grid step size and grid points.

    Call Signature:
        RectangularMesh(dX, iMax, dY=None, jMax=None)

    Parameters
    ----------
    dX: float
        Grid step size along X-axis.
    iMax : int
        Number of grid points along X-axis within the domain.
    dY: float
        Grid step size along Y-axis. Value required for 2D applications.
    jMax : int
        Number of grid points along Y-axis within the domain. Value
        required for 2D applications.

    Returns
    -------
    X: ndarray[float], <=2d
        Returns X coordinates at each grid points locations.
    Y: ndarray[float], <=2d
        Returns Y coordinates at each grid points locations. Returns 0 for
        1D applications.
    """
    print("Uniform rectangular grid generation in cartesian\
 coordinate system: Completed.")
    if isinstance(dY, float) and isinstance(jMax, int):
        X = np.zeros((iMax, jMax), dtype="float")
        Y = np.zeros((iMax, jMax), dtype="float")
        for i in range(0, iMax):
            for j in range(0, jMax):
                X[i][j] = i*dX
                Y[i][j] = j*dY
        return X, Y
    else:
        X = np.zeros(iMax, dtype="float")
        for i in range(0, iMax):
            X[i] = i*dX
        return X


def StructuredMesh1D(ClustLoc, ClustOpt=True, Alpha=None, Beta=None,
                     CfgClsObj=None, **mesh_kwargs):
    """Return a rectangular uniform/non-uniform rectangular mesh.

    Documentation incomplete. This routine is under construction.
    """
    print("Calculating X locations of all grid points within\
 the mesh.")
    dXi = 1.0
    if CfgClsObj is None:
        dX = mesh_kwargs.get("dX", None)
        iM = mesh_kwargs.get("iM", None)
    else:
        dX = CfgClsObj.dX
        iM = CfgClsObj.iM
    Xi = np.zeros((iM), dtype="float")
    for i in range(0, iM):
        Xi[i] = i*dXi
    if ClustOpt is False:
        X = RectangularMesh(dX, iM)
    else:
        X = _getMesh1D(ClustLoc, Xi, dX, iM, Alpha, Beta)
    return X


def _getMesh1D(clust_loc, Xi, dX, iM, alpha, beta):
    """Get X and Y coordinate locations within the user specified grid."""
    from .backend.mesh1d import meshing_func
    X = meshing_func(clust_loc, Xi, dX, iM, alpha, beta)
    return X


def StructuredMesh(GeomTemplate, ClustOpt=True, CfgClsObj=None,
                   **mesh_kwargs):
    """Return a rectangular uniform/non-uniform rectangular mesh.

    Documentation incomplete. This routine is under construction.
    """
    print("Calculating X and Y locations of all grid points within\
 the mesh.")
    rect_grid_types = ["flat-plate", "duct", "cavity"]
    if CfgClsObj is None:
        dX = mesh_kwargs.get("dX", None)
        dY = mesh_kwargs.get("dY", None)
        iM = mesh_kwargs.get("iM", None)
        jM = mesh_kwargs.get("jM", None)
    else:
        dX = CfgClsObj.dX
        dY = CfgClsObj.dY
        iM = CfgClsObj.iM
        jM = CfgClsObj.jM
    if GeomTemplate.lower() in rect_grid_types and ClustOpt is False:
        X, Y = RectangularMesh(dX, iM, dY, jM)
    else:
        X, Y = getMesh(iM, jM, GeomTemplate, ClustOpt, CfgClsObj,
                       **mesh_kwargs)
    return X, Y


def getMesh(iM, jM, geo_temp, clust_option, CfgObject, **mesh_kw):
    """Get X and Y coordinate locations within the user specified grid."""
    from .backend.mesh2d import meshing_func
    X, Y = meshing_func(iM, jM, geo_temp, clust_option, CfgObject,
                        **mesh_kw)
    return X, Y


def CalcMeshMetrics(X, Y):
    """Calculate metrics and Jacobian of the transformation."""
    from .backend.meshmetrics import metrics_2d
    XiX, XiY, EtaX, EtaY, JJ = metrics_2d(X, Y)
    print("Grid metrics and Jacobian evaluation: Completed.")
    return XiX, XiY, EtaX, EtaY, JJ


def CalcMeshMetrics1D(X):
    """Calculate metrics and Jacobian of the transformation."""
    from .backend import meshmetrics
    XiX, JJ = meshmetrics.metrics_1d(X)
    print("Grid metrics and Jacobian evaluation: Completed.")
    return XiX, JJ


def PlotMeshMetrics(XiX, XiY, EtaX, EtaY, x=None, y=None):
    """Plot metrics data."""
    from .backend.plotmetrics import plot_metrics_2d
    plot_metrics_2d(XiX, XiY, EtaX, EtaY, x, y)


def SaveMeshMetrics(X, Y, XiX, XiY, EtaX, EtaY, JJ, MetFName):
    """Save metrics data to an output file."""
    from .backend.writefiles import save_metrics_2d
    save_metrics_2d(X, Y, XiX, XiY, EtaX, EtaY, JJ, MetFName)


def SaveMeshMetrics1D(X, XiX, JJ, MetFName):
    """Save metrics data to an output file."""
    from .backend.writefiles import save_metrics_1d
    save_metrics_1d(X, XiX, JJ, MetFName)


def PlotMesh(X, Y):
    """Plot grid data."""
    from .backend.plotmetrics import plot_grid
    plot_grid(X, Y)


def CalcTimeStep(CFL, diff, conv, dX, dY, Dimension, Model):
    """Return the time step size in the numerical approximation.

    Call Signature:
        CalcTimeStep(CFL, diff, conv, dX, dY, Dimension, Model)

    Parameters
    ----------
    CFL: float
        In this program, CFL is treated as the
            diffusion number for diffusion equations, and
            Courant number for the convection equations.
            Caution: This is not a true numerical definition of CFL though.
    diff: float
        Physics specific coefficient in the diffusion model.
        For example, kinematic viscosity or thermal diffusivity.
    conv: float
        Physics specific coefficient in the convection model.
        For example, speed of sound in the first-order linear wave eq.
    dX: float
        Grid step size along X-axis.
    dY: float
        Grid step size along Y-axis. Value required for 2D applications.
    Dimension: str
        Dimension of the domain. Allowed inputs are "1D" or "2D".
    Model: str
        Model of the governing equation. To see available options for this
        parameter, type the following command on your terminal
            python fetchoption.py "model"

    Returns
    -------
    TimeStep: float
       Time step in the model equation.
    """
    print("Calculating time step size for the simulation: Completed.")
    # ************** DIFFUSION EQN. ******************
    if Model.upper() == "DIFFUSION":
        dX2 = dX*dX
        if Dimension.upper() == "1D":
            TimeStep = CFL*dX2/diff
            return TimeStep
        elif Dimension.upper() == "2D":
            dY2 = dY*dY
            TimeStep = CFL*(1.0/((1/dX2) + (1/dY2)))/diff
            return TimeStep
    # ************** FIRST-ORDER WAVE EQN. *****************
    elif Model.upper() == "FO_WAVE":
        if Dimension.upper() == "1D":
            TimeStep = CFL*dX/conv
            return TimeStep
    # ************** BURGERS EQN. *****************
    elif Model.upper() == "INV_BURGERS":
        if Dimension.upper() == "1D":
            TimeStep = CFL*dX
            return TimeStep
    elif Model.upper() == "VISC_BURGERS":
        if Dimension.upper() == "1D":
            dX2 = dX*dX
            TimeStep = CFL*dX2
            return TimeStep


def CalcMaxSteps(State, nMax, dT, simTime):
    """Return the max iteration/time steps for the program to run.

    Call Signature:
        CalcMaxSteps(State, nMax, dT, simTime)

    Parameters
    ----------
    State: str
        State at which the final solution is desired. It can be
        steady-state or transient.
        To obtain solution at several intermediate time steps before
        convergence, use transient option and provide the time in
        configuration file at which the solution is desired. The
        program will calculate when to stop the solution.
        Available inputs are "STEADY" or "TRANSIENT"
    nMax: int
        Maximum number of iterations until which the program must seek
        convergence. If convergence is not achieved after going thtough
        nMax steps, the program will stop solving any further.
    dT: float
        Time step in the discretized equation. The value is auto calculated
         by the program from the CFL value during the configuration step.
    simTime: float
        Intermediate time before convergence at which numerical solution
        is required.

    Returns
    -------
    MaxSteps: int
        Maximum iteration/time steps for the program to run.
    """
    print("Calculating maximum iterations/steps for the simulation:\
 Completed.")
    if State.upper() == "TRANSIENT":
        if not simTime > 0.0:  # simulation time can't be negative
            raise InvalidValueError("SIM_TIME", simTime)
        try:
            MaxSteps = int(round(simTime/dT))
            return MaxSteps
        except dT:
            raise Exception("No time step provided.")
    elif State.upper() == "STEADY":
        MaxSteps = nMax
        return MaxSteps
