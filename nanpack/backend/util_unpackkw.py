"""Unpack keyword arguments kwargs."""


def converg_kwargs(**kwargs):
    """Unpack keyword arguments related to monitor convergence."""
    ndisplay = kwargs.get("nDisplay", 10)
    return ndisplay


def writesolution_kwargs(**kwargs):
    """Unpack keyword arguments related to write-solution-to-file."""
    nwrite = kwargs.get("nWrite", 10)
    fname = kwargs.get("FileName", "solution.dat")
    nmax = kwargs.get("nMax", None)
    dx = kwargs.get("dX", None)
    dy = kwargs.get("dY", None)
    var_dict = {
        "nWrite": nwrite,
        "FileName": fname,
        "nMax": nmax,
        "dX": dx,
        "dY": dy
        }
    return var_dict


def write1dsolution_kwargs(**kwargs):
    """Unpack keyword arguments related to write-solution-in-1D-format."""
    fname = kwargs.get("FileName1D", "solution2d_to_1d.dat")
    dx = kwargs.get("dX", None)
    dy = kwargs.get("dY", None)
    ax = kwargs.get("Direction", None)
    ax_loc = kwargs.get("Locations", None)
    var_dict = {
        "FileName": fname,
        "dX": dx,
        "dY": dy,
        "ax_loc": ax_loc,
        "axis": ax
        }
    return var_dict


def writeconv_kwargs(**kwargs):
    """Unpack keyword arguments related to write-convergence-hist-to-file ."""
    ndisplay = kwargs.get("nDisplay", 10)
    fname = kwargs.get("HistFName", "history.dat")
    nmax = kwargs.get("nMax", None)
    state = kwargs.get("State", None)
    crit = kwargs.get("ConvCriteria", 0.01)
    stime = kwargs.get("SimTime", None)
    var_dict = {
        "nDisplay": ndisplay,
        "FileName": fname,
        "nMax": nmax,
        "State": state,
        "ConvCriteria": crit,
        "SimTime": stime
        }
    return var_dict


def configuration_kwargs(**kwargs):
    """Unpack keyword arguments related to simulation configuration."""
    Model = kwargs.get("Model")
    conv = kwargs.get("ConvCoeff")
    dX = kwargs.get("dX")
    dY = kwargs.get("dY")
    iM = kwargs.get("iM")
    jM = kwargs.get("jM")
    var_dict = {
        "Model": Model,
        "conv": conv,
        "dX": dX,
        "dY": dY,
        "iM": iM,
        "jM": jM
        }
    return var_dict


def hypmodel_kwargs(**kwargs):
    """Unpack keyword arguments related to hyperbolic solvers."""
    return
