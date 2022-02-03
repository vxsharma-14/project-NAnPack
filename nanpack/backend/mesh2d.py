"""."""
import numpy as np
from .util_meshkwargs import unpack_clust_kwargs
from .util_meshkwargs import unpack_geometric_kwargs
from .mesh_internalflow import duct, cavity
from .mesh_flatplate import flateplate
from .mesh_otype import ogrid_airfoil, ogrid_cylinder
from .mesh_bluntbody import blunt_body_cone, blunt_body_ellipt
from .exceptions import MeshingInputError, GeometryTemplateError


def meshing_func(iM, jM, geo_type, clust_opt, CfgObject, **mesh_kw):
    """Return a meshing function requested by the user."""
    geo_temp = {
        "blunt-body-cone": bluntcone_mesh,
        "blunt-body-ellipse": bluntellip_mesh,
        # "C-grid": _cgrid_grid,
        "cavity": cavity_mesh,
        "duct": duct_mesh,
        "flat-plate": flatplate_mesh,
        "o-grid-airfoil": o_airfoil_mesh,
        "o-grid-cylinder": o_cyl_mesh,
        # "Wind-tunnel-airfoil": _wt_airfoil_grid,
        # "Wind-tunnel-cylinder": _wt_cyl_grid,
        }
    SelectedGeometry = geo_temp.get(geo_type, undefined_geom)
    Xi, Eta = curvilinear_coordinates(iM, jM)

    return SelectedGeometry(clust_opt, Xi, Eta, iM, jM, CfgObject,
                            **mesh_kw)


def curvilinear_coordinates(iM, jM):
    """Return grid points on the curvilinear coord. Xi, Eta."""
    dXi = 1.0
    dEta = 1.0
    Xi = np.zeros((iM, jM), dtype="float")
    Eta = np.zeros((iM, jM), dtype="float")
    for i in range(0, iM):
        for j in range(0, jM):
            Xi[i][j] = i*dXi
            Eta[i][j] = j*dEta

    return Xi, Eta


def undefined_geom():
    """Return an exception when an invalid geomtery is requested.

    When a user requests a geometry which is not available in the
    inbuilt geometry library, an exception is raised.

    Returns
    -------
    None.

    """
    raise GeometryTemplateError()


def bluntcone_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = unpack_geometric_kwargs(**mesh_kw)
    length = geom["length"]
    angle = geom["cangle"]
    radius = geom["cradius"]
    major = geom["major_out"]
    minor = geom["minor_out"]
    i1loc = geom["i1_location"]
    if length is None:
        raise MeshingInputError("Length", "key not entered in the function\
 arguments")

    x, y = blunt_body_cone(length, angle, radius, major, minor,
                           Xi, Eta, iM, jM, clust_opt, beta, i1loc)
    return x, y


def bluntellip_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = unpack_geometric_kwargs(**mesh_kw)
    major1 = geom["major_out"]
    major2 = geom["major_in"]
    minor1 = geom["minor_out"]
    minor2 = geom["minor_in"]
    i1loc = geom["i1_location"]
    x, y = blunt_body_ellipt(major1, major2, minor1, minor2, Xi, Eta,
                             iM, jM, clust_opt, beta, i1loc)
    return x, y


def cgrid_mesh(clust_opt, Xi, Eta, dX, dY, iM, jM, config_obj, **mesh_kw):
    # clust_opt option be ON or OFF in this type of mesh.
    return


def duct_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    """Return a mesh in a duct with clustering along X or Y walls."""
    alpha, beta = unpack_clust_kwargs(**mesh_kw)
    alpha_x = alpha["alphaX"]
    alpha_y = alpha["alphaY"]
    beta_x = beta["betaX"]
    beta_y = beta["betaY"]
    # clust_opt option must always be ON in this type of mesh.
    if beta_x is None and beta_y is None:
        raise MeshingInputError("BetaX or BetaY", "value not entered")
    elif beta_x is not None and beta_y is not None:
        raise MeshingInputError("BetaX, BetaY", "only 1 must be entered.")
    geom = unpack_geometric_kwargs(**mesh_kw)
    if config_obj is not None:
        dX = config_obj.dX
        dY = config_obj.dY
        A = config_obj.Length
        B = config_obj.Height
    else:
        dX = geom["dx"]
        dY = geom["dy"]
        A = geom["length"]
        B = geom["height"]
    x = Xi.copy()
    y = Eta.copy()
    x, y = duct(A, B, x, y, dX, dY, Xi, Eta, iM, jM, alpha_x, alpha_y,
                beta_x, beta_y)

    return x, y


def cavity_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    """Return a mesh in a cavity with clustering at all walls."""
    alpha, beta = unpack_clust_kwargs(**mesh_kw)
    alpha_x = alpha["alphaX"]
    alpha_y = alpha["alphaY"]
    beta_x = beta["betaX"]
    beta_y = beta["betaY"]
    # clust_opt option must always be ON in this type of mesh.
    if beta_x is None or beta_y is None:
        raise MeshingInputError("BetaX, BetaY", "both must be entered")
    geom = unpack_geometric_kwargs(**mesh_kw)
    if config_obj is not None:
        A = config_obj.Length
        B = config_obj.Height
    else:
        A = geom["length"]
        B = geom["height"]
    if A is None:
        raise MeshingInputError("A:", "Domain length not entered")
    if B is None:
        raise MeshingInputError("B:", "Domain height not entered")
    x = Xi.copy()
    y = Eta.copy()
    x, y = cavity(A, B, x, y, iM, jM, Xi, Eta, alpha_x, alpha_y,
                  beta_x, beta_y)

    return x, y


def flatplate_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    # clust_opt option must always be ON in this type of mesh.
    """Return a mesh in a cavity with clustering at all walls."""
    _, beta = unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    # clust_opt option must always be ON in this type of mesh.
    if beta is None:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = unpack_geometric_kwargs(**mesh_kw)
    if config_obj is not None:
        dX = config_obj.dX
        dY = config_obj.dY
        A = config_obj.Length
        B = config_obj.Height
    else:
        dX = geom["dx"]
        dY = geom["dy"]
        A = geom["length"]
        B = geom["height"]
    x = Xi.copy()
    y = Eta.copy()
    x, y = flateplate(A, B, x, y, dX, dY, iM, jM, Xi, Eta, beta)

    return x, y


def o_cyl_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = unpack_geometric_kwargs(**mesh_kw)
    orad = geom["radius_out"]
    irad = geom["radius_in"]
    x = Xi.copy()
    y = Xi.copy()
    x, y = ogrid_cylinder(x, y, orad, irad, Xi, Eta, iM, jM,
                          clust_opt, beta)

    return x, y


def o_airfoil_mesh(clust_opt, Xi, Eta, iM, jM, config_obj, **mesh_kw):
    # clust_opt option may be ON or OFF in this type of mesh.
    _, beta = unpack_clust_kwargs(**mesh_kw)
    beta = beta["beta"]
    if beta is None and clust_opt is True:
        raise MeshingInputError("Beta", "clustering value not provided")
    geom = unpack_geometric_kwargs(**mesh_kw)
    ch = geom["chord"]
    thick = geom["thickness"]
    rad = geom["radius_out"]
    x = Xi.copy()
    y = Xi.copy()
    x, y = ogrid_airfoil(x, y, rad, ch, thick, Xi, Eta, iM, jM,
                         clust_opt, beta)

    return x, y


def wt_airfoil_mesh():
    return


def wt_cyl_mesh():
    return
