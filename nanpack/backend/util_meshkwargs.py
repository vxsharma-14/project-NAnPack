#   ***********************************************************************
#
#   FILE         utils_meshkwargs.py
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


def unpack_clust_kwargs(**grid_kw):
    alpha = grid_kw.get("Alpha", None)
    beta = grid_kw.get("Beta", None)
    alphaX = grid_kw.get("AlphaX", 0.5)
    alphaY = grid_kw.get("AlphaY", 0.5)
    betaX = grid_kw.get("BetaX", None)
    betaY = grid_kw.get("BetaY", None)
    alpha = _pack_alpha_into_dict(alphaX, alphaY, alpha)
    beta = _pack_beta_into_dict(betaX, betaY, beta)
    return alpha, beta


def unpack_geometric_kwargs(**grid_kw):
    lgt = grid_kw.get("Length", None)
    hgt = grid_kw.get("Height", None)
    dx = grid_kw.get("dX", None)
    dy = grid_kw.get("dY", None)
    chd = grid_kw.get("Chord", 1.0)
    th = grid_kw.get("Thickness", 0.2)
    orad = grid_kw.get("outRad", 3.0)
    irad = grid_kw.get("inRad", 1.0)
    cang = grid_kw.get("cAngle", 5.0)
    crad = grid_kw.get("cRadius", 1.0)
    omajor = grid_kw.get("outMajor", 4.0)
    ominor = grid_kw.get("outMinor", 3.5)
    imajor = grid_kw.get("inMajor", 3.0)
    iminor = grid_kw.get("inMinor", 2.0)
    i1loc = grid_kw.get("i1Location", "below-stagnation")
    geo = _pack_geom_into_dict(lgt, hgt, dx, dy, chd, th, orad, irad, cang,
                               crad, omajor, ominor, imajor, iminor, i1loc)
    return geo


def _pack_alpha_into_dict(alphaX, alphaY, alpha):
    """Return a dictionary for alpha."""
    alp = {
        "alphaX": alphaX,
        "alphaY": alphaY,
        "alpha": alpha
        }

    return alp


def _pack_beta_into_dict(betaX, betaY, beta):
    """Return a dictionary for beta."""
    bet = {
        "betaX": betaX,
        "betaY": betaY,
        "beta": beta
        }

    return bet


def _pack_geom_into_dict(lgt, hgt, dx, dy, ch, th, orad, irad, cang, crad,
                         omaj, ominr, imaj, iminr, i1loc):
    """Return a dictionary for geometric data."""
    geo = {
        "length": lgt,
        "height": hgt,
        "dx": dx,
        "dy": dy,
        "chord": ch,
        "thickness": th,
        "radius_out": orad,
        "radius_in": irad,
        "cangle": cang,
        "cradius": crad,
        "major_out": omaj,
        "minor_out": ominr,
        "major_in": imaj,
        "minor_in": iminr,
        "i1_location": i1loc
        }
    return geo
