#   ***********************************************************************
#
#   FILE         mesh1d.py
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

from .exceptions import MeshingInputError
from .util_clustering import calculate_beta1, calculate_eta1


def meshing_func(clust_loc, Xi, dX, iM, alpha, beta):
    """Return a meshing function requested by the user."""
    if beta is None:
        raise MeshingInputError("Beta", "clustering value not provided")
    clust_region = {
        "left": clustering_left_wall,
        "right": clustering_right_wall,
        "both": clustering_both_walls,
        "middle": clustering_mid_axis,
        }
    SelectedFunction = clust_region.get(clust_loc)

    return SelectedFunction(Xi, dX, iM, alpha, beta)


def clustering_left_wall(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = calculate_beta1(beta)
    for i in range(0, iM):
        # Grid clustering near X = 0.0
        xi1 = calculate_eta1(Xi[i], Xi[-1])
        x[i] = A * (((beta+1.0) - (beta-1.0)*(beta1**xi1))
                    / ((beta1**xi1)+1.0))
    return x


def clustering_right_wall(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = calculate_beta1(beta)
    alpha = 0.0
    for i in range(0, iM):
        # Grid clustering near X = A
        xi1 = calculate_eta1(Xi[i], Xi[-1], False)
        x[i] = get_clust_wall(xi1, alpha, beta, beta1, A)
    return x


def clustering_both_walls(Xi, dX, iM, alpha, beta):
    # clust_opt option must always be ON in this type of mesh.
    x = Xi.copy()
    A = dX * (iM-1)
    beta1 = calculate_beta1(beta)
    for i in range(0, iM):
        # Grid clustering equally distributed between X = 0 and X = A
        xi1 = calculate_eta1(Xi[i], Xi[-1], False)
        x[i] = get_clust_wall(xi1, alpha, beta, beta1, A)
    return x


def clustering_mid_axis(clust_opt, Xi, dX, iM):
    return


def get_clust_wall(metric1, alpha, beta, beta1, dom_size):
    A = dom_size
    met_e = (metric1-alpha) / (1.0-alpha)
    axis = A * (((2.0*alpha+beta)*(beta1**met_e)
                + 2.0*alpha - beta)
                / ((2.0*alpha+1.0)*(1.0+(beta1**met_e))))

    return axis
