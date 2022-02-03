"""Not a public module."""
#   ***********************************************************************
#
#   FILE         mesh_otype.py
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

import math
from .util_clustering import clustering_parameter


def ogrid_airfoil(x, y, radius, chord, thick, Xi, Eta, iM, jM, clust_opt,
                  beta):
    if iM % 2 == 0:
        imid = int(iM/2 + 1)
    else:
        imid = int((iM-1)/2 + 1)
    # distance between airfoil leading edge and # the outer boundary.
    del_d = radius - chord/2.0

    # Define airfoil surface
    dx_afc = chord / ((iM+1)/2-1)  # grid steps size on the airfoil chord
    for i in range(0, imid):
        x[i][0] = del_d + dx_afc*(imid-i-1)
        x_af = x[i][0] - del_d  # x-locations on airfoil surface
        xaf2 = x_af * x_af
        xaf3 = x_af * x_af * x_af
        y[i][0] = -(thick/0.2) * (0.2969 * math.sqrt(x_af)
                                  - 0.126 * x_af
                                  - 0.3516 * xaf2
                                  + 0.2843 * xaf3
                                  - 0.1015 * xaf2 * xaf2)

    for i in range(imid, iM-1):
        x[i][0] = del_d + dx_afc*(i+1-imid)
        x_af = x[i][0] - del_d
        xaf2 = x_af * x_af
        xaf3 = x_af * x_af * x_af
        y[i][0] = (thick/0.2) * (0.2969 * math.sqrt(x_af)
                                 - 0.126 * x_af
                                 - 0.3516 * xaf2
                                 + 0.2843 * xaf3
                                 - 0.1015 * xaf2 * xaf2)
    x[-1][0] = x[0][0]
    y[-1][0] = y[0][0]

    # Define outer circular boundary
    dtheta = 2.0 * math.pi / (iM-1)  # calculate angular step size
    for i in range(0, imid):
        x[i][-1] = radius + radius*math.cos(dtheta*(i))
        y[i][-1] = -radius*math.sin(dtheta*(i))

    for i in range(imid, iM):
        x[i][-1] = radius + radius*math.cos(dtheta*(iM-i-1))
        y[i][-1] = radius*math.sin(dtheta*(iM-i-1))

    x, y = interior_of_ogrid(iM, jM, imid, x, y, clust_opt, Eta, beta)
    return x, y


def ogrid_cylinder(x, y, rad_o, rad_i, Xi, Eta, iM, jM, clust_opt, beta):
    if iM % 2 == 0:
        imid = int(iM/2 + 1)
    else:
        imid = int((iM-1)/2 + 1)
    dtheta = 2.0 * math.pi / (iM-1)  # calculate angular step size
    for i in range(0, imid):
        x[i][0] = rad_o + rad_i*math.cos(dtheta*(i))
        y[i][0] = - rad_i*math.sin(dtheta*(i))

    for i in range(imid, iM-1):
        x[i][0] = rad_o + rad_i*math.cos(dtheta*(iM-i-1))
        y[i][0] = rad_i*math.sin(dtheta*(iM-i-1))
    x[-1][0] = x[0][0]
    y[-1][0] = y[0][0]

    # Define outer circular boundary
    for i in range(0, imid):
        x[i][-1] = rad_o + rad_o*math.cos(dtheta*(i))
        y[i][-1] = -rad_o*math.sin(dtheta*(i))

    for i in range(imid, iM):
        x[i][-1] = rad_o + rad_o*math.cos(dtheta*(iM-i-1))
        y[i][-1] = rad_o*math.sin(dtheta*(iM-i-1))

    x, y = interior_of_ogrid(iM, jM, imid, x, y, clust_opt, Eta, beta)

    return x, y


def interior_of_ogrid(iM, jM, imid, x, y, copt, Eta, beta):
    """Return x, y at all interior grid points."""
    S = Eta.copy()
    # Calculate the distance, delI between j = jM and j = 1
    # Calculate the angle alpI of each line i from j = 1 to j = j jM
    delI = []
    alpha = []
    for i in range(0, iM):
        delxx = x[i, 0] - x[i, -1]
        delyy = y[i, -1] - y[i, 0]
        delI.append(math.sqrt(delxx*delxx + delyy*delyy))
        alpha.append(math.atan2(abs(delyy), abs(delxx)))
    # Calculate interior grid point locations using
    #  clustering option = True or False.
    for i in range(0, iM):
        for j in range(0, jM):
            if copt is True:
                S[i][j] = clustering_parameter(Eta[i, j], Eta[i, -1],
                                               beta, delI[i])
            else:
                S[i][j] = delI[i]*(j)/(jM-1)

    if imid % 2 == 0:
        i_lower_mid = int(imid/2)
        i_upper_mid = imid + int((iM-imid-1)/2)
    else:
        i_lower_mid = int((imid+1)/2)
        i_upper_mid = imid + int((iM-imid)/2)

    for j in range(1, jM-1):
        for i in range(0, iM-1):
            if i < i_lower_mid:
                x[i][j] = x[i][0] + S[i][j]*math.cos(alpha[i])
                y[i][j] = y[i][0] - S[i][j]*math.sin(alpha[i])
            elif i_lower_mid <= i < imid:
                x[i][j] = x[i][0] - S[i][j]*math.cos(alpha[i])
                y[i][j] = y[i][0] - S[i][j]*math.sin(alpha[i])
            elif imid <= i < i_upper_mid:
                x[i][j] = x[i][0] - S[i][j]*math.cos(alpha[i])
                y[i][j] = y[i][0] + S[i][j]*math.sin(alpha[i])
            elif i_upper_mid <= i < iM:
                x[i][j] = x[i][0] + S[i][j]*math.cos(alpha[i])
                y[i][j] = y[i][0] + S[i][j]*math.sin(alpha[i])
        x[-1][j] = x[0][j]
        y[-1][j] = y[0][j]

    return x, y
