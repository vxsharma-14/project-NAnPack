"""Not a public module."""
#   ***********************************************************************
#
#   FILE         boundary.py
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


def ReadBCfromFile(BCFileName):
    """Read user specified boundary conditions from the file.

    Call signature;

        ReadBCfromFile(BCFileName)

    Parameters
    ----------
    InFileName: list

        The string value representing the file to be read for
        simulation inputs.

    Inlet: list

        Boundary settings at the inlet boundary.

    Wall: list

        Boundary settings at the wall boundary.

    Farfld: list

        Boundary settings at the far-field boundary.

    Outlet: list

        Boundary settings at the outlet boundary.
    """
    import configparser
    from .checkconfig import CheckBCSections
    from .exceptions import InputFileError

    print("Reading boundary configurations from the file:")
    print(f'"{BCFileName}"')
    print('***********************************************************')

    config = configparser.ConfigParser()
    dataset = config.read(BCFileName)

    if dataset:
        print('SUCCESS: Boundary conditions configuration file parsing.')
    else:
        raise InputFileError("FileNotFound", BCFileName)

    # *********** CHECK IF ALL SECTIONS EXIST *********
    print('Checking whether all sections are included in config file.')
    CheckBCSections(config, BCFileName)

    # ************ INLET BC ***************
    InletAxis1 = config['INLET']['AXIS_1']
    Ain1 = float(config['INLET']['A1'])
    Bin1 = float(config['INLET']['B1'])
    InletAxis2 = config['INLET']['AXIS_2']
    if InletAxis2.lower() != 'none':
        Ain2 = float(config['INLET']['A2'])
        Bin2 = float(config['INLET']['B2'])
    else:
        Ain2 = 'none'
        Bin2 = 'none'
    InBCType = config['INLET']['BC_TYPE']
    Uin = float(config['INLET']['U'])
    print('Accessing boundary conditions at INLET: Completed.')

    # ************ WALL BC **************
    WallAxis1 = config['WALL']['AXIS_1']
    Aw1 = float(config['WALL']['A1'])
    Bw1 = float(config['WALL']['B1'])
    WallAxis2 = config['WALL']['AXIS_2']
    if WallAxis2.lower() != 'none':
        Aw2 = float(config['WALL']['A2'])
        Bw2 = float(config['WALL']['B2'])
    else:
        Aw2 = 'none'
        Bw2 = 'none'
    WallBCType = config['WALL']['BC_TYPE']
    Uw = float(config['WALL']['U'])
    print('Accessing boundary conditions at WALL: Completed.')

    # ************* FAR-FIELD BC ***************
    FarfldAxis1 = config['FAR-FIELD']['AXIS_1']
    Aff1 = float(config['FAR-FIELD']['A1'])
    Bff1 = float(config['FAR-FIELD']['B1'])
    FarfldAxis2 = config['FAR-FIELD']['AXIS_2']
    if FarfldAxis2.lower() != 'none':
        Aff2 = float(config['FAR-FIELD']['A2'])
        Bff2 = float(config['FAR-FIELD']['B2'])
    else:
        Aff2 = 'none'
        Bff2 = 'none'
    FarfldBCType = config['FAR-FIELD']['BC_TYPE']
    Uff = float(config['FAR-FIELD']['U'])
    print('Accessing boundary conditions at FAR-FIELD: Completed.')

    # *********** OUTLET BC **************
    OutletAxis1 = config['OUTLET']['AXIS_1']
    Ao1 = float(config['OUTLET']['A1'])
    Bo1 = float(config['OUTLET']['B1'])
    OutletAxis2 = config['OUTLET']['AXIS_2']
    if OutletAxis2.lower() != 'none':
        Ao2 = float(config['OUTLET']['A2'])
        Bo2 = float(config['OUTLET']['B2'])
    else:
        Ao2 = 'none'
        Bo2 = 'none'
    OutBCType = config['OUTLET']['BC_TYPE']
    Uo = float(config['OUTLET']['U'])
    print('Accessing boundary conditions at OUTLET: Completed.')

    print('Reading and acessing BC from config file: Completed.')
    Inlet = [InletAxis1, Ain1, Bin1, InletAxis2, Ain2, Bin2, InBCType, Uin]
    Wall = [WallAxis1, Aw1, Bw1, WallAxis2, Aw2, Bw2, WallBCType, Uw]
    Farfld = [FarfldAxis1, Aff1, Bff1, FarfldAxis2, Aff2, Bff2,
              FarfldBCType, Uff]
    Outlet = [OutletAxis1, Ao1, Bo1, OutletAxis2, Ao2, Bo2, OutBCType, Uo]

    print('***********************************************************')

    return Inlet, Wall, Farfld, Outlet


def BCatAxis(U, axis, Ubc, A, B, dX, dY):
    """Assign boundary conditions along the input axis."""
    if axis in ['X-lo', 'X-hi']:
        delta = dY
    elif axis in ['Y-lo', 'Y-hi']:
        delta = dX

    ijA = int(A/delta)
    ijB = int(B/delta)

    if axis == 'X-lo':
        U[0, ijA:ijB+1] = Ubc
    elif axis == 'X-hi':
        U[-1, ijA+1:ijB] = Ubc
    elif axis == 'Y-lo':
        U[ijA:ijB+1, 0] = Ubc
    elif axis == 'Y-hi':
        U[ijA+1:ijB, -1] = Ubc

    return U


def BC2D(U, BC, delX, delY):
    """Assign boundary condition at the walls, inlet and outlet.

    This function uses the conditions specified in the boundary
    configuration file which is provided as an argument to this function.

    Call signature;

        BC(U, BC, delX, delY)

    Parameters
    ----------
    U: 2D array

       Dependent variable.

    BC: list

        Boundary conditions settings obtained from boundary configuration
        file.

     delX: float

        Grid step size along X-axis.

    delY: float

        Grid step size along Y-axis.

    Returns
    -------
    U: 2D array

        Updated boundary conditions for the dependent variable.
    """
    Inlet, Wall, Farf, Outlet = BC
    InletAxis1, Ain1, Bin1, InletAxis2, Ain2, Bin2, InBCType, Uin = Inlet
    WallAxis1, Aw1, Bw1, WallAxis2, Aw2, Bw2, WallBCType, Uw = Wall
    FarfAxis1, Aff1, Bff1, FarfAxis2, Aff2, Bff2, FarfBCType, Uff = Farf
    OutletAxis1, Ao1, Bo1, OutletAxis2, Ao2, Bo2, OutBCType, Uo = Outlet

    # ---------------------------------------------------------------------
    #                         INLET BOUNDARY CONDITIONS
    # ---------------------------------------------------------------------
    U = BCatAxis(U, InletAxis1, Uin, Ain1, Bin1, delX, delY)
    if InletAxis2.lower() != 'none':
        U = BCatAxis(U, InletAxis2, Uin, Ain2, Bin2, delX, delY)
    # ---------------------------------------------------------------------
    #                         WALL BOUNDARY CONDITIONS
    # ---------------------------------------------------------------------
    U = BCatAxis(U, WallAxis1, Uw, Aw1, Bw1, delX, delY)
    if WallAxis2.lower() != 'none':
        U = BCatAxis(U, WallAxis2, Uw, Aw2, Bw2, delX, delY)
    # ---------------------------------------------------------------------
    #                     FAR-FIELD BOUNDARY CONDITIONS
    # ---------------------------------------------------------------------
    U = BCatAxis(U, FarfAxis1, Uff, Aff1, Bff1, delX, delY)
    if FarfAxis2.lower() != 'none':
        U = BCatAxis(U, FarfAxis2, Uff, Aff2, Bff2, delX, delY)
    # ---------------------------------------------------------------------
    #                          OUTLET CONDITIONS
    # ---------------------------------------------------------------------
    U = BCatAxis(U, OutletAxis1, Uo, Ao1, Bo1, delX, delY)
    if OutletAxis2.lower() != 'none':
        U = BCatAxis(U, OutletAxis2, Uo, Ao2, Bo2, delX, delY)

    return U
