#**************************************************************************
def BC2D(U, BC, delX, delY):
    '''Assigns boundary condition at the walls, inlet and outlet using the
    conditions specified in the boundary configuration file.

    Call signature;

        BC(U, BC, delX, delY)

    Parameters
    ----------

    U: 2D array

       Dependent variable.

    BC: boundary values obtained from boundary configuration file.
    '''    
    Inlet, Wall, Farf, Outlet = BC
    print('Proceeding to assign boundary conditions.')
    InletAxis1, Ain1, Bin1, InletAxis2, Ain2, Bin2, InBCType, Uin = Inlet
    WallAxis1, Aw1, Bw1, WallAxis2, Aw2, Bw2, WallBCType, Uw = Wall
    FarfAxis1, Aff1, Bff1, FarfAxis2, Aff2, Bff2, FarfBCType, Uff = Farf
    OutletAxis1, Ao1, Bo1, OutletAxis2, Ao2, Bo2, OutBCType, Uo = Outlet

    #----------------------------------------------------------------------
    #                         INLET BOUNDARY CONDITIONS
    #----------------------------------------------------------------------
    U = BCatAxis(U, InletAxis1, Uin, Ain1, Bin1, delX, delY)
    print(f'Inlet boundary conditions assigned at axis: {InletAxis1}.')
    if InletAxis2.lower() != 'none':
        U = BCatAxis(U, InletAxis2, Uin, Ain2, Bin2, delX, delY)
        print(f'Inlet boundary conditions assigned at axis: {InletAxis2}.')
    #----------------------------------------------------------------------
    #                         WALL BOUNDARY CONDITIONS
    #----------------------------------------------------------------------
    U = BCatAxis(U, WallAxis1, Uw, Aw1, Bw1, delX, delY)
    print(f'Wall boundary conditions assigned at axis: {WallAxis1}.')
    if WallAxis2.lower() != 'none':
        U = BCatAxis(U, WallAxis2, Uw, Aw2, Bw2, delX, delY)
        print(f'Wall boundary conditions assigned at axis: {WallAxis2}.')
    #----------------------------------------------------------------------
    #                     FAR-FIELD BOUNDARY CONDITIONS
    #----------------------------------------------------------------------
    U = BCatAxis(U, FarfAxis1, Uff, Aff1, Bff1, delX, delY)
    print(f'Far-field boundary conditions assigned at axis: {FarfAxis1}.')
    if FarfAxis2.lower() != 'none':
        U = BCatAxis(U, FarfldAxis2, Uff, Aff2, Bff2, delX, delY)
        print(f'Far-field boundary conditions assigned at axis: {FarfAxis2}.')
    #----------------------------------------------------------------------
    #                          OUTLET CONDITIONS
    #----------------------------------------------------------------------
    U = BCatAxis(U, OutletAxis1, Uo, Ao1, Bo1, delX, delY)
    print(f'Outlet boundary conditions assigned at axis: {OutletAxis1}.')
    if OutletAxis2.lower() != 'none':
        U = BCatAxis(U, OutletAxis2, Uo, Ao2, Bo2, delX, delY)
        print(f'Outlet boundary conditions assigned at axis: {OutletAxis2}.')
    
    print('Boundary configuration: Completed.')

    return U

#**************************************************************************
def BCatAxis(U, axis, Ubc, A, B, dX, dY):
    '''Configure boundary conditions at inlet, outlet, wall and
    far-field boundaries.
    '''
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
