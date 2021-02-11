
def test_curvgrid:
    '''Test curvilinear grid generation routine.'''
    import nanpack.grid as grid
    x, y = grid.CurvilinearGrid(0.1, 21, 0.05, 41)


