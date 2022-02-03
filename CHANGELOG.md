# Version updates.

### 1.0.0-alpha5
**Patches**  
- Fixed docstring text-- Call Signature in `.postprocess.Plot1DResults` function.  
- Changed function arguments in `.postprocess.Plot1DResults()` - added `dataFiles` parameter to required argument, 
removed `dataFiles` from kwargs.  
- Changed `POISSONS` to `LAPLACE` in equation model definitions wherever required. Entering `POISSONS` as the 
  simulation model now will produce Exception error.
- Changed the way preprocessing is done. Users will not be effected by this change.

**Minor**  
- Function `.preprocess.DiffusionNumbers()`, instead of returning fixed two values as previously, now returns one 
value for 1D and two values for 2D applications.
- Function `.meshing.RectangularMesh()`, instead of returning fixed two values as previously, now returns one value 
  for 1D and two values for 2D applications.

**Major**   
- Function names changes -  
  - `.grid.ComputeGridPoints` to `.meshing.CalcGridPoints`
  - `.grid.ComputeGridSteps` to `.meshing.CalcGridStepsize`
  - `.grid.RectangularGrid` to `.meshing.RectangularMesh`
  - Module name `.grid` is changed to `.meshing`.
  - Module name `.secondaryfunctions` is changed to `nanpack.utils`.
- Changed paramters in `.postprocess.WriteSolutionToFile`  
      Previously allowed paramters: U, n, nWrite, nMax, OutFileName, dX, dY  
      Currently allowed: n, U, CfgClsObj, **kwargs
- Changed paramters in `.postprocess.WriteConvHistToFile`  
      Previously allowed paramters: CfgClsObj, n, Error, HistFName=None  
      Currently allowed: n, Error, CfgClsObj, **kwargs
- Changed paramters in `.postprocess.WriteSolutionIn1DFormat`  
      Previously allowed paramters: CfgClsObj, U  
      Currently allowed: U, CfgClsObj, **kwargs  
- New features
  - Added `.postprocess.Plot2DResults` function for 2D plotting - single plot or multiple plots
  - Added module `.mixedpdesolvers` for solution of viscous Burgers equation.
  - Added features to solve non-dimensionalized form of the model equations.
  - Added `.preprocess.NondimensionalizeTime()` module.
  - Added `.preprocess.NondimensionalizeGrid()` module.
  - Added function argument `Scaling` in function `.preprocess.DiffusionNumbers()` to incorporate non-dimensionalization of the diffusion numbers.
  - Added `.postprocess.DimensionalizeSolution` module.
  - Added functions `.meshing.StructuredMesh` to generate 2D non-uniform mesh in a physical coordinate system (x, y).
  - Added functions `.meshing.StructuredMesh1D` to generate 1D non-uniform mesh in a physical coordinate system (x).
  - Added function `,meshing.CalcMeshMetrics` to calculate metrics and Jacobian of the transformation values from 
    physical mesh data.
  - Added function `,meshing.CalcMeshMetrics1D` to calculate metrics and Jacobian of the transformation values from 
    physical 1D mesh data.
  - Added function `.meshing.SaveMeshMetrics` to write metrics data to file.
  - Added function `.meshing.SaveMeshMetrics1D` to write 1D metrics data to file.
  - Added function `.meshing.PlotMesh` to plot 2D mesh.
  - Added function `.meshing.PlotMeshMetrics` to plot metrics data.
  
### 1.0.0-alpha4
made changes to the modules- added/deleted//merged modules
added toutorials
added function documentations etc.
This is not an accurate update log because a lot of changes were made in the overall structure of the package and 
since this is a pre-release version, not all changes are recorded.

### 1.0.0-alpha3
made changes to the setup.py file to correctly package the components.

### version 1.0.0-alpha2
1. setup.py file not correctly configured in the alpha1 version.
2. edited/changed documentations