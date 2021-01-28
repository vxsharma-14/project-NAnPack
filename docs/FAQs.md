# Frequently Asked Questions

### What is NAnPack?

NAnPack project is a collection of numerical solvers developed in Python for learning and teaching in scientific computing. It is an open-source package released under the MIT License.

This package is developed in a simple programming structure to promote ease-of-use and to accomadate all users within the scientific community having different coding skills.

This code has been developed by a researcher for researchers.

### What is included in the package?

Presently, the version 1.0.0-alpha1 (pre-release) of NAnPack is capable of computationally solving the 1D and 2D diffusion equation, wave equation, Poissons equation, and inviscid and viscous Bergers equation on structured, uniform grids using several finite difference methods. It includes the complete life-cycle stages for the numerical experimentation of a physical processes which includes pre-processing, running iterable/time-dependent simulations, and post-processing.

The modules for analytical solutions are also included in the package which can be used for initial benchmarking and validations which is a very important step in scientific computing.

This version also includes the module for grid transformation but it is not yet implemented within the numerical schemes.

*Please note that this is a pre-release version and things may not work as it is still under testing. If you come acros bugs, please report it through our Github repo.*

### What are the future plans for this package?

The first release version will include the finite difference numerical schemes for full Navier-Stokes equations on the curvilinear coordinate system that is capable of solving flow around airfoils.

If the project finds more contributors that can contribute to the development of this package, we are very ambitious in the plans for future versions. For example-  

* Regarding physical processes - In the subsequent versions, our aim is to include numerical methods for

    * Simulating electromagnetics applications by solving Maxwells equations,
    * Simulating chemically reacting hypersonic flow around airfoils.  
    

* Regarding numerical methods - We also aim to include the popular finite volume (FV) methods, immersed boundary method (IBM), unstructured grid generation procedures.

### Why NAnPack is build on Python and not C++ or Fortran?

**The primary objective of this project is to have a scientific computing package that can be easily used by everyone, be it a student, or a faculty member and that will ensure less time being spent on development and more time dedicated to research.**

Keeping this in mind it was decided to build NAnPack on Python language.
Here are some key factors that went into thinking:

* Firstly, its open-source.
* Python is easy to learn becuase of its simple syntax.
* Installation of Python and its components is super-easy.
* The pre-exisiting libraries and modules in Python save developers from "reinventing the wheel".
* All processing stages of scientific computations can be executed in a few lines of codes without requiring external third-party softwares/programs.   
For example, data visualization can be done using Python libraries without having to pay to use commercial visualization softwares. 
* Jupyter Notebook built for Python is a versatile browser based tool which makes the development, documentation, and code testing comparatively easy and requires less time.
* There are tons of tutorials, blog posts and support available for Python developers.

### Why use NAnPack and not other available commercial or open-source softwares?

* NAnPack is an open-source package, written in a simple programming structure that can be easily understood by anyone with minimum programming experience who are interested in learning and developing a strong background in scientific computing.
* Most academicians prefer Fortran and thus NAnPack will make the transition from Fortran to Python easy for them.
* Our codes can be used as a foundation to develop and run complex simulations.
* In our knowledge, there are no open-source softwares available that support learning for beginners in a way that we aim to do.