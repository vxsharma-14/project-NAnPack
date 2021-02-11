# Installation (v1.0.0-alpha4)

## For Windows OS

### I. Requirements

The package requires Python 3, which may be downloaded from the [Python homepage](https://www.python.org/downloads/). It comes with the integrated environment IDLE and Python Shell. The other popular development environment recommended for Python is [Jupyter Notebook](https://jupyter.org/) which is also open-source.

### II. Approach 1 - Installing NAnPack from source

One straightforward way of installing without having to use Git is downloading zip files from the GitHub repository.
1. Visit GitHub project page, [link here](https://github.com/vxsharma-14/project-NAnPack).
2. Download ZIP in the target directory `/path/to/myproject` and unzip the contents.
3. On your terminal/command window, change the directory to the root of the unzipped NAnPack installation directory where "setup.py" is located `cd /path/to/myproject`.
4. Install nanpack using the following command.  
`python setup.py install`

This process will ensure that that you have downloaded the required configuration files located in the `./input/` folder.

### III. Approach 2 - Installing NAnPack using PIP

NAnPack is uploaded on Python Package Index (PyPI) repository and thus it can be easily installed by entering the following on your terminal:

    pip install nanpack

To get the configuration files, you may have to dig into the directory where all python packages are installed and copy the input folder to the target directory for your projects.

If you don't have PIP installed, first read 'Check PIP' section and then continue from here.

### IV. Check Installation

To check correct installation of the package, run the following tests on your command window/terminal

1. Test nanpack installation - `python -m nanpack.tests.test_nanpackinstall`
2. Test required third-party packages - `python -m nanpack.tests.test_thirdpartyinstalls`. If this test fails, proceed to the section V.
3. Run an example case to test everything is working - `python -m nanpack.tests.test_run`.

The detailed outputs from these tests can be [found here](/running-tests.ipynb).

*If Test#2 is passed, skip the below sections V and VI.*

### V. Installing other Python packages 

Following are the required additional third-party packages to ensure correct functioning of NAnPack - NumPy and Matplotlib.

Whether or not these packages are installed on your system can be checked by entering on the command window:

    pip show <package-name>

If they are not already installed, type the following in the command window:

    pip install numpy  
    pip install matplotlib

### VI. Check PIP

After you have downloaded the Python environment, we will use PIP to install packages/modules. PIP is the package manager for Python modules which is included by default with Python 3.4 or above. First check whether PIP is installed correctly by typing the following command in the command window and enter

    C:\Users>pip --version

The output should be similar to as shown below

    <pip 20.2.4 from c:\users\owner\appdata\local\programs\python\python37-32\lib\site-packages\pip (python 3.7)>

If it does not work, check that the Python directory is included in your system environment PATH variable or re-install Python or try installing PIP.






