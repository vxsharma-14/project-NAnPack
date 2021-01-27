## Installation

### For Windows OS

NAnPack has been designed with the intention to be used by anyone with the minimum computational skills to have a widespread user base. 

#### Dependencies

Since this package is written in Python 3 language, which may be download from the [Python homepage](https://www.python.org/downloads/). It comes with the integrated environment IDLE and Python Shell. The other popular development environment recommended for Python is [Jupyter Notebook](https://jupyter.org/) which is also open-source.

#### Installing useful Python packages

##### Check PIP

After you have downloaded the Python environment, we will use PIP to install packages/modules. PIP is the package manager for Python modules which is included by default with Python 3.4 or above. First check whether PIP is installed correctly by tyoing the following command in the command window. and enter

`C:\Users>pip --version`

The output should be similar to as shown below

`<pip 20.2.4 from c:\users\owner\appdata\local\programs\python\python37-32\lib\site-packages\pip (python 3.7)>`

If it does not work, check that the Python directory is included in your system environment PATH variable or re-install Python or try installing PIP.

##### Install NumPy package

Type the below command in your command window and enter. If the NumPy package is already installed in your system, its details will be displayed.

`pip show numpy`

If it is not pre-installed, type the below command and enter to install Numpy.

`pip install numpy`

###### List of required packages

Following the instructions mentioned above for installing Numpy package and intall other required packages given below.

*matplotlib* - for plotting simulation data  
`pip install matplotlib`

*math* - required for evaluating trignotmetric functions  
`pip install math`





