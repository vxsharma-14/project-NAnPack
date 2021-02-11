#   ***********************************************************************
#
#   FILE         _plots.py
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

def _Plot1D(dataFiles, uAxis, legend, markers, useFileCol,
            title, xlbl, ylbl):
    """Plot results along 1D axis."""
    import matplotlib.pyplot as plt
    import numpy as np

    print("Preparing data to plot results...")
    countFiles = len(dataFiles)

    # Assign font family and create axis for plotting
    plt.rc('font', family='sans-serif', size=10)
    fig, ax = plt.subplots(dpi=150)

    data = {}  # Initialize a dictionary

    # Use dictionary to create multiple line plots from saved files
    for i in range(countFiles):
        data[f"x{i}"] = np.loadtxt(dataFiles[i], unpack=True, skiprows=3,
                                   usecols=0)
        data[f"u{i}"] = np.loadtxt(dataFiles[i], unpack=True, skiprows=3,
                                   usecols=useFileCol)

    # Start plotting
    print("Plotting 1D results")
    for i in range(countFiles):
        if markers is not None:
            mrkr = markers[i]
        else:
            mrkr = None
        if uAxis == "Y":
            ax.plot(data[f"x{i}"], data[f"u{i}"], linewidth=0.5,
                    label=legend[i], marker=mrkr, markersize=3,
                    markevery=3)
        elif uAxis == "X":
            ax.plot(data[f"u{i}"], data[f"x{i}"], linewidth=0.5,
                    label=legend[i], marker=mrkr, markersize=3,
                    markevery=3)
    # Format and customize plot
    plt.grid(which='both', axis='both', color='lightgrey', linestyle=':',
             linewidth=0.5)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.title(title, fontsize=10)
    ax.legend(fontsize=9)
    # plt.tight_layout()
    plt.show()
