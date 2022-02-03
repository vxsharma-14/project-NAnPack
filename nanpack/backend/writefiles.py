#   ***********************************************************************
#
#   FILE         writefiles.py
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

from .util_unpackkw import writesolution_kwargs, write1dsolution_kwargs


def save_metrics_2d(x, y, xix, xiy, etax, etay, jj, f_name):
    iM, jM = x.shape
    outf = open(f_name, "w")
    print(f"IMAX= {iM}, JMAX= {jM}", file=outf)
    print(f'{"X":^11} {"Y":^11} {"XiX":^11} {"XiY":^11} {"EtaX":^11}\
 {"EtaY":^11} {"JJ":^11}', file=outf)
    for i in range(0, iM):
        for j in range(0, jM):
            print(f"{x[i][j]:>11.6f} {y[i][j]:>11.6f} {xix[i][j]:>11.6f}\
 {xiy[i][j]:>11.6f} {etax[i][j]:>11.6f} {etay[i][j]:>11.6f}\
 {jj[i][j]:>11.6f}", file=outf)

    outf.close()


def save_metrics_1d(x, xix, jj, f_name):
    iM, = x.shape
    outf = open(f_name, "w")
    print(f"IMAX= {iM}", file=outf)
    print(f'{"X":^11} {"XiX":^11} {"JJ":^11}', file=outf)
    for i in range(0, iM):
        print(f"{x[i]:>11.6f} {xix[i]:>11.6f} {jj[i]:>11.6f}",
              file=outf)

    outf.close()


def save_solution(CfgClsObj, n, U, **kwargs):
    """Write 1D or 2D results to a specified file."""
    if CfgClsObj is not None:
        nMax = CfgClsObj.nMax
        nWrite = CfgClsObj.nWrite
        dX = CfgClsObj.dX
        dY = CfgClsObj.dY
        fName = kwargs.get("FileName", None)
        if fName is None:
            fName = CfgClsObj.OutFileName
    else:
        kw = writesolution_kwargs(**kwargs)
        nMax = kw["nMax"]
        nWrite = kw["nWrite"]
        fName = kw["FileName"]
        dX = kw["dX"]
        dY = kw["dY"]

    if n % nWrite == 0 or n == nMax:
        OutFile = open(fName, "w")
        shapeU = U.shape  # Obtain shape for Dimension
        if len(shapeU) == 1:  # Dimension = 1D
            iM, = shapeU
            save_solution_1d(iM, dX, U, OutFile)
        elif len(shapeU) == 2:  # Dimension = 2D
            iM, jM = shapeU
            save_solution_2d(iM, jM, dX, dY, U, OutFile)
        OutFile.close()


def save_solution_2d(iMax, jMax, dX, dY, U, OutFile):
    """Save the 2d solution to a specified file."""
    print(f"IMAX= {iMax}, JMAX= {jMax}", file=OutFile)
    print(f"Grid Steps= {dX:>6.4f} {dY:>6.4f}",
          file=OutFile)
    print(f'{"X":^14} {"Y":^14} {"U":^14}', file=OutFile)
    for i in range(0, iMax):
        for j in range(0, jMax):
            print(f"{dX*i:>12.8e} {dY*j:>12.8e} {U[i][j]:>12.8e}",
                  file=OutFile)


def save_solution_1d(iMax, dX, U, OutFile):
    """Save the 1d solution to a specified file."""
    print(f"IMAX= {iMax}", file=OutFile)
    print(f"Grid Step= {dX:>6.4f}", file=OutFile)
    print(f'{"X":^14} {"U":^14}', file=OutFile)
    for i in range(0, iMax):
        print(f"{dX*i:>12.8e} {U[i]:>12.8e}", file=OutFile)


def save_2d_solution_as_1d(CfgClsObj, U, **kwargs):
    """Write 2d solution in 1d format to the specified file."""
    if CfgClsObj is not None:
        dX = CfgClsObj.dX
        dY = CfgClsObj.dY
        nodes = CfgClsObj.nodes  # locations on the axis
        axis = CfgClsObj.PrintNodesDir  # coordinate axis
        fName = kwargs.get("FileName1D", None)
        if fName is None:
            fName = CfgClsObj.Out1DFName
    else:
        kw = write1dsolution_kwargs(**kwargs)
        fName = kw["FileName"]
        dX = kw["dX"]
        dY = kw["dY"]
        nodes = kw["ax_loc"]
        axis = kw["axis"]
    shapeU = U.shape
    iMax, jMax = shapeU
    line = []  # initialize a var to store linewise data
    # Format line # 1 in output file which print location points
    locations = ["{0:^12.2f}".format(li) for li in nodes]
    # Format line # 2 in output file which prints dashes
    dashes = ["{0:^12}".format('------') for li in nodes]
    OutFile = open(fName, "w")
    if axis.upper() == "X":
        save_2d_along_x(U, dX, dY, jMax, nodes, locations, dashes,
                        line, OutFile)

    elif axis.upper() == "Y":
        save_2d_along_y(U, dX, dY, iMax, nodes, locations, dashes,
                        line, OutFile)
    OutFile.close()

    return


def save_2d_along_x(U, dX, dY, jM, nodes, locations, dashes, line, OutF):
    """Save 2d results at various x nodes as a function of y."""
    iX = [int(i/dX) for i in nodes]
    print(f'{"Y":<5} {"X=":^2}', *locations, sep=' ', file=OutF)
    print(f'{"---":<8}', *dashes, sep=' ', file=OutF)
    for j in range(0, jM):
        del line[:]  # Delete elements in the list for next Y location
        for i in iX:
            # Read output at different X locations at a constant Y and
            # store in one line
            line.append(U[i][j])
            # Convert "line" list elements to float & format the output
            lines = ["{0:12.8f}".format(float(li)) for li in (line)]
        print(f"{dY*j:<8.4f}", *lines, sep=' ', file=OutF)


def save_2d_along_y(U, dX, dY, iM, nodes, locations, dashes, line, OutF):
    """Save 2d results at various y nodes as a function of x."""
    jY = [int(j/dY) for j in nodes]
    print(f'{"X":<5} {"Y=":^2}', *locations, sep=' ', file=OutF)
    print(f'{"---":<8}', *dashes, sep=' ', file=OutF)
    for i in range(0, iM):
        del line[:]  # Delete elements in the list for next X location
        for j in jY:
            # Read output at different Y locations at a constant X and
            # store in one line
            line.append(U[i][j])
            # Convert "line" list elements to float & format the output
            lines = ["{0:12.8f}".format(float(li)) for li in (line)]
        print(f"{dX*i:<8.4f}", *lines, sep=' ', file=OutF)


def save_convergence_hist(CfgClsObj, n, error, **kwargs):
    """Write 1D or 2D results to a specified file."""
    from .util_unpackkw import writeconv_kwargs
    if CfgClsObj is not None:
        nDisplay = CfgClsObj.nDisplay
        nMax = CfgClsObj.nMax
        State = CfgClsObj.State
        ConvCrit = CfgClsObj.ConvCrit
        totTime = CfgClsObj.totTime
        fName = kwargs.get("HistFName", None)
        if fName is None:
            fName = CfgClsObj.HistFileName
    else:
        kw = writeconv_kwargs(**kwargs)
        nDisplay = kw["nDisplay"]
        fName = kw["FileName"]
        nMax = kw["nMax"]
        State = kw["State"]
        ConvCrit = kw["ConvCriteria"]
        totTime = kw["SimTime"]
    # Write convergence data to file
    if n == 0 or n == 1:
        HistFile = open(fName, "w")
        print(f'{"ITER":>7} {"ERROR":>15}', file=HistFile)
        print(f'{"----":>7} {"-----":>15}', file=HistFile)
    if n % nDisplay == 0 or n == nMax or error < ConvCrit:
        HistFile = open(fName, "a")
        print(f'{n:>7} {error:>15.8f}', file=HistFile)
        HistFile.close()
        write_conv_msg(n, nMax, error, fName, State, ConvCrit, totTime)


def write_conv_msg(n, nMax, Error, fName, State, ConvCrit, totTime):
    """Write convergence status message at the end of conv. file."""
    if State.upper() == "STEADY":
        if n == nMax and Error > ConvCrit:
            msg = write_steadyst_notconv_msg(nMax)
            save_msg(msg, fName)
        elif Error < ConvCrit:
            msg = write_steadyst_conv_msg(n, ConvCrit)
            save_msg(msg, fName)
    elif State.upper() == "TRANSIENT" and n == nMax:
        msg = write_transient_conv_msg(nMax, totTime)
        save_msg(msg, fName)


def save_msg(msg, fileName):
    """Save convergence message to the file."""
    HistFile = open(fileName, "a")
    print(msg, file=HistFile)
    print("Writing convergence log file: Completed.")
    print("Files saved:")
    print(f'"{fileName}".')
    HistFile.close()


def write_steadyst_notconv_msg(nMax):
    """Return the convergence status message for writing to file."""
    PrintMsg = f"\nSTATUS: NOT CONVERGED\nMAX. ITERATIONS={nMax}"
    print(PrintMsg)
    print()
    return PrintMsg


def write_steadyst_conv_msg(n, ConvCrit):
    """Return the convergence status message for writing to file."""
    PrintMsg = f"\nSTATUS: CONVERGED SOLUTION OBTAINED AT\
\nITERATIONS={n},\nMAX. ERROR= {ConvCrit}"
    print(PrintMsg)
    print()
    return PrintMsg


def write_transient_conv_msg(nMax, totTime):
    """Return the convergence status message for writing to file."""
    PrintMsg = f"\nSTATUS: SOLUTION OBTAINED AT\nTIME LEVEL=\
 {totTime} s.\nTIME STEPS= {nMax}"
    print(PrintMsg)
    print()
    return PrintMsg
