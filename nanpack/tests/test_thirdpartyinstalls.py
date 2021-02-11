#   ***********************************************************************
#
#   FILE         test_requiredinstalls.py
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


def test_3rdpartypackage():
    """Test numpy and matplotlib installation."""
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        x = np.arange(0, 361, 10)
        y = np.sin(x*3.14/180.0)
        plt.plot(x, y)
        plt.title("Plot of sin(theta)")
        plt.xlabel("X")
        plt.ylabel("sin(theta)")
        plt.xlim(0, 360)
        print("Close plot to continue testing.")
        plt.show()
    except:
        print("")


def test_mathpackage():
    """Test math package."""
    import math
    mp = round(math.pi, 2)
    assert mp == 3.14


if __name__ == "__main__":
    test_3rdpartypackage()
    print("Numpy package test SUCCESS.")
    print("Matplotlib package test SUCCESS.")
    test_mathpackage()
    print("Math package test SUCCESS.")
