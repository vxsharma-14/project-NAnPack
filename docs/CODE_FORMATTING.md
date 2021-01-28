# Code Formatting Guidelines

Please read these guidelines before coding a new function or a module. This is not an exhaustive list of rules. For other formatting guidlines refer to [Style Guide for Python](https://www.python.org/dev/peps/pep-0008/).

### General instructions

* *Any new module must include the docstring text at the top as mentioned in this [file](https://github.com/vxsharma-14/project-NAnPack/blob/main/docs/mod_top_text.txt). Please copy the text and edit only the file name, author name, and the version number*

* *Any new function must include the docstring text in the format as given in this [file](https://github.com/vxsharma-14/project-NAnPack/blob/main/docs/f_descript_format.txt). Please add function description, arguments, argument type, any default values, return parameters and any other detail that might be  helpful for the user. Please see other existing functions for the docstring example. All information must be included in the specified format.*

* *Restrict maximum characters in the code in each line to <75.*

* *If possible, break down large functions into smaller functions.*

* *Add comments whererver possible for the better understanding of the flow of the code; include any description or sources of equations with equation numbers in the comments.*

* *All other comments (single or multi-line) must start with a # and not triple quotes.*

* *Separate the multi-line codes in a logical and consistent manner. For example, break the codes before the math symbols (+,-,*,/) or at the paranthesis().*

### Naming Conventions

**Filename or Module name**  

*Use a logical module classifier and use the name in all small letters and without any symobols. *  

Note - If the function fits into one of the existing modules category, don't make a new module*

**Function Name**  

*Use a concise logical function name that can be related to the the task that function is performing.  Each word in the function name must be capitalized as in Camel case. For example, a function that solves a Poissons equation using Point-Gauss Seidel method can be named as* `PointGaussSeidel()`.

*Use abbreviated function names whenever possible, however, the function name must be explanatory.*

**Variable Name**

*Use a concise and logical variable name that are easily identifiable. Use word capitalization similar to function names. For example, a variable that is required to contain the output file name may be declared as* `OutFileName`.

*Single letter variables must be specified by small letters such as* `i`, `j`, `n`.

*Restrict the use of underscores in a function name or variable name to a minimum (strictly no hyphens at all).*




