# CODE FORMAT GUIDELINES

Please read these guidelines before coding a new function or a module.

### General instructions

* *Any new module must include the docstring text as given in this [file](moduletoptext.txt) on top of that module. Please copy the text and edit only the file name, author name, and the version number*

* *Any new function must include the docstring text in the format giver in this [file](FunctionDescriptionFormat.txt). Please add function description, arguments, argument type, any default values and any other detail that might be  helpful for the user. Please see other existing function for the docstring example. All information must be included in the specified format.*

* *Restrict maximum characters in the code in each line to <75.*

* *If possible, break down large functions into smaller functions.*

* *Add comments whererver possible for better understanding of the flow of the code, include any description or sources of equations with equation numbers in the comments.*

* *All other comments (single or multi-line) must start with a # and not triple quotes.*

* *Separate the multi-line codes in a logical and consistent manner. For example, break the codes at the math symbols (+,-,*,/) or at the paranthesis().*

### Naming Conventions

**Filename or Module name**  

*Use a logical module classifier and use the classifier as the module name in all small letters and without any symobols. *  

Note - If the function fits into one of the existing modules category, don't make a new module*

**Function Name**  

*Use a concise logical function name that can be related to the the task that function is performing.  Each word in the function names must be Capitalized. For example, a function that solves a Poissons equation using Point-Gauss Seidel method can be named as* `PointGaussSeidel()`.

*Use abbreviated function names whenever possible, however, the function name must be explanatory.*

**Variable Name**

*Use a concise and logical variable name that are easily identifiable. Use word Capitalization similar to function names. For example, a variable that is required to contain the output file name may be declared as* `OutFileName`.

*Single letter variables must be specified by small letters such as* `i`, `j`, `n`.

*Restrict the use of underscores in a function name or variable name to a minimum (strictly no hyphens at all).*




