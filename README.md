# FPKriSten
## Description
`FPKriSten` is a tool for computation of roundoff error bounds using Bernstein expansion.

If the roundoff error of a program implementing a polynomail function f(x)  is given by r(x,e) s.t.:

			r(x,e) = l(x,e)+h(x,e)

with `l(x,e)` the part of `r(x,e)` linear w.r.t `e` and `h(x,e)` the part  of `r(x,e)` non-linear in `e`.

Then `FPKriSten` gives an upper bound of `|l(x,e)|` for `x` laying inside a box, and `e` enclose by a given epsilon. 

Thus, the semantic of the handled programs is:

- polynomial functions taken into isolation

ex: x1^3 + (3/4)*x1*x2^2
- inputs laying inside a box

ex: x1 = [-1 1] and x2 = [-2 2]

### Programs representation
FPBern(a) handles input files with .ini extension (this is mandatory) with the following structure:

 + OPTIONS
- name = `name of the program`
- precision = `machine precision as 2^(-precision)`
- nbvars = `dimension of x`
- nberrors = `dimension of e`
 + Programs
- function = ` polynomials function with the operation +,*,- and ^ (and / in the coefficients)`
- input_bl = `lower bounds of the input values`
- input_bu = `upper bounds of the input values`

## Installation instructions
### Prerequisites
FPBern(a) is implemented in C++. Thus, a C++ compiler is required.
FPBern(a) was tested on Ubuntu 14.04 LTS.

Moreover, FPBern(a) relies on three external libraries:
- GiNaC (GiNaC is Not a CAS), for the symbolic manipulation of polynomials
- GLPK (GNU Linear Programming Kit), for solving linear programming problems
- boost (boost/property_tree/ptree.hpp and boost/property_tree/ini_parser.hpp and boost/lexical_cast.hpp)

###Download
FPBern(a) is maintained as a GitHub repository at the address https://github.com/roccaa/FPBern.git

It can be obtained either by typing the shell command:

$ git clone https://github.com/roccaa/FPBern.git

or by downloading the ZIP archive at https://github.com/roccaa/FPBern.git

###Installation
To install from the source type:

	$ make

This creates a binary called FPBern in /bin

To run FPBern, move to /bin and launch the binary with the command:

	$ ./FPBern file_name1 file_name2 ...    (without the .ini)

To run a set of benchmarks, launch the command

 	$ ./FPBern rigidbody1 rigidbody2 kepler0 kepler1 kepler2 sineTaylor sineOrder3 sqroot himmilbeau





