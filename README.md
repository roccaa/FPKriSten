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
`FPKriSten` is use in Matlab script files. The script is separated in two phases: 

+ The building and parsing phase (done on separate script files for the benchmarks)
+ The computation of the roundoff error.

PHASE1: Building

(1) declare the system:

- `complex_sparse = 0;` , currently not used but necessary to avoid error
- `nvar` = number of variable `x`;
- `name` = put the name of the program;
- `nparam` = number of variable `e`;
- `[str,vars]  = build_sdpvar(nvar,nparam);eval(str);vars = eval(vars);` build the variables with Yalmip parser
- `q` = here put the `l(x,e)` program with the symbols `x_i` with `i=1,...,nvar` for `x`, and `x_j` with `j = nvar+1,...,nvar+nparam`  for `e`;
- interval=[[one line for each interval of `x` and `e` (in order)]];
- `qsdp = box_norm(q,vars,interval);` compute the normalization
- `[powers,coefficients] = getexponentbase(qsdp,vars);` polynomial parsing with Yalmip
- `p = str2double(sdisplay(powers));` power matrix
- `c = str2double(sdisplay(coefficients));` coefficients matrix
- `n = nvar+nparam;` total number of variables `(x,e)`
- `[I,J] = build_box_sparcity(nvar,nparam);`I and J sparcity pattern specific to this problem
- `system_info = [n complex_sparse];` build an option arry (mandatory)
- `G = create_unitBox(n);` build a unit box (mandatory)

Build the path and the models files (if you want):

	path = [name '/' name];
	mkdir(name);
	dlmwrite([path '_g.dat'],G);
	dlmwrite([path '_s.dat'],system_info);
	dlmwrite([path '_c.dat'],c);
	dlmwrite([path '_p.dat'],p);
	
(2) compute the roundoff error:

- read the file product by the above script with `[F,I,J ,G ,n,d,k] = read_examples(program_name,"MAX");`
- execute `[ bound,build_time,solving_time ] = solve_examples( F,G,I,J,d,k,tag);`, where `bound` will be the roundoff error proportionnal to `epsilon` the machine precision.


## Installation instructions
### Prerequisites
FPKriSten is implemented in `Matlab2015a`, thus it is not guaranteed it will work on a previous (or later) version.

Moreover, FPKriSten relies on two external libraries for Matlab:

- `Yalmip` for the poynomial parsing (Free)
- `Cplex` IBM LP solver (Free for Academics) (added for a short period of time in the repository for linux 64)


### Download
`FPKriSten` is maintained as a GitHub repository at the address https://github.com/roccaa/FPKriSten.git

It can be obtained either by typing the shell command:

$ git clone https://github.com/roccaa/FPKriSten.git

or by downloading the ZIP archive at https://github.com/roccaa/FPKriSten.git

### Benchmarks

In the matlab command window, add to path the `Cplex`( in the  `matlab/x86-64_linux` directory if you are with linux 64 ) and `Yalmip` librairies.

Add to path the directories `Round_error` and `SBSOS` (dot not `SBSOS` from another source as the code files as been modified in this present version):

	$ addpath(genpath(Round_error/));addpath(genpath(SBSOS/));
	
Go to `Round_error` directory and build the models:

	$ cd Round_error/;
	$ build_all;

Return to `Round_error/` as the current directory has changed, and execute the `do_all` script:

	$ cd ..
	$ do_all:
	
Show the results summary:	

	$ result


