BALLISTIC (0.1.0 version)
=========================

A software to benchmark ballistic models.

AUTHORS
-------

* Javier Burguete Tolosa (jburguete@eead.csic.es)

INTRODUCTION
------------

This software benchmark some ballistic models comparing numerical versus
analytical solutions.

Four 

PARABOLIC TRAJECTORY
_____________________

The program solves the movement equation:

[!equation 1](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D%2C%5C%5C%0D%0A%5Cvec%7Bg%7D%3D%280%2C%5C%3B0%2C%5C%3B-g%29%2C)

with g the gravity constant. The analytical solution is:

[!equation 2](http://latex.codecogs.com/svg.latex?%5Cdot%7B%5Cvec%7Br%7D%7D%3D%5Cdot%7B%5Cvec%7Br%7D%7D_0%2B%5Cvec%7Bg%7D%5C%2Ct%2C%5C%5C%0D%0A%5Cvec%7Br%7D%3D%5Cvec%7Br%7D_0%2B%5Cdot%7B%5Cvec%7Br%7D%7D_0%5C%2Ct%2B%5Cfrac12%5C%2C%5Cvec%7Bg%7D%5C%2Ct%5E2.)

1ST RESISTANCE MODEL
____________________

In this case the movement equation is:

[!equation 3](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D-%5Clambda%5C%2C%5Cleft%28%5Cdot%7B%5Cvec%7Br%7D%7D-%5Cvec%7Bw%7D%5Cright%29%2C%5C%5C%0D%0A%5Cvec%7Bw%7D%3D%5Cleft%28w_x%2C%5C%3Bw_y%5C%3B0%5Cright%29%2C)

with \lambda a resistance coefficient and w_x and w_y the wind velocity vector
components. The analytical solution is:
[!equation 4](http://latex.codecogs.com/svg.latex?http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D-%5Clambda%5C%2C%5Cleft%28%5Cdot%7B%5C)

FILES
-----

* README.md: Readme file.
* LICENSE: BSD type license.
* TODO: List of tasks TO DO.
* \*.h: Header files.
* \*.c: Source files.
* Doxyfile: configuration file to generate doxygen documentation.
* Makefile: to build the tests.

BUILDING THIS LIBRARY ON OTHER PROGRAMS
---------------------------------------

REQUIRED LIBRARIES AND UTILITIES
________________________________

Mandatory:
* [gcc](https://gcc.gnu.org) or [clang](http://clang.llvm.org) (to compile the
  source code).
* [make](http://www.gnu.org/software/make) (to build the executable file).
* [pkg-config](http://www.freedesktop.org/wiki/Software/pkg-config) (to find the
  libraries to compile).
* [glib](https://developer.gnome.org/glib) (extended utilities of C to work with
  data, lists, mapped files, regular expressions, using multicores in shared
  memory machines, ...).
* [gettext](http://www.gnu.org/software/gettext) (to work with different
  locales and languages).

Optional to build documentation:
* [doxygen](http://www.stack.nl/~dimitri/doxygen) (standard comments format to
  generate documentation).
* [latex](https://www.latex-project.org/) (to build the PDF manuals).

BUILDING INSTRUCTIONS
______________________

1. Load the last library version:
> $ git clone https://github.com/jburguete/ballistic

2. Access to the program directory:
> $ cd ballistic/0.1.0

3. Build:
> $ make

MAKING REFERENCE MANUAL INSTRUCTIONS (latex/refman.pdf file)
------------------------------------------------------------

Execute on a terminal:
> $ cd ballistic/0.1.0
>
> $ doxygen
>
> $ cd latex
>
> $ make
