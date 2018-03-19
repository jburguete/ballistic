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
[!img](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D%2C%5Cquadd%5Cvec%7Bg%7D%3D%280%2C%5C%3B0%2C%5C%3B-g%29)

with g the gravity constant. The analytical solution is:
[!img](http://latex.codecogs.com/svg.latex?%5Cdot%7B%5Cvec%7Br%7D%7D%3D%5Cdot%7B%5Cvec%7Br%7D%7D_0%2B%5Cveg%7Bg%7D%5C%2Ct%5C%5C%0D%0A%5Cvec%7Br%7D%3D%5Cvec%7Br%7D_0%2B%5Cdot%7B%5Cvec%7Br%7D%7D_0%5C%2Ct%2B%5Cvec%7Bg%7D%5C%2Ct%5E2.)

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
