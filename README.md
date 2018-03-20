BALLISTIC (0.1.2 version)
=========================

A software to benchmark ballistic models.

AUTHORS
-------

* Javier Burguete Tolosa (jburguete@eead.csic.es)

MOVEMENT EQUATIONS
------------------

This software benchmark some ballistic models comparing numerical versus
analytical solutions.

Four equations with analytical solution are implemented.

PARABOLIC TRAJECTORY
_____________________

The program solves the movement equation:

![equation 1](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D%2C%5C%5C%0D%0A%5Cvec%7Bg%7D%3D%280%2C%5C%3B0%2C%5C%3B-g%29%2C)

with g the gravity constant. The analytical solution is:

![equation 2](http://latex.codecogs.com/svg.latex?%5Cdot%7B%5Cvec%7Br%7D%7D%3D%5Cdot%7B%5Cvec%7Br%7D%7D_0%2B%5Cvec%7Bg%7D%5C%2Ct%2C%5C%5C%0D%0A%5Cvec%7Br%7D%3D%5Cvec%7Br%7D_0%2B%5Cdot%7B%5Cvec%7Br%7D%7D_0%5C%2Ct%2B%5Cfrac12%5C%2C%5Cvec%7Bg%7D%5C%2Ct%5E2.)

1ST RESISTANCE MODEL
____________________

In this case the movement equation is:

![equation 3](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D-%5Clambda%5C%2C%5Cleft%28%5Cdot%7B%5Cvec%7Br%7D%7D-%5Cvec%7Bw%7D%5Cright%29%2C%5C%5C%0D%0A%5Cvec%7Bw%7D%3D%5Cleft%28w_x%2C%5C%3Bw_y%5C%3B0%5Cright%29%2C)

with λ a resistance coefficient and w<sub>x</sub> and w<sub>y</sub> the wind
velocity vector components. The analytical solution is:

![equation 4](http://latex.codecogs.com/svg.latex?%5Cdot%7B%5Cvec%7Br%7D%7D%3D%5Cdot%7B%5Cvec%7Br%7D%7D_0%5C%2C%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29%2B%5Cleft%28%5Cvec%7Bw%7D%2B%5Cfrac%7B%5Cvec%7Bg%7D%7D%7B%5Clambda%7D%5Cright%29%5C%2C%5Cleft%5B1-%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29%5Cright%5D%2C%5C%5C%0D%0A%5Cvec%7Br%7D%3D%5Cvec%7Br%7D_0%2B%5Cleft%28%5Cvec%7Bw%7D%2B%5Cfrac%7B%5Cvec%7Bg%7D%7D%7B%5Clambda%7D%5Cright%29%5C%2Ct%2B%5Cfrac%7B%5Cdot%7B%5Cvec%7Br%7D%7D_0-%5Cvec%7Bw%7D-%5Cvec%7Bg%7D%2F%5Clambda%7D%7B%5Clambda%7D%5C%2C%5Cleft%5B1-%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29%5Cright%5D.%0D%0A)

2ND RESISTANCE MODEL
____________________

The movement equations are:

![equation 5](http://latex.codecogs.com/svg.latex?%5Cddot%7Bx%7D%3D-%5Clambda%5C%2C%5Cleft%7C%5Cdot%7Bx%7D-w_x%5Cright%7C%5C%2C%5Cleft%28%5Cdot%7Bx%7D-w_x%5Cright%29%2C%5C%5C%0D%0A%5Cddot%7By%7D%3D-%5Clambda%5C%2C%5Cleft%7C%5Cdot%7By%7D-w_y%5Cright%7C%5C%2C%5Cleft%28%5Cdot%7By%7D-w_y%5Cright%29%2C%5C%5C%0D%0A%5Cddot%7Bz%7D%3D-g-%5Clambda%5C%2C%5Cleft%7C%5Cdot%7Bz%7D%5Cright%7C%5C%2C%5Cdot%7Bz%7D%2C)

The analytical solution is:

![equation 6](https://latex.codecogs.com/gif.latex?%5Cdot%7Bx%7D%3Dw_x&plus;%5Cfrac%7B%5Cdot%7Bx%7D_0-w_x%7D%7B1&plus;%5Clambda%5C%2C%5Cleft%7C%5Cdot%7Bx%7D_0-w_x%5Cright%7C%5C%2Ct%7D%2C%5C%5C%20%5Cdot%7By%7D%3Dw_y&plus;%5Cfrac%7B%5Cdot%7By%7D_0-w_y%7D%7B1&plus;%5Clambda%5C%2C%5Cleft%7C%5Cdot%7By%7D_0-w_y%5Cright%7C%5C%2Ct%7D%2C%5C%5C%20%5Cdot%7Bz%7D_0%5Cleq%200%5CRightarrow%5Cquad%5Cdot%7Bz%7D%3D%5Csqrt%7B%5Cfrac%7Bg%7D%7B%5Clambda%7D%7D%5C%2C%20%5Cfrac%7B%5Cdot%7Bz%7D_0%5C%2C%5Ccosh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%20-%5Csqrt%7Bg/%5Clambda%7D%5C%2C%5Csinh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%7D%20%7B%5Csqrt%7Bg/%5Clambda%7D%5C%2C%5Ccosh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%20-%5Cdot%7Bz%7D_0%5C%2C%5Csinh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%7D%2C%5C%5C%20%5Cdot%7Bz%7D_0%3E0%2C%5C%3B%20t%5Cleq%5Cfrac%7B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%7D%20%7B%5Csqrt%7Bg%5C%2C%5Clambda%7D%7D%5CRightarrow%5Cquad%20%5Cdot%7Bz%7D%3D%5Csqrt%7B%5Cfrac%7Bg%7D%7B%5Clambda%7D%7D%20%5C%2C%5Ctan%5Cleft%5B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Cfrac%7B%5Clambda%7D%7Bg%7D%7D%5Cright%29%20-%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%5D%2C%5C%5C%20%5Cdot%7Bz%7D_0%3E0%2C%5C%3B%20t%3E%5Cfrac%7B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%7D%20%7B%5Csqrt%7Bg%5C%2C%5Clambda%7D%7D%5CRightarrow%5Cquad%20%5Cdot%7Bz%7D%3D-%5Csqrt%7B%5Cfrac%7Bg%7D%7B%5Clambda%7D%7D%5C%2C%5Ctanh%5Cleft%5B%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%20-%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Cfrac%7B%5Clambda%7D%7Bg%7D%7D%5Cright%29%5Cright%5D.)

![equation 6b](http://latex.codecogs.com/svg.latex?x%3Dx_0%2Bw_x%5C%2Ct%2B%5Cfrac%7B%5Cdot%7Bx%7D_0-w_x%7D%7B%5Clambda%5C%2C%5Cleft%7C%5Cdot%7Bx%7D_0-w_x%5Cright%7C%7D%0D%0A%5C%2C%5Cln%5Cleft%281%2B%5Clambda%5C%2C%5Cleft%7C%5Cdot%7Bx%7D_0-w_x%5Cright%7C%5C%2Ct%5Cright%29%2C%5C%5C%0D%0Ay%3Dy_0%2Bw_y%5C%2Ct%2B%5Cfrac%7B%5Cdot%7By%7D_0-w_y%7D%7B%5Clambda%5C%2C%5Cleft%7C%5Cdot%7By%7D_0-w_y%5Cright%7C%7D%0D%0A%5C%2C%5Cln%5Cleft%281%2B%5Clambda%5C%2C%5Cleft%7C%5Cdot%7By%7D_0-w_y%5Cright%7C%5C%2Ct%5Cright%29%2C%0D%0A)

![equation 6c](http://latex.codecogs.com/svg.latex?https://latex.codecogs.com/gif.latex?%5Cdot%7Bz%7D_0%5Cleq%200%5CRightarrow%5Cquad%20z%3Dz_0-%5Cfrac%7B%5Cln%5Cleft%5B%5Ccosh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%20-%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5C%2C%20%5Csinh%5Cleft%28%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%29%5Cright%5D%7D%7B%5Clambda%7D%2C%5C%5C%20%5Cdot%7Bz%7D_0%3E0%2C%5C%3B%20t%5Cleq%5Cfrac%7B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%7D%20%7B%5Csqrt%7Bg%5C%2C%5Clambda%7D%7D%5CRightarrow%5Cquad%20z%3Dz_0&plus;%5Cfrac%7B1%7D%7B%5Clambda%7D%5C%2C%5Cln%5Cleft%5C%7B%5Cfrac%7B%5Ccos%5Cleft%5B%20%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29-%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%5Cright%5D%7D%20%7B%5Ccos%5Cleft%5B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%5Cright%5D%7D%5Cright%5C%7D%2C%5C%5C%20%5Cdot%7Bz%7D_0%3E0%2C%5C%3B%20t%3E%5Cfrac%7B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%7D%20%7B%5Csqrt%7Bg%5C%2C%5Clambda%7D%7D%5CRightarrow%5Cquad%20z%3Dz_0-%5Cfrac%7B%5Cln%5Cleft%5C%7B%20%5Ccos%5Cleft%5B%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%5Cright%5D%20%5C%2C%5Ccosh%5Cleft%5B%5Csqrt%7Bg%5C%2C%5Clambda%7D%5C%2Ct%20-%5Carctan%5Cleft%28%5Cdot%7Bz%7D_0%5C%2C%5Csqrt%7B%5Clambda/g%7D%5Cright%29%5Cright%5D%5Cright%5C%7D%7D%7B%5Clambda%7D.)

FORCED MODEL
____________

In this case the movement equation is:

![equation 7](http://latex.codecogs.com/svg.latex?%5Cddot%7B%5Cvec%7Br%7D%7D%3D%5Cvec%7Bg%7D%2B%5Cvec%7Bw%7D%5C%2C%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29)

The analytical solution is:

![equation 8](http://latex.codecogs.com/svg.latex?%5Cdot%7B%5Cvec%7Br%7D%7D%3D%5Cdot%7B%5Cvec%7Br%7D%7D_0%2B%5Cvec%7Bg%7D%5C%2Ct%2B%5Cfrac%7B%5Cvec%7Bw%7D%7D%7B%5Clambda%7D%5Cleft%5B1-%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29%5Cright%5D%2C%5C%5C%0D%0A%5Cvec%7Br%7D%3D%5Cvec%7Br%7D_0%2B%5Cleft%28%5Cdot%7B%5Cvec%7Br%7D%7D_0%2B%5Cfrac%7B%5Cvec%7Bw%7D%7D%7B%5Clambda%7D%5Cright%29%5C%2Ct%2B%5Cfrac12%5C%2C%5Cvec%7Bg%7D%5C%2Ct%5E2%2B%5Cfrac%7B%5Cvec%7Bw%7D%7D%7B%5Clambda%5E2%7D%5C%2C%5Cleft%5B1-%5Cexp%5Cleft%28-%5Clambda%5C%2Ct%5Cright%29%5Cright%5D.)

FILES
-----

* README.md: Readme file.
* LICENSE: BSD type license.
* TODO: List of tasks TO DO.
* 0.1.2/\*.h: Header files.
* 0.1.2/\*.c: Source files.
* 0.1.2/Doxyfile: configuration file to generate doxygen documentation.
* 0.1.2/Makefile: to build the executable program.
* 0.1.2/case\*: different example cases input files.
* 0.1.2/script: a shell script to generate some convergence plots.
* 0.1.2/plot: a GNUPlot input file to do the plots.

BUILDING THE EXECUTABLE FILE
----------------------------

REQUIRED LIBRARIES AND UTILITIES
________________________________

Mandatory:
* [gcc](https://gcc.gnu.org) or [clang](http://clang.llvm.org): to compile the
  source code.
* [make](http://www.gnu.org/software/make): to build the executable file.
* [pkg-config](http://www.freedesktop.org/wiki/Software/pkg-config): to find the
  libraries to compile.
* [gsl](http://www.gnu.org/software/gsl): to generate random numbers.

Optional:
* [doxygen](http://www.stack.nl/~dimitri/doxygen): to generate documentation.
* [latex](https://www.latex-project.org): to build the PDF manual.
* [gnuplot](http://gnuplot.info): to do the graphical plots.

BUILDING INSTRUCTIONS
______________________

1. Load the last program version:
> $ git clone https://github.com/jburguete/ballistic

2. Access to the program directory:
> $ cd ballistic/0.1.2

3. Build:
> $ make

MAKING REFERENCE MANUAL INSTRUCTIONS (latex/refman.pdf file)
------------------------------------------------------------

Execute on a terminal:
> $ cd ballistic/0.1.2
>
> $ doxygen
>
> $ cd latex
>
> $ make
