Diffusion1D
===========

Computes diffusion equation in 1D using 2 different methods

Diffusion Equation
------------------
The diffusion equation takes the form

∂u/∂t = D∂(∂u/∂x)/∂x

where we assume that D is spatially independent. In order to solve this computationally, we discretize the variable u in the [typical ways](https://me.ucsb.edu/~moehlis/APC591/tutorials/tutorial5/node3.html). We then solve the system using two difference methods.

Method 1
--------
Standard Crank-Nicolson scheme. Due to it being unconditionally stable, the algorithm is generally the go-to method for solving parabolic PDEs such as the diffusion equation. As such, we use it as a baseline for testing alternative methods.

Method 2
--------
Super Time-Stepping scheme of Alexiades 1996. This is an explicit scheme that uses [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials) to relax the time-step so that fewer time steps are required as compared to following `D*dt/dx/dx <= 1` rule typical of explicit methods. 

For more details, see the original paper proposing this method: http://www4.ncsu.edu/eos/users/g/gremaud/WWW/sup.pdf

Test Case
=========
This program is primarily used as a test case to see how accurate and fast the two methods compare to each other (others can be added easily-enough). Given a sine-wave for the initial condition, 

`u(x,0) = sin(pi*x)`

the solution at some time `t` later is

`u(x,t) = sin(pi*x) * exp(-pi*pi*t)`

The two methods iterate a discretized form of the first equation and evolve it in time from `t=0` to `t=0.1` seconds.

Output
======
The output file follows the [curve](http://www.visitusers.org/index.php?title=Reading_curve_data) format, which is easily readable in visit.

Currently, only the initial function and the final function are output to file. Users could modify it slightly to output at regular intervals, but I'm not interested in that, only the initial and final states.


Compilation
===========

This is compiled by entering `make` in the command line. It does require a Fortran compiler. And probably a terminal. And VisIt so you can view the solutions.


Disclaimer
----------
I do not take responsibility for any damages that occur to your computer through using or abusing my code.
