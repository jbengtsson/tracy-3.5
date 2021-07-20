tracy-3.5
=========

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics
------------------------------------------------------------------------

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal, by the author,
with care taken for the software architecture and resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying beam dynamics model.


Contributions
-------------
The symplectic integrator for RADIA kick maps was implemented by Laurent
Nadolski, SOLEIL, 2002.

The original Pascal library/code was machine translated to C (with p2c) by
Michael Boege, SLS, 1998.


Requirements
------------

   GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
   GNU autoconf/automake environment and libtool.
   "Numerical Recipes in C": http://www.nr.com.

To install
----------

First create environment variable $TRACY_LINK e.g.:

   export TRACY_LIB=$HOME/git_repos/tracy-3.5

then:

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/tracy-3.5.git
   cd tracy-3.5
   ./make_tracy.sh
