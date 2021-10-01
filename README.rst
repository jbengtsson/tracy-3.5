Thor
====

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics
------------------------------------------------------------------------

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal as a "*beam dynamics library*",
by the author, 1990; with care taken for the software architecture & resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying "*beam dynamics model*".
The resulting C code, see below, has now been re-factored by introducing a C++ "*beam line class*";
i.e., to recover the transparency & simplicity of the original model & code.


Contributions
-------------
* The symplectic integrator for RADIA kick maps was implemented by Laurent Nadolski, SOLEIL, 2002.

* The original Pascal library/code was machine translated to C by Michael Boege, SLS, 1998 (with p2c).

* Python interface::

  Intial demo/prototype & guidelines by Jan Chrin, PSI, 2017.
  
  Guidelines & automated regression testing prototyped by Pierre Schnizer.


Requirements
------------
* GNU C/C++ compiler.
* GNU Python.
* GNU Scientific Library (GSL): https://www.gnu.org/software/gsl.
* Armadillo (for linear algebra): http://arma.sourceforge.net.


To install
----------
.. code:: shell

  mkdir git_repos
  cd git_repos
  git clone git@github.com:jbengtsson/tracy-3.5_scsi.git
  cd tracy-3.5_scsi/python
  python3 setup.py build

To run the regression tests
---------------------------
.. code:: shell

python3 setup.py nosetests

To run the demo/test program
------------------------
.. code:: shell

python3 examples/tst.py


Obsolete
--------

* GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
* GNU autoconf/automake environment and libtool.
* GNU scientific library see https://www.gnu.org/software/gsl/

First create environment variable $TRACY_LINK e.g.:

.. code:: shell

   export TRACY_LIB=$HOME/git_repos/tracy-3.5

then:

.. code:: shell

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/tracy-3.5.git
   cd tracy-3.5
   ./make_tracy.sh
