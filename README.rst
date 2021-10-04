Thor
====

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for Charged Particle Beam Dynamics
------------------------------------------------------------------------

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal as a *beam dynamics library*,
by the author, 1990; with care taken for the software architecture & resulting records/modules
– akin to *objects* although not explicitly supported by the grammar –  to reflect the structure of the mathematical objects describing
the underlying *beam dynamics model*.

The resulting C code, see below, has now been re-factored by introducing a C++ "*beam line class*";
i.e., to recover the transparency & simplicity of the original model & code.

Nota Bene: Although the beam dynamics model had to be replaced & the model/code re-architectured/structured as a Pascal beam dynamics libary, the code was named Tracy-2, i.e., after a precursor Tracy::

H\. Nishimura "*TRACY, A Tool for Accelerator Design and Analysis*" `EPAC 1988.`_

.. _`EPAC 1988.`: https://accelconf.web.cern.ch/e88/PDF/EPAC1988_0803.PDF

i.e., based on the linearized quadratic Hamiltonian for linear optics design. However, it was prototyped/implemented by extending the *standard procedures & functions* for the *Pascal-S compiler/interpreter* by N. Wirth::

N. Wirth *PASCAL-S A Subset and its Implementation* `Institut für Informatik (1975).`_

.. _`Institut für Informatik (1975).`: http://pascal.hansotten.com/uploads/pascals/PASCAL-S%20A%20subset%20and%20its%20Implementation%20012.pdf

which we liked we named it Tracy-2.

Thor

Contributions
-------------
* The symplectic integrator for RADIA kick maps was implemented by Laurent Nadolski, SOLEIL, 2002.

* The original Pascal library/code was machine translated to C by Michael Boege, SLS::

  M\. Böge, J. Chrin *A CORBA Based Client-Server Model for Beam Dynamics Applications* `ICALEPS 1999.`_

  .. _`ICALEPS 1999.`: https://accelconf.web.cern.ch/ica99/papers/mc1p61.pdf

  with `p2c.`_

  .. _`p2c.`: http://users.fred.net/tds/lab/p2c/historic/daves.index-2012Jul25-20-44-55.html

* Similarly, James Rowland utilized the C version to implement a *virtual accelerator* interfaced to EPICS as a *Virtual Input Output Controller* (VIOC)::

  M\. Heron, J. Rowland, et al *Progress on the Implementation of the DIAMOND Control System* `ICALEPCS 2005.`_

  .. _`ICALEPCS 2005.` https://accelconf.web.cern.ch/ica05/proceed-ings/pdf/P1_018.pdf

* Python interface::

  Intial demo/prototype & guidelines by Jan Chrin, PSI, 2017.
  
  J\. Chrin *Channel Access from Cython (and other Cython use cases)* `EPICS Collaboration Meeting 2017.`_
  
  .. _`EPICS Collaboration Meeting 2017.`: https://indico.esss.lu.se/event/889/contributions/7038/attachments/6800/9762/Cython_EpicsTM_Oct2017_Barcelona.pdf

  Guidelines & automated regression testing bootstrapped by Pierre Schnizer.


Requirements
------------
* (GNU compatible) C/C++ compiler
* GNU autoconf/automake environment and libtool.
* GNU Scientific Library (GSL): https://www.gnu.org/software/gsl.
* Armadillo (for linear algebra): http://arma.sourceforge.net.
* Python https://www.python.org/ for the python interface

The library uses the range checking inmplementation of e.g. `std::vector` as
provided by GNU C++; thus its dependency on the GUN compiler collections

To install
----------

Setup of repository
~~~~~~~~~~~~~~~~~~~

Dowload the repository and checkout the proper branch. Here it's assumed you
will use the directoy `git_repos/tracy-3.5` in your home directory for the
tracy code tree.

For this use the following commands to create the directoy `git_repos`
and to clone the tree into the tracy-3.5 directory.

.. code:: shell

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/tracy-3.5.git
   cd tracy-3.5

Then select the proper tree by

.. code:: shell

   git checkout tracy-3.5_scsi



C++ library
~~~~~~~~~~~


First create environment variable $TRACY_LIB. This will be the prefix where the
built library and include files will be installed later on e.g:

.. code:: shell

   export TRACY_LIB=$HOME/git_repos/tracy-3.5


To build the library use:

.. code:: shell

   cd tracy-3.5
   libtoolize
   ./bootstrap
   ./configure --prefix=$TRACY_LIB
   make
   make install

Please note: using the dynamic library in non standard location will require
proper set up of the environment later on (e.g. adding the directory where the
library is located to `LD_LIBRARY_PATH` environment variable).


Python interface
~~~~~~~~~~~~~~~~

The python interface is based on https://github.com/pybind/pybind11. Building this interface
requires to select the proper directory

.. code:: shell

  cd git_repos
  cd tracy-3.5/python

Install proper dependencies

.. code:: shell

    pip3 install -r requirements.txt


And build the extension e.g.

.. code:: shell

    python3 setup.py build
    python3 setup.py install

For further details of the build system see https://pypi.org/project/setuptools/


To run the regression tests
---------------------------

All regression tests can be run using

.. code:: shell

    pip3 install nose
    python3 setup.py nosetests

To run the demo/test program
----------------------------


.. code:: shell

    python3 examples/tst.py
