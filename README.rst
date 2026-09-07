Tracy-3.5
=========

The symplectic integrator for realistic modeling of magnetic lattices for ring-based synchrotrons was initially implemented as a *Pascal module/beam dynamics software library*, by the author 1990, as an *on-line model* to guide the ALS commissioning. In particular, care was taken for the software architecture & resulting records/modules – akin to *objects* although not explicitly supported by the artificial language grammar – to reflect the *structure of the mathematical objects* describing the underlying *beam dynamics model*.

Hence, the code was also benchmarked & calibrated as part of the ALS commissioning:

  J\. Bengtsson, M. Meddahi *Modeling of Beam Dynamics and Comparison with Measurements for the Advanced Light Source (ALS)* `EPAC 1994.`_

  .. _`EPAC 1994.`: https://accelconf.web.cern.ch/e94/PDF/EPAC1994_1021.PDF


Remark: Although the entire *beam dynamics model* had to be replaced & the model/code/"approach" re-architectured & structured – for a reusable approach – as a *Pascal beam dynamics libary* (standard practise in software engineering), the code was named *Tracy-2*, i.e., inspired by the, somewhat archaic demo/prototype/concept *Tracy*:

  H\. Nishimura *TRACY, A Tool for Accelerator Design and Analysis* `EPAC 1988`_

  .. _`EPAC 1988`: https://accelconf.web.cern.ch/e88/PDF/EPAC1988_0803.PDF

for which the *beam dynamics model* was based on the *linearized quadratic Hamiltonian*:

  .. image:: images/H_2.png

for *linear optics design*. I.e., for a *bare lattice* with *mid-plane symmetry*.

E.g. by not having figured out/mastered how to pass records (structures in C) as function/procedure variables – vs. scalars only – for the Pascal-S compiler/interpreter to the beam dynamics library. The API was rather poor/sloppy. I.e., not scalable and thus ill suited to cope with the complexity of a realistic model. As expressed by Forest in the title of:

  E\. Forest *A Hamiltonian-Free Description of Single Particle Dynamics for Hopelessly Complex Periodic Systems* `J. Math. Phys. 31 (1990).`_

  .. _`J. Math. Phys. 31 (1990).`: http://dx.doi.org/10.1063/1.528795%7D

Hence, the one thing we did find useful for a realistic on-line model – having already implemented an on-line model as a sci fellow for LEAR, CERN, in the late 1980s and before that having worked as a teaching assistent at the *dept. of Software Engineering, Lund Inst. of Tech, Sweden* (next to *MAX Lab*) while pursuing a MsSci EE – and adopted for ALS. Was the implementation of the beam dynamics model as an *extension of the standard procedures & functions* for the *Pascal-S compiler/interpreter* by N. Wirth (implemented/coded in it's native grammar); architected as a Pascal software library/module:

  M\. Rees, D\. Robson *Practical Compiling with Pascal-S* `(Addison-Wesley, 1988).`_

  .. _`(Addison-Wesley, 1988).`: https://books.google.com/books?id=hLomAAAAMAAJ

  S\. Pemberton, M\. Daniels *The P4 Compiler and Interpreter* `(1982).`_

  .. _`(1982).`: https://homepages.cwi.nl/~steven/pascal/book/pascalimplementation.html

  N\. Wirth *PASCAL-S: A Subset and its Implementation* `Institut für Informatik, ETH, Zürich (1975).`_

  .. _`Institut für Informatik, ETH, Zürich (1975).`: http://pascal.hansotten.com/uploads/pascals/PASCAL-S%20A%20subset%20and%20its%20Implementation%20012.pdf

  *Pascal-P6* https://sourceforge.net/projects/pascal-p6.


Contributions
-------------
* The symplectic integrator for *RADIA kick maps*:

    P\. Elleaume *A New Approach to the Electron Beam Dynamics in Undulators and Wigglers* `EPAC 1992.`_

    .. _`EPAC 1992.`: https://accelconf.web.cern.ch/e92/PDF/EPAC1992_0661.PDF

  was implemented by Laurent Nadolski, SOLEIL, 2002.

* The original *Pascal library/code* was machine translated to C and re-used to implement a *model server* for the SLS commissioning:

    M\. Böge *Update on TRACY-2 Documentation* `SLS Tech Note SLS-TME-TA-1999-0002 (1999).`_

    .. _`SLS Tech Note SLS-TME-TA-1999-0002 (1999).`: http://ados.web.psi.ch/slsnotes/tmeta9902.pdf

    M\. Böge, J. Chrin *A CORBA Based Client-Server Model for Beam Dynamics Applications* `ICALEPCS 1999.`_

    .. _`ICALEPCS 1999.`: https://accelconf.web.cern.ch/ica99/papers/mc1p61.pdf

  with `p2c.`_

    .. _`p2c.`: http://users.fred.net/tds/lab/p2c/historic/daves.index-2012Jul25-20-44-55.html

* Similarly, James Rowland re-used the C version to implement a *Virtual Accelerator* interfaced to EPICS as a *Virtual Input Output Controller* (VIOC):

    M\. Heron, J. Rowland, et al *Progress on the Implementation of the DIAMOND Control System* `ICALEPCS 2005.`_

    .. _`ICALEPCS 2005.`: https://accelconf.web.cern.ch/ica05/proceed-ings/pdf/P1_018.pdf

* Besides, a subset of the internal *numerical engine* was manually translated to C and re-used for:

    A\. Terebilo *Accelerator Toolbox for MATLAB* `SLAC-PUB-8732 (2001).`_

    .. _`SLAC-PUB-8732 (2001).`: http://www-public.slac.stanford.edu/sciDoc/docMeta.aspx?slacPubNumber=SLAC-PUB-8732


However, unfortunately, the effort required to understand how the linear optics is obtained from the symplectic integrater, by utilising Power Series Algebra (TPSA),
for a self-consistent approach was not pursued – a major shortcut – i.e., the resulting model/code is not self-consistent.


Requirements
------------

  * GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
  * GNU autoconf/automake environment and libtool.
  * *GNU Scientific Library* (GSL): https://www.gnu.org/software/gsl.
    Located automatically via ``gsl-config``; see *Troubleshooting* below.
  * *Numerical Recipes in C*: http://numerical.recipes.
  * *Thor*, for the scripts under ``projects/src/ptc`` only.

E.g. on macOS::

   brew install autoconf automake libtool gcc gsl

or on Debian/Ubuntu::

   sudo apt install build-essential gfortran autoconf automake libtool libgsl-dev

Installation
------------

The build has two stages: the Tracy libraries at the top level, then the
scripts under ``projects/``. Both are needed to run anything.

1. Clone the repository::

      git clone git@github.com:jbengtsson/tracy-3.5.git
      cd tracy-3.5

2. Set the three paths the build needs. ``TRACY_LIB`` is this tree;
   ``NUM_REC`` is your Numerical Recipes tree; ``THOR_LIB`` is your Thor
   tree::

      export TRACY_LIB=$PWD
      export NUM_REC=$HOME/Software/numrec
      export THOR_LIB=$HOME/Software/thor-2.0

   If any is unset or wrong, ``configure`` stops with a message naming it,
   rather than failing later with a confusing link error.

3. Build and install the libraries into ``./tracy/lib``::

      ./make_tracy-3.5.sh

4. Build the scripts::

      cd projects
      ./make_projects.sh

   The executables land in ``projects/src`` and ``projects/src/ptc``.

Both scripts pass any extra arguments through to ``configure``, and both
print a summary of the compilers, versions and library locations they
actually used.

Troubleshooting
---------------

**GSL is not found.** The build asks ``gsl-config`` – installed alongside
GSL – where its headers and libraries are. If it is not on your ``PATH``,
point at the install prefix directly, i.e. the directory holding ``bin``,
``include`` and ``lib``::

   ./make_tracy-3.5.sh --with-gsl=DIR

**Rebuilding after moving or renaming the tree.** The generated Makefiles
and libtool scripts hold absolute paths, so a moved tree will fail with
``Library not loaded`` or a ``config/missing`` error. The reliable reset is::

   git clean -xfd

That deletes every file in the tree that git does not track: ``-x`` includes
the ones ``.gitignore`` covers, ``-d`` recurses into directories, and ``-f``
is required before git will delete anything. Since every build artefact here
is gitignored, this removes the entire build and leaves your tracked sources
untouched. Run ``git clean -xfdn`` first – ``-n`` lists what would go without
deleting it – and note that any untracked file you wanted to keep, such as a
lattice or an output file you have not committed, will be deleted too.

**aclocal cannot find the m4 directory.** It needs the libtool macros in
``m4/``, which are generated by ``libtoolize`` (``glibtoolize`` on macOS).
Both build scripts do this for you; if you run ``./bootstrap`` by hand, run
``libtoolize`` first.
