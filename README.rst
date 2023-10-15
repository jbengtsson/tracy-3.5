Tracy-3.5
=======

The symplectic integrator for realistic modeling of magnetic lattices for ring-based synchrotrons was initially implemented as a *Pascal module/beam dynamics software library*, by the author 1990, as an *on-line model* to guide the ALS commissioning. In particular, care was taken for the software architecture & resulting records/modules – akin to *objects* although not explicitly supported by the artificial language grammar – to reflect the *structure of the mathematical objects* describing the underlying *beam dynamics model*.

Hence, the code was also benchmarked & calibrated as part of the ALS commissioning:

  J\. Bengtsson, M. Meddahi *Modeling of Beam Dynamics and Comparison with Measurements for the Advanced Light Source (ALS)* `EPAC 1994.`_

  .. _`EPAC 1994.`: https://accelconf.web.cern.ch/e94/PDF/EPAC1994_1021.PDF

Malika Meddahi was a postdoc at ALS who did her thesis at CERN (for which she received the *Prix Daniel Guinier* in France that year):

  *Effets faisceau-faisceau dans le collisionneur protons-antiprotons du SPS* `CERN-SL-91-30-BI (1991)`_

  .. _`CERN-SL-91-30-BI (1991)`: http://cds.cern.ch/record/223301)

The resulting C code, see below, has now been re-factored by introducing a C++ *beam line class*; i.e., to recover the transparency & simplicity of the original *beam dynamics model*.

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





tracy-3.5

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics.

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal, by the author,
with care taken for the software architecture and resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying beam dynamics model.


The symplectic integrator for RADIA kick maps was implemented by Laurent
Nadolski, SOLEIL, 2002.

The original Pascal library/code was machine translated to C (with p2c) by
Michael Boege, SLS, 1998.


Requirements:

   GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
   GNU autoconf/automake environment and libtool.
   "Numerical Recipes in C": http://numerical.recipes.

To install:

First create environment variable $TRACY_LINK e.g.:

   export TRACY_LIB=$HOME/git_repos/tracy-3.5

then:

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/tracy-3.5.git
   cd tracy-3.5
   ./make_tracy.sh
