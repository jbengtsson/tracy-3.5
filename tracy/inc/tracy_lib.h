/* Tracy-3

   J. Bengtsson, BNL 2007

   NO   1   link to the linear TPSA (nv_tps = 1)
       >1   link to Berz' TPSA

*/

#ifndef TRACY_LIB_H
#define TRACY_LIB_H

#define _GLIBCXX_DEBUG          1
#define _GLIBCXX_DEBUG_PEDANTIC 1

// C standard library.
#include <stdio.h>
#include <stddef.h>
#include <setjmp.h>
#include <time.h>
#include <memory.h>

// C++ standard library.
#include <cstdlib>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>

#include <armadillo>
#include <gsl/gsl_linalg.h>

#include <vector>


using namespace std;


// Tracy-3
#include "field.h"

#if NO_TPSA == 1
  // linear TPSA
  #include "tpsa_lin.h"
  #include "tpsa_lin_pm.h"
#else
  // interface to M. Berz' TPSA
  #include "tpsa_for.h"
  #include "tpsa_for_pm.h"
#endif

#include "ety.h"
#include "eigenv.h"

#include "tracy.h"

#include "t2elem.h"
#include "t2cell.h"
#include "t2ring.h"
// #include "sigma_track.h"

#include "t2lat.h"
#include "rdmfile.h"
#include "prtmfile.h"

#include "set_errors.h"
#include "lsoc.h"

// #include "orb_corr.h"
// #include "param.h"
// #include "dynap.h"

// #include "physlib.h"
// #include "fft.h"

#include "radia2tracy.h"
// #include "modnaff.h"
// #include "naffutils.h"
// #include "complexeheader_naff.h"
// #include "soleillib.h"

// #include "nsls-ii_lib.h"


// Truncated Power Series Algebra (TPSA)
extern const int nv_tps, nd_tps, iref_tps;
extern int       no_tps, ndpt_tps;
extern double    eps_tps;

extern LatticeType lat;

#endif
