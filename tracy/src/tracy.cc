
#include "tracy_lib.h"

#include "field.cc"

#if NO == 1
// linear TPSA
#include "tpsa_lin.cc"
#include "tpsa_lin_pm.cc"
#else
// interface to M. Berz' TPSA
#include "tpsa_for_pm.cc"
#endif

#include "mathlib.cc"

#include "ety.cc"
#include "eigenv.cc"

#include "LEGO.cc"

#include "t2lat.cc"
#include "t2elem.cc"
#include "t2cell.cc"
#include "t2ring.cc"

#include "pascalio.cc"

#include "lsoc.cc"

#include "prtmfile.cc"
#include "rdmfile.cc"

#include "fft.cc"

#include "physlib.cc"

#include "naffutils.cc"

#include "modnaff.cc"
#include "radia2tracy.cc"
#include "soleillib.cc"

#include "nsls-ii_lib.cc"
#include "orb_corr.cc"
#include "param.cc"
#include "dynap.cc"


// Truncated Power Series Algebra (TPSA)
const int     nv_tps   = ss_dim, // no of variables
  nd_tps   = 3,      // no of degrees of freedom
  iref_tps = 0;      /* file with resonances to be excluded from
			the map normal form: fort.7 */

double eps_tps  = 1e-25;  // floating point truncation


// Instantiate templates.

template class ss_vect<double>;
template class ss_vect<tps>;


template void p_rot(double, ss_vect<double> &);
template void p_rot(double, ss_vect<tps> &);


template void get_B2(const double, const double [], const ss_vect<double> &,
		     double &, double &);
template void get_B2(const double, const tps [], const ss_vect<tps> &,
		     tps &, tps &);

template void radiate(ss_vect<double> &, const double, const double,
		      const double []);
template void radiate(ss_vect<tps> &, const double, const double,
		      const tps []);

template void radiate_ID(ss_vect<double> &, const double, const double &);
template void radiate_ID(ss_vect<tps> &, const double, const tps &);

template void Drift(const double, ss_vect<double> &);
template void Drift(const double, ss_vect<tps> &);

template void bend_fringe(const double, ss_vect<double> &);
template void bend_fringe(const double, ss_vect<tps> &);

template void EdgeFocus(const double, const double, const double,
			ss_vect<double> &);
template void EdgeFocus(const double, const double, const double,
			ss_vect<tps> &);

template void quad_fringe(const double, ss_vect<double> &);
template void quad_fringe(const double, ss_vect<tps> &);


template void thin_kick(const int, const double [], const double, const double,
			const double, ss_vect<double> &);
template void thin_kick(const int, const double [], const double, const double,
			const double, ss_vect<tps> &);


template void SplineInterpolation2(double &, double &, double &, double &,
				   CellType &, bool &);
template void SplineInterpolation2(tps &, tps &, tps &, tps &,
				   CellType &, bool &);

template void spline(const double [], const double [], int const,
		     double const, const double, double []);
template void spline(const double [], const tps [], int const,
		     double const, const double, tps []);

template void splint(const double[], const double [], const double [],
		     const int, const double &, double &);
template void splint(const double[], const double [], const double [],
		     const int, const tps &, tps &);
template void splint(const double[], const tps [], const tps [],
		     const int, const tps &, tps &);

template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const double &, const double &, double &);
template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const tps &, const tps &, tps &);

// template void Insertion_Pass(CellType &, ss_vect<double> &);
// template void Insertion_Pass(CellType &, ss_vect<tps> &);


double d_sign(double a, double b)
{
  double x;

  x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}

int P_eof(FILE *f)
{
  register int ch;

  if (feof(f)) return 1;
  if (f == stdin) return 0; /* not safe to look-ahead on the keyboard! */
  ch = getc(f);
  if (ch == EOF) return 1;
  ungetc(ch, f);

  return 0;
}


/* Check if at end of line (or end of entire file). */

int P_eoln(FILE *f)
{
  register int ch;

  ch = getc(f);
  if (ch == EOF) return 1;
  ungetc(ch, f);
  return (ch == '\n');
}
