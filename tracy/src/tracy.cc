/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

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

double        eps_tps  = 1e-25;  // floating point truncation


// instantiate templates

template class ss_vect<double>;

template class ss_vect<tps>;


template void GtoL(ss_vect<double> &, Vector2 &, Vector2 &,
		   const double, const double, const double);

template void GtoL(ss_vect<tps> &, Vector2 &, Vector2 &,
		   const double, const double, const double);

template void LtoG(ss_vect<tps> &, Vector2 &, Vector2 &,
		   double, double, double);

template void LtoG(ss_vect<double> &, Vector2 &, Vector2 &,
		   double, double, double);

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


template void Drift_Pass(CellType &, ss_vect<double> &);

template void Drift_Pass(CellType &, ss_vect<tps> &);

template void thin_kick(const int, const double [], const double, const double,
			const double, ss_vect<double> &);

template void thin_kick(const int, const double [], const double, const double,
			const double, ss_vect<tps> &);

template void Mpole_Pass(CellType &, ss_vect<double> &);

template void Mpole_Pass(CellType &, ss_vect<tps> &);

template void Marker_Pass(CellType &, ss_vect<double> &);

template void Marker_Pass(CellType &, ss_vect<tps> &);

template void Cav_Focus(const double L, const double delta, const bool entrance,
            ss_vect<double> &ps);

template void Cav_Focus(const double L, const double delta, const bool entrance,
            ss_vect<tps> &ps);

template void Cav_Focus(const double L, const tps delta, const bool entrance,
            ss_vect<tps> &ps);

template void Cav_Pass(CellType &, ss_vect<double> &);

template void Cav_Pass(CellType &, ss_vect<tps> &);

template void Wiggler_pass_EF(const elemtype &elem, ss_vect<double> &x);

template void Wiggler_pass_EF(const elemtype &elem, ss_vect<tps> &x);

template void Wiggler_pass_EF2(int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<double> &x);

template void Wiggler_pass_EF2(int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<tps> &x);

template void Wiggler_pass_EF3(const elemtype &elem, ss_vect<double> &x);

template void Wiggler_pass_EF3(const elemtype &elem, ss_vect<tps> &x);

template void Wiggler_Pass(CellType &, ss_vect<double> &);

template void Wiggler_Pass(CellType &, ss_vect<tps> &);

template void FieldMap_Pass(CellType &, ss_vect<double> &);

template void FieldMap_Pass(CellType &, ss_vect<tps> &);

template void sol_pass(const elemtype &, ss_vect<double> &);

template void sol_pass(const elemtype &, ss_vect<tps> &);

template void Solenoid_Pass(CellType &, ss_vect<double> &);

template void Solenoid_Pass(CellType &, ss_vect<tps> &);

template void LinearInterpolation2(double &, double &, double &, double &,
				   double &, CellType &, bool &, int);

template void LinearInterpolation2(tps &, tps &, tps &, tps &, tps &,
				   CellType &, bool &, int);

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

template void Insertion_Pass(CellType &, ss_vect<double> &);

template void Insertion_Pass(CellType &, ss_vect<tps> &);


template void Elem_Pass(const long, ss_vect<double> &);

template void Elem_Pass(const long, ss_vect<tps> &);


template void Cell_Pass(const long, const long, ss_vect<double> &, long &);

template void Cell_Pass(const long, const long, ss_vect<tps> &, long &);


/* Global variable used through the code */
globvalrec globval;

statusrec status;
bool trace, traceID;
bool cellconcat;

/* Random stuff */
long rseed0, rseed;
double normcut_;

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
