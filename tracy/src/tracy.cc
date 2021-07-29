/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#include "tracy_lib.h"

#include "field.cc"

#if NO_TPSA == 1
  // linear TPSA
  #include "tpsa_lin.cc"
  #include "tpsa_lin_pm.cc"
#else
  // interface to M. Berz' TPSA
  #include "tpsa_for_pm.cc"
#endif

// #include "mathlib.cc"
#include "ety.cc"
#include "eigenv.cc"

#include "t2elem.cc"
#include "t2cell.cc"
#include "t2ring.cc"
// #include "sigma_track.cc"

#include "t2lat.cc"
#include "prtmfile.cc"
#include "rdmfile.cc"

#include "set_errors.cc"
#include "lsoc.cc"

// #include "orb_corr.cc"
// #include "param.cc"
// #include "dynap.cc"

// ALS.
// #include "physlib.cc"
// #include "fft.cc"

// Soleil.
#include "radia2tracy.cc"
// #include "naffutils.cc"
// #include "modnaff.cc"
// #include "soleillib.cc"

// NSLS-II.
// #include "nsls-ii_lib.cc"

#if 0
#include "../../python/src/tracy_py.cc"
#endif

str80
  finame,               // input data file.
  foname,               // output data file.
  fname;                // temp file name.

FILE
  *fi,                  // lattice input  file.
  *fo,                  // lattice output file.
  *psin[maxincl+1],     // program input file.
  *psout,               // program output file.
  *prr[maxfil-2];       // prr[1] : input, prr[2] : output.


int
  Fnum_Cart,
  n_iter_Cart;

double
  u_Touschek,           // argument for Touschek D(ksi)
  chi_m;                // argument for IBS D(ksi)

// IBS (Bjorken-Mtingwa)
double a_IBS, b_IBS, c_IBS, a_k_IBS, b_k_IBS;

// for IBS
int    i_, j_;
double **C_;

// ss_vect<tps> map;
// MNF_struct MNF;

// Truncated Power Series Algebra (TPSA)
const int
  nv_tps   = ps_dim,    // No of variables.
  nd_tps   = 3,         // No of degrees of freedom.
  iref_tps = 0;         /* File with resonances to be excluded from the map
			   normal form: fort.7. */

double eps_tps = 1e-25; // Floating point truncation.


// Instantiate templates.

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

template void p_rot(ConfigType &conf, double, ss_vect<double> &);
template void p_rot(ConfigType &conf, double, ss_vect<tps> &);


template void get_B2(const double, const double [], const ss_vect<double> &,
		     double &, double &);
template void get_B2(const double, const tps [], const ss_vect<tps> &,
		     tps &, tps &);

template void radiate(ConfigType &conf, ss_vect<double> &, const double,
		      const double, const double []);
template void radiate(ConfigType &conf, ss_vect<tps> &, const double,
		      const double, const tps []);

template void radiate_ID(ConfigType &conf, ss_vect<double> &,
			 const double, const double &);
template void radiate_ID(ConfigType &conf, ss_vect<tps> &,
			 const double, const tps &);

template void Drift(ConfigType &conf, const double, ss_vect<double> &);
template void Drift(ConfigType &conf, const double, ss_vect<tps> &);

template void bend_fringe(ConfigType &conf, const double, ss_vect<double> &);
template void bend_fringe(ConfigType &conf, const double, ss_vect<tps> &);

template void EdgeFocus(ConfigType &conf, const double, const double,
			const double, ss_vect<double> &);
template void EdgeFocus(ConfigType &conf, const double, const double,
			const double, ss_vect<tps> &);

template void quad_fringe(ConfigType &conf, const double, ss_vect<double> &);
template void quad_fringe(ConfigType &conf, const double, ss_vect<tps> &);


template void thin_kick(ConfigType &conf, const int, const MpoleArray &,
			const double, const double, const double,
			ss_vect<double> &);
template void thin_kick(ConfigType &conf, const int, const MpoleArray &,
			const double, const double, const double,
			ss_vect<tps> &);

template void Cav_Focus(const double L, const double delta, const bool entrance,
			ss_vect<double> &ps);
template void Cav_Focus(const double L, const double delta, const bool entrance,
			ss_vect<tps> &ps);
template void Cav_Focus(const double L, const tps delta, const bool entrance,
			ss_vect<tps> &ps);

template void Wiggler_pass_EF(ConfigType &conf, const ElemType *elem,
			      ss_vect<double> &x);
template void Wiggler_pass_EF(ConfigType &conf, const ElemType *elem,
			      ss_vect<tps> &x);

template void Wiggler_pass_EF2(ConfigType &conf, int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<double> &x);
template void Wiggler_pass_EF2(ConfigType &conf, int nstep, double L,
			       double kxV, double kxH, double kz,
			       double BoBrhoV, double BoBrhoH, double phi,
			       ss_vect<tps> &x);

template void Wiggler_pass_EF3(ConfigType &conf, ElemType *Cell,
			       ss_vect<double> &x);
template void Wiggler_pass_EF3(ConfigType &conf, ElemType *Cell,
			       ss_vect<tps> &x);

template void sol_pass(ConfigType &conf, const ElemType *, ss_vect<double> &);
template void sol_pass(ConfigType &conf, const ElemType *, ss_vect<tps> &);

template void LinearInterpolation2(double &, double &, double &, double &,
				   double &, ElemType *, bool &, int);
template void LinearInterpolation2(tps &, tps &, tps &, tps &, tps &,
				   ElemType *, bool &, int);

template void SplineInterpolation2(double &, double &, double &, double &,
				   ElemType *, bool &);
template void SplineInterpolation2(tps &, tps &, tps &, tps &,
				   ElemType *, bool &);

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

int no_tps = 1;


double d_sign(double a, double b)
{
  double x;

  x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}

// Pascal file I/O (legacy).

int P_eof(FILE *f)
{
  int ch;

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
  int ch;

  ch = getc(f);
  if (ch == EOF) return 1;
  ungetc(ch, f);
  return (ch == '\n');
}


// C++ file I/O.

void file_rd(std::ifstream &inf, const __gnu_debug::string &file_name)
{
  inf.open(file_name.c_str(), std::ios::in);
  if (!inf.is_open()) {
    cout << "File not found: " << file_name << "\n";
    exit_(-1);
  }
}


void file_wr(std::ofstream &outf, const __gnu_debug::string &file_name)
{
  outf.open(file_name.c_str(), std::ios::out);
  if (!outf.is_open()) {
    cout << "Could not create file: " << file_name << "\n";
    exit_(-1);
  }
}


void file_rd(std::ifstream &inf, const char file_name[])
{
  inf.open(file_name, std::ios::in);
  if (!inf.is_open()) {
    printf("File not found: %s\n", file_name);
    exit_(-1);
  }
}


void file_wr(std::ofstream &outf, const char file_name[])
{
  outf.open(file_name, std::ios::out);
  if (!outf.is_open()) {
    printf("Could not create file: %s\n", file_name);
    exit_(-1);
  }
}


// C file I/O.

FILE* file_read(const char file_name[])
{
  FILE *fp;
  fp = fopen(file_name, "r");
  if (fp == NULL) {
    printf("File not found: %s\n", file_name);
    exit_(-1);
  }
  return(fp);
}


FILE* file_write(const char file_name[])
{
  FILE *fp;
  fp = fopen(file_name, "w");
  if (fp == NULL) {
    printf("Could not create file: %s\n", file_name);
    exit_(-1);
  }
  return(fp);
}

void t2init(void)
{
//  iniranf(0); /* initialise le generateur aleatoire: graine 0 */

//  fprintf(stdout,"pi = %20.16e \n",pi);

//  daini((long)no_, (long)nv_, 0);

//  lieini((long)no_, (long)nv_, (long)nd2_);
}


// Matlab BS
void exit_(int exit_code)
{

  printf("fatal error, <ret> to continue "); std::cin.ignore(1, '\n');

  exit(exit_code);
}


double xabs(long n, ss_vect<double> &x)
{
  long    i;
  double  sum;

  sum = 0.0;
  for (i = 0; i < n; i++)
    sum += sqr(x[i]);

  return sqrt(sum);
}


void prt_name_ascii(string &name)
{
  int i;

  printf("  %-8s (", name.c_str());
  for (i = 0; i < (int)name.length(); i++)
    printf(" %3d", (int)name[i]);
  printf(" )\n");
}


long LatticeType::ElemIndex(const __gnu_debug::string &name)
{
  long        i;
  __gnu_debug::string name1 = name;

  const bool prt = false;

  for (i = 0; i < (int)name1.length(); i++)
    name1[i] = tolower(name1[i]);

  if (prt) {
    printf("\nElemIndex:\n     ");
    prt_name_ascii(name1);
    printf("\n");
  }

  for (i = 1; i <= (int)elemf.size(); i++) {
    if (prt) {
      printf("  %3ld", i);
      prt_name_ascii(elemf[i-1].ElemF->Name);
    }

    if (name1 == elemf[i-1].ElemF->Name) break;
  }

  if (name1 != elemf[i-1].ElemF->Name) {
    printf("ElemIndex: undefined element\n");
    exit_(1);
  }

  return i;
}
