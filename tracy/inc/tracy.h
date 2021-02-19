/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

#ifndef TRACY_H
#define TRACY_H

// Defined in math.h.
//#define M_PI   3.14159265358979323846  // pi

#ifndef LONG_MAX
# define LONG_MAX   ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN   (~LONG_MAX)
#endif

#define S_SIZE 200  // max size for file name of a lattice file

#define NameLength      150  // maximum length of identifiers (e.g. file names)
#define SymbolLength    15   // maximum length of element name

#define blankname       "               "

#define maxincl         5
#define maxfil          10
#define bigvectmax      4096

// Physics constants
const double  c0    = 2.99792458e8;             // speed of light in vacuum
const double  q_e   = 1.602e-19;                // electron charge
const double  m_e   = 0.51099906e6;             // electron rest mass [eV/c^2]
const double  mu_0  = 4.0*M_PI*1e-7;            // permittivity of free space
const double  eps_0 = 1.0/(sqr(c0)*mu_0);       // permeability of free space
const double  r_e   = q_e/(4.0*M_PI*eps_0*m_e); // classical electron radius
const double  h_bar = 6.58211899e-16;           /* reduced Planck constant
						   [eV s] */
typedef char str80[80];

typedef char alfa_[NameLength];

typedef long   iVector2[2];
typedef double Vector2[2];
typedef double Vector3[3];

#define fitvectmax      200
typedef long   fitvect[fitvectmax];

extern bool  stable;
extern bool  ErrFlag;
extern bool  trace, traceID;

extern double  Fdrift1, Fkick1, Fdrift2, Fkick2, crad, cfluc;

extern str80  finame,   /* input  data file  */
              foname,   /* output data file */
              fname;    /* temp file name */

extern FILE  *fi,       /* lattice input  file  */
             *fo,       /* lattice output file */
             *psin[],   /* program input file */
             *psout,              /* program output file*/
             *prr[];   /* prr[1] : input, prr[2] : output */

extern bool reverse_elem;

extern int P_eof(FILE *f);

extern int P_eoln(FILE *f);

extern void GDiag(int n, double C, Matrix &A, Matrix &Ainv, Matrix &R,
		  Matrix &M, double &Omega, double &Yalphac);
extern void NormEigenVec(Matrix &Vr, Matrix &Vi, double *wr, double *wi,
			 Matrix &t6a);

extern void t2init(void);

extern void prt_gcmat(int bpm, int corr, int plane);

extern void gcmat(int bpm, int corr, int plane);

extern void lsoc(int niter, int bpm, int corr, int plane);

/**** same as asctime in C without the \n at the end****/
char *asctime2(const struct tm *timeptr);
struct tm* GetTime();

#endif
