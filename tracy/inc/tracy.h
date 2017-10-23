/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


#ifndef LONG_MAX
# define LONG_MAX ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN (~LONG_MAX)
#endif

#define S_SIZE 200  // max size for file name of a lattice file

#define blankname    "               "

#define maxincl      5
#define maxfil       10
#define bigvectmax   4096

typedef char   str80[80];

extern bool  stable;
extern bool  ErrFlag;
extern bool  trace, traceID;

extern str80  finame,   /* input  data file  */
              foname,   /* output data file */
              fname;    /* temp file name */

extern FILE  *fi,       /* lattice input  file  */
             *fo,       /* lattice output file */
             *psin[],   /* program input file */
             *psout,              /* program output file*/
             *prr[];   /* prr[1] : input, prr[2] : output */

extern int P_eof(FILE *f);

extern int P_eoln(FILE *f);

/**** same as asctime in C without the \n at the end****/
char *asctime2(const struct tm *timeptr);
struct tm* GetTime();
