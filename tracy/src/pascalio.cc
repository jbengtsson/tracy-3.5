/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

 pascalio.c -- PASCAL i/o routines and initialization routines

*/


str80 finame;   /* input  data file  */
str80 foname;   /* output data file */
str80 fname;   /* temp file name */

FILE *fi;   /* lattice input  file  */
FILE *fo;   /* lattice output file */
FILE *psin[maxincl + 1];   /* program input file */
FILE *psout;   /* program output file*/
FILE *prr[maxfil - 2];   /* prr[1] : input, prr[2] : output */

bool ErrFlag;


void t2init(void)
{
//  iniranf(0); /* initialise le generateur aleatoire: graine 0 */

//  fprintf(stdout,"pi = %20.16e \n",pi);

  cellconcat = false;

//  daini((long)no_, (long)nv_, 0);

//  lieini((long)no_, (long)nv_, (long)nd2_);

  globval.MatMeth = false;
}


// Matlab BS
void exit_(int exit_code)
{

  printf("fatal error, <ret> to continue "); cin.ignore(1, '\n');

  exit(exit_code);
}
