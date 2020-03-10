/* Tracy-2

   J. Bengtsson  CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#ifndef EIGENV_H
#define EIGENV_H

void geigen(int n, Matrix &fm, Matrix &Vre, Matrix &Vim,
	    psVector &wr, psVector &wi);

#endif
