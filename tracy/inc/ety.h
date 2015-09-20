/* Tracy-2

   J. Bengtsson  CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

void ETY(int n, int low, int high, Matrix &a, Vector &ort);

void ETYT(int n, int low, int high, Matrix &a, Vector &ort,
                 Matrix &z);

void ety2(int n, int low, int high, Matrix &h, Vector &wr,
	  Vector &wi, Matrix &z, int &ierr);


