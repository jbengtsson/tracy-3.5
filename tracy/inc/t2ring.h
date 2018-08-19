/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -                                         */

void GetNu(Vector2 &nu, Matrix &M);

void Cell_GetABGN(Matrix &M,
		  Vector2 &alpha, Vector2 &beta, Vector2 &gamma, Vector2 &nu);

void Cell_Geteta(long i0, long i1, bool ring, double dP);

void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP);

void Ring_Getchrom(double dP);

void Ring_GetTwiss(bool chroma, double dP);
