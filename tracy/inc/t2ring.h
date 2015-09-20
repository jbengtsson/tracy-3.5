/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

void GetNu(Vector2 &nu, Matrix &M);

void Cell_Twiss(long i0, long i1, ss_vect<tps> &Ascr, bool chroma, bool ring,
		double dP);

void Cell_Twiss_M(long i0, long i1, Matrix &Ascr, bool chroma, bool ring,
		  double dP);

void Cell_GetABGN(Matrix &M,
		  Vector2 &alpha, Vector2 &beta, Vector2 &gamma, Vector2 &nu);

void Cell_Geteta(long i0, long i1, bool ring, double dP);

void Ring_Getchrom(double dP);

void Ring_GetTwiss(bool chroma, double dP);


void Ring_Fittune(Vector2 &nu, double eps, iVector2 &nq, long qf[], long qd[],
                  double dkL, long imax);

void Ring_Fitchrom(Vector2 &ksi, double eps, iVector2 &ns,
		   long sf[], long sd[], double dkpL, long imax);

void Ring_FitDisp(long pos, double eta, double eps, long nq, long q[],
                  double dkL, long imax);


void TransTrace(long i0, long i1, Vector2 &alpha, Vector2 &beta, Vector2 &eta,
                Vector2 &etap, Vector &codvect);
