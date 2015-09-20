/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

 /**** Protypes ****/
void SetErr(void);
void InducedAmplitude(long spos);
void Hfonction(long pos, double dP,Vector2 H);
void Hcofonction(long pos, double dP,Vector2 H);
void Get_Disp_dp(void);
void read_corrh(void);
void set_vectorcod(Vector codvector[], double dP);
void SetDecapole(void);

/* Tracking */
void Phase(double x,double xp,double y, double yp,double energy, double ctau, long Nbtour);
void Phase2(long pos, double x,double xp,double y, double yp,double energy, double ctau,
            long Nbtour);
void PhasePoly(long pos, double x0,double px0, double z0, double pz0, double delta0,
               double ctau0, long Nbtour);
void Check_Trac(double x, double px, double y, double py, double dp);
void PhasePortrait(double x0,double px0,double z0, double pz0, double delta0, double ctau,
                          double end, long Nb, long Nbtour, int num);
void PhasePortrait2(long pos,double x0,double px0,double z0, double pz0, double delta0, double ctau,
                          double end, long Nb, long Nbtour, int num);
void Multipole(void);
void MomentumAcceptance(long deb, long fin, double ep_min, double ep_max, long nstepp,
                        double em_min, double em_max, long nstepm);
void Trac_Tab(double x, double px, double y, double py, double dp,
            long nmax, long pos, long *lastn, long *lastpos, FILE *outf1, double Tx[][NTURN]);
void SetSkewQuad(void);
void TracCO(double x, double px, double y, double py, double dp, double ctau,
                 long nmax, long pos, long *lastn, long *lastpos, FILE *outf1);
void Dyna(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
               double energy, bool diffusion);
                            
/* Frequency map analysis */
void NuDp(long Nb, long Nbtour, double emax);
void fmap(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
                 double energy, bool diffusion, bool matlab);
void fmapdp(long Nbx, long Nbe, long Nbtour, double xmax, double emax,
              double z, bool diffusion, bool matlab);
void Nu_Naff(void);
void NuDx(long Nbx, long Nbz, long Nbtour, double xmax, double ymax,
                 double energy);

/* Vacuum chamber */
void DefineCh(void);
void Enveloppe(double x, double px, double y, double py,
                      double dp, double nturn);
void ChamberOn(void);

/* Longitudinal Hamiltonian*/
void PhaseLongitudinalHamiltonien(void);
void PassA(double *phi, double delta0, double step);
void PassB(double phi0, double *delta, double step);
double Hsynchrotron(double phi, double delta);

/* Miscelleneous */ 
void Enveloppe2(double x, double px, double y, double py,
                      double dp, double nturn);
void Phase3(long pos, double x,double px,double y, double py,double energy,
            double ctau, long Nbtour);
double EnergySmall(double *X, double irho);
double EnergyDrift(double *X);
void getA4antidamping();
void fmapfull(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion);
void spectrum(long Nbx, long Nbz, long Nbtour, double xmax, double zmax,
              double energy, bool diffusion);
