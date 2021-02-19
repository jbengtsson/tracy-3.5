/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

#ifndef SOLEILLIB_H
#define SOLEILLIB_H

/**** Protypes ****/
void SetErr(LatticeType &lat);
void InducedAmplitude(LatticeType &lat, long spos);
void Hfonction(LatticeType &lat, long pos, double dP,Vector2 H);
void Hcofonction(LatticeType &lat, long pos, double dP,Vector2 H);
void Get_Disp_dp(LatticeType &lat);
void read_corrh(void);
void set_vectorcod(LatticeType &lat, psVector codvector[], double dP);
void SetDecapole(void);

/* Tracking */
void Phase(LatticeType &lat, double x,double xp,double y, double yp,
	   double energy, double ctau, long Nbtour);
void Phase2(LatticeType &lat, long pos, double x,double xp,double y, double yp,
	    double energy, double ctau, long Nbtour);
void PhasePoly(LatticeType &lat, long pos, double x0,double px0, double z0,
	       double pz0, double delta0, double ctau0, long Nbtour);
void Check_Trac(LatticeType &lat, double x, double px, double y, double py,
		double dp);
void PhasePortrait(LatticeType &lat, double x0,double px0,double z0, double pz0,
		   double delta0, double ctau, double end, long Nb, long Nbtour,
		   int num);
void PhasePortrait2(long pos,double x0,double px0,double z0, double pz0,
		    double delta0, double ctau,
		    double end, long Nb, long Nbtour, int num);
void Multipole(LatticeType &lat);
void MomentumAcceptance(LatticeType &lat, long deb, long fin, double ep_min,
			double ep_max, long nstepp, double em_min,
			double em_max, long nstepm);
void Trac_Tab(LatticeType &lat, double x, double px, double y, double py,
	      double dp, long nmax, long pos, long *lastn, long *lastpos,
	      FILE *outf1, double Tx[][NTURN]);
void SetSkewQuad(LatticeType &lat);
void TracCO(LatticeType &lat, double x, double px, double y, double py,
	    double dp, double ctau, long nmax, long pos, long *lastn,
	    long *lastpos, FILE *outf1);
void Dyna(LatticeType &lat, long Nbx, long Nbz, long Nbtour, double xmax,
	  double zmax, double energy, bool diffusion);
                            
/* Frequency map analysis */
void NuDp(LatticeType &lat, long Nb, long Nbtour, double emax);
void fmap(LatticeType &lat, long Nbx, long Nbz, long Nbtour, double xmax,
	  double zmax, double energy, bool diffusion, bool matlab);
void fmapdp(LatticeType &lat, long Nbx, long Nbe, long Nbtour, double xmax,
	    double emax, double z, bool diffusion, bool matlab);
void Nu_Naff(void);
void NuDx(LatticeType &lat, long Nbx, long Nbz, long Nbtour, double xmax,
	  double ymax, double energy);

/* Vacuum chamber */
void DefineCh(LatticeType &lat);
void Enveloppe(LatticeType &lat, double x, double px, double y, double py,
	       double dp, double nturn);
void ChamberOn(LatticeType &lat);

/* Longitudinal Hamiltonian*/
void PhaseLongitudinalHamiltonien(void);
void PassA(double *phi, double delta0, double step);
void PassB(double phi0, double *delta, double step);
double Hsynchrotron(double phi, double delta);

/* Miscelleneous */ 
void Enveloppe2(LatticeType &lat, double x, double px, double y, double py,
		double dp, double nturn);
void Phase3(LatticeType &lat, long pos, double x,double px,double y, double py,
	    double energy, double ctau, long Nbtour);
double EnergySmall(double *X, double irho);
double EnergyDrift(double *X);
void getA4antidamping(LatticeType &lat);
void fmapfull(LatticeType &lat, long Nbx, long Nbz, long Nbtour, double xmax,
	      double zmax, double energy, bool diffusion);
void spectrum(LatticeType &lat, long Nbx, long Nbz, long Nbtour, double xmax,
	      double zmax, double energy, bool diffusion);

#endif
