/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

#ifndef LSOC_H
#define LSOC_H

extern int               n_bpm_[2], n_corr_[2];
extern long unsigned int *bpms_[2], *corrs_[2];

void zero_trims(LatticeType &lat);

void prt_gcmat(const int plane);

void gcmat(LatticeType &lat, const int plane);

void gcmat(LatticeType &lat, const int n_bpm, const long int bpms[],
	   const int n_corr, const long int corrs[], const int plane,
	   const bool svd);

void gcmat(LatticeType &lat, const int bpm, const int corr,
	   const int plane);

void lsoc(LatticeType &lat, const int plane, const double scl);

void gtcmat(LatticeType &lat, const int plane);

void gtcmat(LatticeType &lat, const int n_bpm, const long int bpms[],
	    const int n_corr, const long int corrs[], const int plane,
	    const bool svd);

void lstc(LatticeType &lat, const int plane, const double scl);

#endif
