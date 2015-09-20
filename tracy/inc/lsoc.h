/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

extern int                n_bpm_[2], n_corr_[2];
extern long unsigned int  *bpms_[2], *corrs_[2];

void zero_trims(void);

void prt_gcmat(const int plane);

void gcmat(const int plane);

void gcmat(const int n_bpm, const long int bpms[],
	   const int n_corr, const long int corrs[], const int plane,
	   const bool svd);

void gcmat(const int bpm, const int corr, const int plane);

void lsoc(const int plane);

void gtcmat(const int plane);

void gtcmat(const int n_bpm, const long int bpms[],
           const int n_corr, const long int corrs[], const int plane,
           const bool svd);

void lstc(const int plane, const long int lastpos);

