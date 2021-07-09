
#ifndef SET_ERRORS
#define SET_ERRORS


enum error_type
  { bn_des   = 0,
    dbn_des  = 1,
    bnL_des  = 2,
    dbnL_des = 3,
    bnL_sys  = 4,
    bnL_rms  = 5,
    bnr_sys  = 6,
    bnr_rms  = 7 };


void iniranf(const long i);

void newseed(void);

double ranf(void);

void setrancut(const double cut);

double normranf(void);

void CheckAlignTol(LatticeType &lat, const char *OutputFile);

void misalign_rms_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd);

void misalign_sys_elem(LatticeType &lat, const int Fnum, const int Knum,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys);

void misalign_rms_fam(LatticeType &lat, const int Fnum,
		      const double dx_rms, const double dy_rms,
		      const double dr_rms, const bool new_rnd);

void misalign_sys_fam(LatticeType &lat, const int Fnum,
		      const double dx_sys, const double dy_sys,
		      const double dr_sys);

void misalign_rms_type(LatticeType &lat, const int type,
		       const double dx_rms, const double dy_rms,
		       const double dr_rms, const bool new_rnd);

void misalign_sys_type(LatticeType &lat, const int type,
		       const double dx_sys, const double dy_sys,
		       const double dr_sys);

void misalign_rms_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_rms, const double dy_rms,
			  const double dr_rms, const bool new_rnd);

void misalign_sys_girders(LatticeType &lat, const int gs, const int ge,
			  const double dx_sys, const double dy_sys,
			  const double dr_sys);

void set_aper_elem(LatticeType &lat, const int Fnum, const int Knum,
		   const double Dxmin, const double Dxmax,
		   const double Dymin, const double Dymax);

void set_aper_fam(LatticeType &lat, const int Fnum,
		  const double Dxmin, const double Dxmax,
		  const double Dymin, const double Dymax);

void set_aper_type(LatticeType &lat, const int type, const double Dxmin,
		   const double Dxmax, const double Dymin, const double Dymax);

double get_L(LatticeType &lat, const int Fnum, const int Knum);

void set_L(LatticeType &lat, const int Fnum, const int Knum, const double L);
void set_L(LatticeType &lat, const int Fnum, const double L);
void set_dL(LatticeType &lat, const int Fnum, const int Knum, const double dL);
void set_dL(LatticeType &lat, const int Fnum, const double dL);

//------------------------------------------------------------------------------

void get_bn_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			const int n, double &bn, double &an);
void get_bnL_design_elem(LatticeType &lat, const int Fnum, const int Knum,
			 const int n, double &bnL, double &anL);

void set_bn(LatticeType &lat, error_type err, const int Fnum, const int Knum,
	    const int n, const double bn, const double an, const bool new_rnd);

void set_bn(LatticeType &lat, error_type err, const int Fnum, const int n,
	    const double bn, const double an, const bool new_rnd);

//------------------------------------------------------------------------------

void set_bnL_design_type(LatticeType &lat, const int type, const int n,
			 const double bnL, const double anL);

void set_bnL_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL);
void set_bnr_sys_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr);

void set_bnL_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnL, const double anL, const bool new_rnd);

void set_bnr_rms_type(LatticeType &lat, const int type, const int n,
		      const double bnr, const double anr,
		      const bool new_rnd);

#endif
