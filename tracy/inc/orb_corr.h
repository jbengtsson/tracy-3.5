#ifndef ORB_CORR_H
#define ORB_CORR_H

class orb_corr_type {
private:
  double                **A, **Ai, *w, **U, **V, *b, *bb, *x, *xx, eps, hcut, vcut;

public:
  bool                  hor, periodic;
  int                   m, n;
  std::vector<long int> bpms, corrs;

  void alloc(const long int i0, const long int i1, const long int i2,
	     const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names,
	     const bool hor, const bool periodic, const double hcut, const double vcut, const double eps);
  void alloc(const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names,
	     const bool hor, const bool periodic, const double hcut, const double vcut, const double eps);
  void dealloc(void);

  void get_trm_mat(void);
  void get_orm_mat(void);
  void svd_decomp(void);
  void prt_svdmat(void);
  void solve(const double scl, const double h_maxkick, const double v_maxkick, const long n_bits) const;
  void clr_trims(void);
};


void codstat(double mean[], double sigma[], double xmax[], const long lastpos,
	     const bool all, const std::vector<long int> &bpms);

void thread_beam(const int n_cell, const string &Fam_name,
		 const std::vector<string> &bpm_Fam_names,
		 const std::vector<string> corr_Fam_names[],
		 const int n_orbit, const double scl,
		 const double h_maxkick, const double v_maxkick, const long n_bits, const double hcut, const double vcut);

void cod_ini(const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> corr_Fam_names[],
	     const double hcut,const double vcut, orb_corr_type orb_corr[]);

bool cod_correct(const int n_orbit, const double scl, const double h_maxkick, const double v_maxkick, const long n_bits,
		 orb_corr_type orb_corr[]);

#endif
