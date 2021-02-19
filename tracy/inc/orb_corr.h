#ifndef ORB_CORR_H
#define ORB_CORR_H

class orb_corr_type {
private:
  double **A, **Ai, *w, **U, **V, *b, *bb, *x, *xx, eps, hcut, vcut;

public:
  bool                  hor, periodic;
  int                   m, n;
  std::vector<long int> bpms, corrs;

  void alloc(LatticeType &lat, const long int i0, const long int i1,
	     const long int i2, const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names, const bool hor,
	     const bool periodic, const double eps);
  void alloc(LatticeType &lat, const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> &corr_Fam_names, const bool hor,
	     const bool periodic, const double eps);
  void dealloc(void);
  void get_trm_mat(LatticeType &lat);
  void get_orm_mat(LatticeType &lat);
  void svd_decomp(void);
  void prt_svdmat(LatticeType &lat);
  void solve(LatticeType &lat, const double scl) const;
  void clr_trims(LatticeType &lat);
};


void codstat(LatticeType &lat, double mean[], double sigma[], double xmax[],
	     const long lastpos, const bool all,
	     const std::vector<long int> &bpms);

void thread_beam(LatticeType &lat, const int n_cell, const string &Fam_name,
		 const std::vector<string> &bpm_Fam_names,
		 const std::vector<string> corr_Fam_names[],
		 const int n_thread, const double scl);

void cod_ini(LatticeType &lat, const std::vector<string> &bpm_Fam_names,
	     const std::vector<string> corr_Fam_names[],
	     orb_corr_type orb_corr[]);

bool cod_correct(LatticeType &lat, const int n_orbit, const double scl,
		 orb_corr_type orb_corr[]);

#endif
