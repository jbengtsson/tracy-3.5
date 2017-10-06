#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

// Initial conditions: alpha, beta, eta, etap.

// Upstream of B20.
// From k_symm.lat.
// const double ic[][2] =
//   {{1.05266, -0.25384}, {0.62733, 5.60502}, {0.06552, 0.0}, {-0.10478, 0.0}};
// Roughly, what's obtained for jb_2.lat:
// const double ic[][2] =
//   {{1.15199, -0.22236}, {0.65878, 5.53043}, {0.03741, 0.0}, {-0.04304, 0.0}};

// Upstream of QF03.
// From k_symm.lat.
// const double ic[][2] =
//   {{-6.00257, 2.31594}, {7.43795, 2.74922}, {0.18750, 0.0}, {0.16825, 0.0}};
// From jb_2.lat.
const double ic[][2] =
  {{-6.57117, 2.40476}, {8.00696, 2.78072}, {0.08540, 0.0}, {0.08387, 0.0}};


const bool qf031 = true;

int loc[10], n, n_strength;

double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max);
double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max);
double get_bn_s(const int Fnum, const int Knum);
void set_bn_s(const int Fnum, const double dbn);
void get_S(void);


struct param_type {
private:

public:
  int                 n_prm;
  double              bn_tol, step;
  std::vector<double> bn_min, bn_max, bn_scl;
  std::vector<int>    Fnum, n;

  void add_prm(const std::string Fname, const int n,
	       const double bn_min, const double bn_max,
	       const double bn_scl);
  void ini_prm(double *bn);
  void set_prm(double *bn) const;
  void prt_prm(double *bn) const;
};


param_type b2_prms;


void param_type::add_prm(const std::string Fname, const int n,
			 const double bn_min, const double bn_max,
			 const double bn_scl)
{
  Fnum.push_back(ElemIndex(Fname.c_str()));
  this->n.push_back(n);
  this->bn_min.push_back(bn_min);
  this->bn_max.push_back(bn_max);
  this->bn_scl.push_back(bn_scl);
  n_prm = Fnum.size();
}


void param_type::ini_prm(double *bn)
{
  int    i;
  double an;

  n_prm = Fnum.size();
  for (i = 1; i <= n_prm; i++) {
    if (n[i-1] > 0)
      // Multipole.
      get_bn_design_elem(Fnum[i-1], 1, n[i-1], bn[i], an);
    else if (n[i-1] == -1)
      // Drift.
      bn[i] = get_L(Fnum[i-1], 1);
    else if (n[i-1] == -2)
      // Placement.
      bn[i] = get_bn_s(Fnum[i-1], 1);
    // Bounded.
    bn[i] = bn_internal(bn[i], bn_min[i-1], bn_max[i-1]);
  }
}


void param_type::set_prm(double *bn) const
{
  int    i;
  double bn_ext;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    if (n[i-1] > 0)
	set_bn_design_fam(Fnum[i-1], n[i-1], bn_ext, 0e0);
    else if (n[i-1] == -1) {
      set_L(Fnum[i-1], bn_ext); get_S();
    } else if (n[i-1] == -2)
      set_bn_s(Fnum[i-1], bn_ext);
  }
}


void param_type::prt_prm(double *bn) const
{
  int    i;
  double bn_ext;

  const int n_prt = 8;

  for (i = 1; i <= n_prm; i++) {
    // Bounded.
    bn_ext = bn_bounded(bn[i], bn_min[i-1], bn_max[i-1]);
    printf(" %9.5f", bn_ext);
    if (i % n_prt == 0) printf("\n");
  }
  if (n_prm % n_prt != 0) printf("\n");
}


double bn_internal(const double bn_bounded,
		   const double bn_min, const double bn_max)
{
  return asin((2e0*(bn_bounded-bn_min))/(bn_max-bn_min)-1e0);
}


double bn_bounded(const double bn_internal,
		  const double bn_min, const double bn_max)
{
  return bn_min + (sin(bn_internal)+1e0)*(bn_max-bn_min)/2e0;
}


void get_s_loc(const int Fnum, const int Knum, int loc[])
{
  partsName name;

  const bool prt = false;

  // Point to multipole.
  loc[1] = Elem_GetPos(Fnum, Knum);
  if (Lattice.Cell[loc[1]-1].Elem.PName[1] == 'u') {
    loc[0] = loc[1] - 1;
    strcpy(name, Lattice.Cell[loc[1]-1].Elem.PName); name[1] = 'd';
    loc[2] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Lattice.Cell[loc[1]-1].Elem.PName[1] == 'd') {
    loc[2] = loc[1] - 1;
    strcpy(name, Lattice.Cell[loc[1]-1].Elem.PName); name[1] = 'u';
    loc[0] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Lattice.Cell[loc[1]+1].Elem.PName[1] == 'd') {
    loc[2] = loc[1] + 1;
    strcpy(name, Lattice.Cell[loc[1]+1].Elem.PName); name[1] = 'u';
    loc[0] = Elem_GetPos(ElemIndex(name), Knum);
  } else if (Lattice.Cell[loc[1]+1].Elem.PName[1] == 'u') {
    loc[0] = loc[1] + 1;
    strcpy(name, Lattice.Cell[loc[1]+1].Elem.PName); name[1] = 'd';
    loc[2] = Elem_GetPos(ElemIndex(name), Knum);
  } else {
    printf("\nget_s_loc: configuration error %s (%d)\n",
	   Lattice.Cell[loc[1]].Elem.PName, loc[1]);
    exit(1);
  }

  if (prt)
    printf("\nget_s_loc: %s %s %s\n",
	   Lattice.Cell[loc[0]].Elem.PName, Lattice.Cell[loc[1]].Elem.PName,
	   Lattice.Cell[loc[2]].Elem.PName);
}


double get_bn_s(const int Fnum, const int Knum)
{
  int    loc[3];
  double ds;

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  ds = Lattice.Cell[loc[0]].Elem.PL;

  if (prt)
    printf("\nget_bn_s:  %s %s(%d) %s %10.3e %10.3e\n",
	   Lattice.Cell[loc[0]].Elem.PName, Lattice.Cell[loc[1]].Elem.PName, Knum,
	   Lattice.Cell[loc[2]].Elem.PName, Lattice.Cell[loc[0]].Elem.PL, Lattice.Cell[loc[2]].Elem.PL);

  return ds;
}


void set_bn_s(const int Fnum, const int Knum, const double ds)
{
  int loc[3];

  const bool prt = false;

  get_s_loc(Fnum, Knum, loc);
  set_L(Lattice.Cell[loc[0]].Fnum, Knum, ds);
  set_L(Lattice.Cell[loc[2]].Fnum, Knum, -ds);

  if (prt)
    printf("\nset_bn_s:  %s %s(%d) %s %10.3e %10.3e\n",
	   Lattice.Cell[loc[0]].Elem.PName, Lattice.Cell[loc[1]].Elem.PName, Knum,
	   Lattice.Cell[loc[2]].Elem.PName, ds, -ds);
}


void set_bn_s(const int Fnum, const double ds)
{
  int k;

  for (k = 1; k <= GetnKid(Fnum); k++)
    set_bn_s(Fnum, k, ds);
}


void get_S(void)
{
  int    j;
  double S;

  S = 0e0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    S += Lattice.Cell[j].Elem.PL; Lattice.Cell[j].S = S;
  }
}


void get_dnu(const int i0, const int i1, const double alpha0[],
	     const double beta0[], const double delta, double dnu[])
{
  long int     lastpos;
  int          k;
  double       m11, m12;
  ss_vect<tps> map;

  map.identity(); map[delta_] += delta;
  Cell_Pass(i0+1, i1, map, lastpos);

  for (k = 0; k < 2; k++) {
    m11 = map[2*k][2*k]; m12 = map[2*k][2*k+1];
    dnu[k] = atan(m12/(beta0[k]*m11-alpha0[k]*m12))/(2e0*M_PI);
    if (m11 < 0e0) dnu[k] += (m12 >= 0)? 0.5e0 : -0.5e0;
    if (dnu[k] < 0e0) dnu[k] += 1e0;
  }
}


void get_dnu_delta(const int i0, const int i1, const double alpha0[],
		const double beta0[], const double delta, double dksi[])
{
  // To evaluate linear chromaticity the linear dispersion for the periodic
  // solution must be known.
  int    k;
  double dnu1[2], dnu0[2];

  get_dnu(loc[0], loc[1], alpha0, beta0, delta, dnu1);
  get_dnu(loc[0], loc[1], alpha0, beta0, -delta, dnu0);

  for (k = 0; k < 2; k++)
    dksi[k] = (dnu1[k]-dnu0[k])/(2e0*delta);
}


void prt_b2(const param_type &b2_prms, const double *b2)
{
  int  k;
  FILE *outf;

  std::string file_name = "match.out";

  outf = file_write(file_name.c_str());

  k = 1;
  fprintf(outf, "QF031: quadrupole, l = 0.217, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "QD041: quadrupole, l = 0.117, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nQ01:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  // k++;
  // fprintf(outf, "Q02:  quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  // 	  ", Method = Meth;\n",
  // 	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "Q03:  quadrupole, l = 0.434, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nEQ01: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ02: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k],
  		     b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  k++;
  fprintf(outf, "\nEQ03: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ04: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  k++;
  fprintf(outf, "EQ05: quadrupole, l = 0.234, k = %8.5f, N = Nquad"
  	  ", Method = Meth;\n",
  	  bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

  if (false) {
    k++;
    fprintf(outf, "\nD_Q01_L  = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_Q02_L  = %8.5f;\n",
    	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_Q03_L  = %8.5f;\n",
    	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    k++;
    fprintf(outf, "\nD_EQ01_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ02_L = %8.5f;\n",
	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    k++;
    fprintf(outf, "\nD_EQ03_L = %8.5f;\n",
    	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ04_L = %8.5f;\n",
    	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
    k++;
    fprintf(outf, "D_EQ05_L = %8.5f;\n",
    	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    // k++;
    // fprintf(outf, "\nD_B10_L  = %8.5f;\n",
    // 	    bn_bounded(b2[k], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));

    // k++;
    // fprintf(outf, "\nU561: drift, L = %8.5f;\n",
    // 	    bn_bounded(b2[18], b2_prms.bn_min[k-1], b2_prms.bn_max[k-1]));
  }

  fclose(outf);
}


void prt_lin_opt(const int loc[])
{
  printf("\n      s    alpha_x  beta_x  eta_x  etap_x  alpha_y  beta_y\n");
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Lattice.Cell[loc[0]].S,
	 Lattice.Cell[loc[0]].Alpha[X_], Lattice.Cell[loc[0]].Beta[X_],
	 Lattice.Cell[loc[0]].Eta[X_], Lattice.Cell[loc[0]].Etap[X_],
	 Lattice.Cell[loc[0]].Alpha[Y_], Lattice.Cell[loc[0]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Lattice.Cell[loc[1]].S,
	 Lattice.Cell[loc[1]].Alpha[X_], Lattice.Cell[loc[1]].Beta[X_],
	 Lattice.Cell[loc[1]].Eta[X_], Lattice.Cell[loc[1]].Etap[X_],
	 Lattice.Cell[loc[1]].Alpha[Y_], Lattice.Cell[loc[1]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Lattice.Cell[loc[2]].S,
	 Lattice.Cell[loc[2]].Alpha[X_], Lattice.Cell[loc[2]].Beta[X_],
	 Lattice.Cell[loc[2]].Eta[X_], Lattice.Cell[loc[2]].Etap[X_],
	 Lattice.Cell[loc[2]].Alpha[Y_], Lattice.Cell[loc[2]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Lattice.Cell[loc[3]].S,
	 Lattice.Cell[loc[3]].Alpha[X_], Lattice.Cell[loc[3]].Beta[X_],
	 Lattice.Cell[loc[3]].Eta[X_], Lattice.Cell[loc[3]].Etap[X_],
	 Lattice.Cell[loc[3]].Alpha[Y_], Lattice.Cell[loc[3]].Beta[Y_]);
  printf("  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
	 Lattice.Cell[loc[4]].S,
	 Lattice.Cell[loc[4]].Alpha[X_], Lattice.Cell[loc[4]].Beta[X_],
	 Lattice.Cell[loc[4]].Eta[X_], Lattice.Cell[loc[4]].Etap[X_],
	 Lattice.Cell[loc[4]].Alpha[Y_], Lattice.Cell[loc[4]].Beta[Y_]);
}


double f_match(double *b2)
{
  static double chi2_ref = 1e30, chi2_prt = 1e30;

  int          i, loc1, loc2;
  double       chi2, dksi[2], L;
  ss_vect<tps> Ascr;

  const int    n_prt     = 50;
  const double eps_delta = 1e-5;

  b2_prms.set_prm(b2);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0], loc[4], Ascr, false, false, 0e0);

  // get_dnu_delta(loc[0], loc[5], ic[0], ic[1], eps_delta, dksi);

  chi2 = 0e0;
  
  // Downstream of 10 degree dipole.
  chi2 += 1e8*sqr(Lattice.Cell[loc[1]].Eta[X_]);
  chi2 += 1e8*sqr(Lattice.Cell[loc[1]].Etap[X_]);
  // chi2 += 1e0*sqr(Lattice.Cell[loc[1]].Beta[Y_]);

  // Center of 1st straight.
  chi2 += 1e5*sqr(Lattice.Cell[loc[2]].Alpha[X_]);
  // chi2 += 1e5*sqr(Lattice.Cell[loc[2]].Alpha[Y_]);
  chi2 += 1e3*sqr(Lattice.Cell[loc[2]].Beta[X_]-9.58);

  // Center of 2nd straight.
  chi2 += 1e5*sqr(Lattice.Cell[loc[3]].Alpha[X_]);
  chi2 += 1e5*sqr(Lattice.Cell[loc[3]].Alpha[Y_]);
  chi2 += 1e3*sqr(Lattice.Cell[loc[3]].Beta[X_]-8.0); 

  for (i = 1; i <= b2_prms.n_prm; i++) {
    loc1 = Elem_GetPos(b2_prms.Fnum[i-1], 1);
    L = Lattice.Cell[loc1].Elem.PL;
    // Need to use internal variable for convergence.
    if (i <= n_strength) {
      chi2 += 1e1*sqr(b2[i]*L*Lattice.Cell[loc1].Beta[X_]);
      chi2 += 1e1*sqr(b2[i]*L*Lattice.Cell[loc1].Beta[Y_]);
     } else {
      chi2 += 1e-10*sqr(b2[i]);
      chi2 += 1e-10*sqr(b2[i]);
   }
  }

  if (chi2 < chi2_ref) {
    n++;

    if (n % n_prt == 0) {
      printf("\n%3d chi2: %12.5e -> %12.5e\n", n, chi2_prt, chi2);
      printf("b2s:\n");
      b2_prms.prt_prm(b2);

      // Downstream of 10 degree dipole.
      printf("\nDownstream of 10 degree dipole:\n");
      printf("eta_x   = %8.5f etap_x  = %8.5f\n",
	     Lattice.Cell[loc[1]].Eta[X_], Lattice.Cell[loc[1]].Etap[X_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     Lattice.Cell[loc[1]].Beta[X_], Lattice.Cell[loc[1]].Beta[Y_]);
      // Center of 1st straight.
      printf("\nCenter of 1st straight:\n");
      printf("alpha_x = %8.5f alpha_y = %8.5f\n",
	     Lattice.Cell[loc[2]].Alpha[X_], Lattice.Cell[loc[2]].Alpha[Y_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     Lattice.Cell[loc[2]].Beta[X_], Lattice.Cell[loc[2]].Beta[Y_]);

      // Center of 2nd straight.
      printf("\nCenter of 2nd straight:\n");
      printf("\nalpha_x = %8.5f alpha_y = %8.5f\n",
	     Lattice.Cell[loc[3]].Alpha[X_], Lattice.Cell[loc[3]].Alpha[Y_]);
      printf("beta_x  = %8.5f beta_y  = %8.5f\n",
	     Lattice.Cell[loc[3]].Beta[X_], Lattice.Cell[loc[3]].Beta[Y_]);

      loc1 = Elem_GetPos(ElemIndex("s_s_1"), 1);
      loc2 = Elem_GetPos(ElemIndex("s_s_1"), 2);
      printf("\nLength of 1st straight: %6.3f m\n", Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);
      loc1 = Elem_GetPos(ElemIndex("s_s_2"), 1);
      loc2 = Elem_GetPos(ElemIndex("s_s_2"), 2);
      printf("Length of 2nd straight: %6.3f m\n", Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);
      loc1 = Elem_GetPos(ElemIndex("s_s_3"), 1);
      loc2 = Elem_GetPos(ElemIndex("s_s_3"), 2);
      printf("Length of 3rd straight: %6.3f m\n", Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);

      prt_b2(b2_prms, b2);

      prtmfile("flat_file.fit");
      prt_lat(loc[0], loc[4], "linlat1.out", globval.bpm, true);
      prt_lat(loc[0]+1, loc[4]-1, "linlat.out", globval.bpm, true, 10);

      chi2_prt = min(chi2, chi2_prt);
    }

    chi2_ref = min(chi2, chi2_ref);
  }

  return chi2;
}


void fit_match(param_type &b2_prms)
{
  int          n_b2, i, j, iter;
  double       *b2, **xi, fret;
  ss_vect<tps> Ascr;

  n_b2 = b2_prms.n_prm;

  b2 = dvector(1, n_b2); xi = dmatrix(1, n_b2, 1, n_b2);

  // Upstream of QF03.
  // loc[0] = Elem_GetPos(ElemIndex("qf03"), 3);
  loc[0] = Elem_GetPos(ElemIndex("qf031"), 1);
  // Upstream of QD04.
  // loc[0] = Elem_GetPos(ElemIndex("qd041"), 1);
  // Upstream of 20 degree dipole.
  // loc[0] = Elem_GetPos(ElemIndex("sb"), 7);

  // Downstream of 10 degree dipole.
  loc[1] = Elem_GetPos(ElemIndex("b10"), 1);
  // Center of 1st straight.
  loc[2] = Elem_GetPos(ElemIndex("ef2"), 4);
  // Center of 2nd straight.
  loc[3] = Elem_GetPos(ElemIndex("ef2"), 16);

  // Upstream of EQ01.
  // loc[4] = Elem_GetPos(ElemIndex("eq01"), 1);
  // End of 1st straight.
  // loc[4] = Elem_GetPos(ElemIndex("eq05"), 1);
  // Downstream of 20 degree dipole.
  // loc[4] = Elem_GetPos(ElemIndex("b20"), 5);
  // Downstream of QF03.
  // loc[4] = Elem_GetPos(ElemIndex("qf03"), 6);
  loc[4] = Elem_GetPos(ElemIndex("qf031"), 4);

  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0], loc[4], Ascr, false, false, 0e0);

  prt_lat(loc[0], loc[4], "linlat1.out", globval.bpm, true);
  prt_lat(loc[0]+1, loc[4]-1, "linlat.out", globval.bpm, true, 10);

  prt_lin_opt(loc);

  b2_prms.ini_prm(b2);

  // Set initial directions (unit vectors).
  for (i = 1; i <= n_b2; i++)
    for (j = 1; j <= n_b2; j++)
      xi[i][j] = (i == j)? 1e0 : 0e0;

  n = 0;
  dpowell(b2, xi, n_b2, 1e-16, &iter, &fret, f_match);

  b2_prms.set_prm(b2);
  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc[0], loc[4], Ascr, false, false, 0e0);

  prt_lat(loc[0], loc[4], "linlat1.out", globval.bpm, true);
  prt_lat(loc[0]+1, loc[4]-1, "linlat.out", globval.bpm, true, 10);

  free_dvector(b2, 1, n_b2); free_dmatrix(xi, 1, n_b2, 1, n_b2);
}


void chk_straights()
{
  int loc1, loc2;
  
  loc1 = Elem_GetPos(ElemIndex("s_s_1"), 1);
  loc2 = Elem_GetPos(ElemIndex("s_s_1"), 2);
  printf("\nLength of 1st straight: %6.3f m\n",
	 Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);
  loc1 = Elem_GetPos(ElemIndex("s_s_2"), 1);
  loc2 = Elem_GetPos(ElemIndex("s_s_2"), 2);
  printf("Length of 2nd straight: %6.3f m\n",
	 Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);
  loc1 = Elem_GetPos(ElemIndex("s_s_3"), 1);
  loc2 = Elem_GetPos(ElemIndex("s_s_3"), 2);
  printf("Length of 3rd straight: %6.3f m\n",
	 Lattice.Cell[loc2].S-Lattice.Cell[loc1].S);
}


int main(int argc, char *argv[])
{
  int          loc0, loc1;
  ss_vect<tps> Ascr;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  no_sxt();

  if (false) {
    chk_straights();
    exit(0);
  }

  if (!qf031)
    loc0 = Elem_GetPos(ElemIndex("qf03"), 3);
  else
    loc0 = Elem_GetPos(ElemIndex("qf031"), 1);
  // loc0 = Elem_GetPos(ElemIndex("qd04"), 7);
  // loc0 = Elem_GetPos(ElemIndex("qd041"), 1);
  // loc1 = Elem_GetPos(ElemIndex("eq05"), 1);
  loc1 = Elem_GetPos(ElemIndex("b20"), 5);
  Ascr = get_A(ic[0], ic[1], ic[2], ic[3]);
  Cell_Twiss(loc0, loc1, Ascr, false, false, 0e0);

  prtmfile("flat_file.dat");
  prt_lat(loc0, loc1, "linlat1.out", globval.bpm, true);
  prt_lat(loc0+1, loc1-1, "linlat.out", globval.bpm, true, 10);
  // exit(0);

  b2_prms.add_prm("qf031", 2, -4.2, 4.2, 1.0);
  b2_prms.add_prm("qd041", 2, -4.2, 4.2, 1.0);

  b2_prms.add_prm("q01",   2, -6.0, 6.0, 1.0);
  // b2_prms.add_prm("q02",   2, -4.2, 4.2, 1.0);
  b2_prms.add_prm("q03",   2, -6.0, 6.0, 1.0);

  b2_prms.add_prm("eq01",  2, -4.2, 4.2, 1.0);
  b2_prms.add_prm("eq02",  2, -4.2, 4.2, 1.0);

  b2_prms.add_prm("eq03",  2, -4.2, 4.2, 1.0);
  b2_prms.add_prm("eq04",  2, -4.2, 4.2, 1.0);
  b2_prms.add_prm("eq05",  2, -4.2, 4.2, 1.0);

  n_strength = 9;

  if (false) {
    // b2_prms.add_prm("q01",  -2,  0.0,   0.05, 1.0);
    // b2_prms.add_prm("q02",  -2,  0.0,   0.05, 1.0);
    // b2_prms.add_prm("q03",  -2,  0.0,   0.05, 1.0);

    // b2_prms.add_prm("eq01", -2,  0.0,   0.05, 1.0);
    // b2_prms.add_prm("eq02", -2,  0.0,   0.05, 1.0);

    // b2_prms.add_prm("eq03", -2, -0.05,  0.05, 1.0);
    // b2_prms.add_prm("eq04", -2,  0.0,   0.05, 1.0);
    // b2_prms.add_prm("eq05", -2,  0.0,   0.05, 1.0);

    // b2_prms.add_prm("b10",  -2, -0.02,  0.02, 1.0);
  }

  // U561 + U562: 2.14.
  // b2_prms.add_prm("u561", -1, 2.14, 2.14, 1.0);

  b2_prms.bn_tol = 1e-6; b2_prms.step = 1.0;

  fit_match(b2_prms);
}
