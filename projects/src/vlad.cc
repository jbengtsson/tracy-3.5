#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


int    n_prm, n_iter, n_prt, Fnums[25];
double f_ref;


int get_ind(const int k)
{
  int index[] = {x_, px_, y_, py_, ct_, delta_};
  return index[k];
}


void chk_bend()
{
  int    k;
  double phi;

  phi = 0e0;
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if ((Cell[k].Elem.Pkind == Mpole) &&
	(Cell[k].Elem.M->n_design == Dip)) {
      phi += Cell[k].Elem.PL*Cell[k].Elem.M->Pirho;
    }
  }
  phi = 180e0*phi/M_PI;
  printf("\nphi = %8.6f\n", phi);
}


void scan_delta(const int n, const double delta)
{
  long int        lastpos;
  int             k;
  double          d;
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "scan_delta.out";

  file_wr(outf, file_name.c_str());

  for (k = 0; k < n; k++) {
    d = (double)k/double(n-1)*delta;
    ps.zero(); ps[delta_] = d; Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << scientific << setprecision(6)
	 << setw(14) << d << setw(14) << ps << endl;
  }
  outf.close();
}


template<typename T>
T trapz(const int n, const double h, const T y[])
{
  int k;
  T   intgrl;

  intgrl = 0e0;
  for (k = 0; k < n-1; k++) {
    intgrl += (y[k+1]+y[k])/2e0*h;
  }

 return intgrl;
}


double f_opt(double prms[])
{
  const int n = 25;

  long int        lastpos;
  int             k;
  double          h, x_max, px_max;
  double          deltas[n], x_delta[n], px_delta[n], ct_delta[n];
  double          x2d_intgrl, px2d_intgrl, sigma, f, T566;
  psVector          b;
  ss_vect<double> ps;

  const int    n_pol = 3;
  const double delta = -5e-2, scl_T566 = 1e-11, T566_ref = -1.125e0;

  n_iter++;

  for (k = 1; k <= n_prm; k++)
    set_bn_design_fam(Fnums[k-1], Sext, prms[k], 0e0);

  x_max = 0e0; px_max = 0e0; h = delta/(n-1);
  for (k = 0; k < n; k++) {
    deltas[k] = k*h;
    ps.zero(); ps[delta_] = deltas[k];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    x_delta[k] = sqr(ps[x_]); px_delta[k] = sqr(ps[px_]); ct_delta[k] = ps[ct_];
    x_max = max(fabs(ps[x_]), x_max); px_max = max(fabs(ps[px_]), px_max);
  }

  x2d_intgrl = trapz(n, h, x_delta); px2d_intgrl = trapz(n, h, px_delta);

  pol_fit(n, deltas, ct_delta, n_pol, b, sigma, false); T566 = b[2];

  f = scl_T566*sqr(T566-T566_ref) + sqr(x2d_intgrl) + sqr(px2d_intgrl);

  if (f < f_ref) {
    prtmfile("flat_file.dat");
    pol_fit(n, deltas, ct_delta, n_pol, b, sigma, true);
    scan_delta(n, delta);
    cout << endl << scientific << setprecision(3)
	 << setw(5) << n_iter << setw(11) << f << " T566 = " << T566
	 << " x_max = " << x_max << " px_max = " << px_max<< endl;
    for (k = 1; k <= n_prm; k++)
      cout << scientific << setprecision(5) << setw(13) << prms[k];
    cout << endl;

    f_ref = f;
  }
  
  return f;
}


void opt_nl_disp(void)
{
  int    i, j, iter;
  double *prms, **xi, fret, an;

  const double ftol  = 1e-5;

  n_prm = 4;

  prms = dvector(1, n_prm); xi = dmatrix(1, n_prm, 1, n_prm);

  i = 0;
  Fnums[i++] = ElemIndex("a3s1"); Fnums[i++] = ElemIndex("a3s2");
  Fnums[i++] = ElemIndex("a3s3"); Fnums[i++] = ElemIndex("a3s4");

  for (i = 0; i < n_prm; i++) {
    get_bn_design_elem(Fnums[i], 1, Sext, prms[i+1], an);
    // prms[i+1] = 81.6515/2.0;
  }

  for (i = 1; i <= n_prm; i++) {
    for (j = 1; j <= n_prm; j++)
      if (i == j)
        xi[i][j] = 0.1e0;
      else
        xi[i][j] = 0e0;
  }

  cout << endl;
  n_iter = 0; n_prt = 1; f_ref = 1e30;
  dpowell(prms, xi, n_prm, ftol, &iter, &fret, f_opt);
  n_prt = 1; f_opt(prms);

  free_dvector(prms, 1, n_prm); free_dmatrix(xi, 1, n_prm, 1, n_prm);
}


void get_nu(const double delta, double nu[])
{
  long int     lastpos;
  ss_vect<tps> A;

  A.zero(); putlinmat(4, Cell[0].A, A); A[delta_] += delta;
  Cell_Pass(0, globval.Cell_nLoc, A, lastpos);
  A = get_A_CS(2, A, nu);
}


void get_ksi(double ksi[])
{
  int    k;
  double nu0[2], nu1[2];

  get_nu(-globval.dPcommon, nu0); get_nu(globval.dPcommon, nu1);
  for (k = 0; k < 2; k++)
    ksi[k] = (nu1[k]-nu0[k])/(2e0*globval.dPcommon);
}


void get_optics(const string &line)
{
  int     k;
  double  nu[2], ksi[2];
  Vector2 alpha, beta, eta, etap;

  cout << endl << line << endl;

  for (k = 0; k < 2; k++) {
    alpha[k] = 0e0; beta[k] = 0e0; eta[k] = 0e0; etap[k] = 0e0;
  }

  if (line.compare("arc3") == 0) {
    // beta[X_] = 5.9942; beta[Y_] = 1.8373;
    // beta[X_] = 4.28064; beta[Y_] = 1.21995; eta[X_] = 1.386;
    // beta[X_] = 3.10119; beta[Y_] = 1.44765; eta[X_] = 0.76;
    // beta[X_] = 0.35751; beta[Y_] = 5.64898; eta[X_] = 0.0698;
    if (false) {
      // Unit cell.
      // beta[X_] = 4.73529; beta[Y_] = 1.55552; eta[X_] = 0.31861;
      beta[X_] = 1.92895; beta[Y_] = 1.01853; eta[X_] = 0.6;
    } else {
      // JB_1.lat.
      // Fine tuned optics.
      // beta[X_] = 6.01442; beta[Y_] = 1.78576;
      if (false) {
	// Dispersion wave.
	beta[X_] = 8.70141; beta[Y_] = 2.00042;
      } else {
	// Matching section.
	// beta[X_]  = 1.52376; beta[Y_]  = 12.34581;
	// alpha[X_] = 3.35393; alpha[Y_] =  5.27033;
	beta[X_]  = 2.70756; beta[Y_]  = 10.31726;
	alpha[X_] = 1.09794; alpha[Y_] =  3.45210;
      }
      // beta[X_] = 5.9942; beta[Y_] = 1.8373;
      // alpha[X_] = 1.22596; alpha[Y_] = -2.53096;
      // beta[X_]  = 0.75827; beta[Y_]  = 6.86095;
      // beta[X_] = 2.67796; beta[Y_] = 1.19590;
      // beta[X_] = 3.61985; beta[Y_] = 4.55008;
      // 8-unit cell structure.
      // alpha[X_] = 0.92093; alpha[Y_] = -0.61452;
      // beta[X_]  = 1.00644; beta[Y_]  =  2.70480;
      // 4-unit cell structure.
      // alpha[X_] = 1.20744; alpha[Y_] = -0.18567;
      // beta[X_]  = 2.84388; beta[Y_]  =  3.90905;
    }
  } else if (line.compare("matcharc3") == 0) {
    alpha[X_] = -2.384; alpha[Y_] = -0.7207;
    beta[X_]  = 18.87;   beta[Y_] = 1.931;
  } else if (line.compare("chicane2arc3combiner2") == 0) {
    alpha[X_] = -0.25; alpha[Y_] = 0.00;
    beta[X_] =  2.50;  beta[Y_] = 5.00;
  } else if (line.compare("disp_wave") == 0) {
    alpha[X_] = 2.19110; alpha[Y_] = -2.77166;
    beta[X_]  = 1.44220;  beta[Y_] =  4.65242;
  } else if (line.compare("linac1") == 0) {
    beta[X_] = 4e0; beta[Y_] = 4e0;
  } else if (line.compare("arc1") == 0) {
    beta[X_] = 1.2413; beta[Y_] = 1.4246;
  }
 
  ttwiss(alpha, beta, eta, etap, 0e0);

  get_nu(0e0, nu); get_ksi(ksi);
  cout << fixed << setprecision(5)
       << "nu  = [" << setw(8) << nu[X_] << ", " << setw(8) << nu[Y_] << "]"
       << endl;
  cout << fixed << setprecision(5)
       << "ksi = [" << setw(8) << ksi[X_] << ", " << setw(8) << ksi[Y_] << "]"
       << endl;
}


double get_disp(const double eta_x0)
{
  int     k;
  double  d, eta_x;
  Vector2 alpha, beta, eta, etap;

  const int    max_iter = 25;
  const double eps = 1e-10;

  for (k = 0; k < 2; k++) {
    alpha[k] = 0e0; beta[k] = 0e0; eta[k] = 0e0; etap[k] = 0e0;
  }

  k = 0; eta_x = eta_x0;
  do {
    k++;
    beta[X_] = 1e0; beta[Y_] = 1e0; eta[X_] = eta_x;
    ttwiss(alpha, beta, eta, etap, 0e0);
    d = Cell[globval.Cell_nLoc].Eta[X_] - eta_x; eta_x += d/2e0;
    cout << fixed << setprecision(5)
	 << setw(3) << k << setw(9) << eta_x << endl;
  } while((fabs(d) > eps) && (k < max_iter));

  return eta_x;
}


void get_twoJ_phi(long int k, const ss_vect<double> &ps,
		  double twoJ[], double phi[])
{
  int             j;
  ss_vect<double> ps_Fl;
  ss_vect<tps>    A;

  const int n = 4;

  A.identity(); putlinmat(n, Cell[k].A, A); ps_Fl = (Inv(A)*ps).cst();
  for (j = 0; j < 2; j++) {
    twoJ[j] = sqr(ps_Fl[2*j]) + sqr(ps_Fl[2*j+1]);
    phi[j] = -atan2(ps_Fl[2*j+1], ps_Fl[2*j]);
    if (phi[j] < 0e0) phi[j] += 2e0*M_PI;
  }
}


void prt_twoJ(const int n, const double A[][2])
{
  long int        lastpos;
  int             j, k;
  double          twoJ[n][globval.Cell_nLoc+1][2], phi[2];
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "twoJ.out";

  cout << endl;
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(3)
	 << "[" << setw(11) << A[j][X_] << ", " << setw(11) << A[j][Y_] << "]"
	 << endl;
    ps.zero(); ps[x_] = A[j][X_]; ps[y_] =  A[j][Y_];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (lastpos != globval.Cell_nLoc)
      cout << endl << "Particle lost at: " << lastpos
	   << "(" << globval.Cell_nLoc << ")" << endl;
    for (k = 0; k <= globval.Cell_nLoc; k++)
      get_twoJ_phi(k, Cell[k].BeamPos, twoJ[j][k], phi);
  }

  file_wr(outf, file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    outf << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << fixed << setprecision(5) << setw(9) << Cell[k].S
	 << setprecision(1) << setw(5) << get_code(Cell[k]);
   for (j = 0; j < n; j++)
      outf << scientific << setprecision(5)
	   << setw(13) << twoJ[j][k][X_] << setw(13) << twoJ[j][k][Y_];
    outf << endl;
  }
  outf.close();
}


void prt_curly_H(const int n, const double delta[])
{
  long int        lastpos;
  int             j, k;
  double          curly_H[n][globval.Cell_nLoc+1][2], phi[2];
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "curly_H.out";

  cout << endl;
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(3) << setw(11) << 1e2*delta[j] << endl;
    ps.zero(); ps[delta_] = delta[j];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    for (k = 0; k <= globval.Cell_nLoc; k++)
      get_twoJ_phi(k, Cell[k].BeamPos, curly_H[j][k], phi);
  }

  file_wr(outf, file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    outf << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << fixed << setprecision(5) << setw(9) << Cell[k].S
	 << setprecision(1) << setw(5) << get_code(Cell[k]);
   for (j = 0; j < n; j++)
      outf << scientific << setprecision(5) << setw(13) << curly_H[j][k][X_];
    outf << endl;
  }
  outf.close();
}


void prt_eta(const int n, const double delta[])
{
  long int        lastpos;
  int             j, k;
  double          eta[n][globval.Cell_nLoc+1], etap[n][globval.Cell_nLoc+1];
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "eta.out";

  cout << endl;
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(2) << setw(10) << 1e2*delta[j] << endl;
    ps.zero(); ps[delta_] = delta[j];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);

    if (lastpos != globval.Cell_nLoc)
      cout << endl << "prt_eta: particle lost at " << lastpos
	   << "(" << globval.Cell_nLoc << ")" << endl;

    for (k = 0; k <= globval.Cell_nLoc; k++) {
      eta[j][k] = Cell[k].BeamPos[x_]; etap[j][k] = Cell[k].BeamPos[px_];
    }
  }

  file_wr(outf, file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    outf << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << fixed << setprecision(5) << setw(9) << Cell[k].S
	 << setprecision(1) << setw(5) << get_code(Cell[k]);
   for (j = 0; j < n; j++)
      outf << scientific << setprecision(5)
	   << setw(13) << eta[j][k] << setw(13) << etap[j][k];
    outf << endl;
  }
  outf.close();
}


void prt_dnu(const int n, const double delta[])
{
  int             j, k, l;
  double          dnu0[globval.Cell_nLoc+1][2];
  double          dnu[n][globval.Cell_nLoc+1][2];
  ofstream        outf;

  const string file_name = "dnu.out";

  ttwiss(Cell[0].Alpha, Cell[0].Beta, Cell[0].Eta, Cell[0].Etap, 0e0);

  for (k = 0; k <= globval.Cell_nLoc; k++)
    for (l = 0; l < 2; l++)
      dnu0[k][l] = Cell[k].Nu[l];

  cout << endl;
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(3) << setw(11) << 1e2*delta[j] << endl;
    ttwiss(Cell[0].Alpha, Cell[0].Beta, Cell[0].Eta, Cell[0].Etap, delta[j]);
    for (k = 0; k <= globval.Cell_nLoc; k++)
      for (l = 0; l < 2; l++)
	dnu[j][k][l] = Cell[k].Nu[l] - dnu0[k][l];
  }

  file_wr(outf, file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    outf << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << fixed << setprecision(5) << setw(9) << Cell[k].S
	 << setprecision(1) << setw(5) << get_code(Cell[k]);
   for (j = 0; j < n; j++)
      outf << scientific << setprecision(5)
	   << setw(13) << dnu[j][k][X_] << setw(13) << dnu[j][k][Y_];
    outf << endl;
  }
  outf.close();
}


void prt_trsvrs(const int n, const double A[][2])
{
  long int        lastpos;
  int             j, k, l, n1;
  double          twoJ[2], phi[2], phi0;
  ss_vect<double> ps, ps_Fl;
  ss_vect<tps>    A0;
  ofstream        outf;

  // End of lattice.
  const int loc = globval.Cell_nLoc, n_phi = 100;
  // End of ARC3.
  // const int loc = Elem_GetPos(ElemIndex("q0h"), 16), n_phi = 100;

  const string file_name = "twoJ_phi.out";

  file_wr(outf, file_name.c_str());

  cout << endl << setw(4) << loc << " " << setw(10) << Cell[loc].Elem.PName
       << fixed << setprecision(5) << setw(9) << Cell[loc].S << endl;

  cout << endl;
  n1 = 0; A0.identity(); putlinmat(4, Cell[0].A, A0);
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(3)
	 << "[" << setw(11) << A[j][X_] << ", " << setw(11) << A[j][Y_] << "]"
	 << endl;
    ps.zero(); ps[x_] = A[j][X_]; ps[y_] = A[j][Y_];
    get_twoJ_phi(0, ps, twoJ, phi);
    for (k = 0; k <= n_phi; k++) {
      n1++;
      phi0 = k*2e0*M_PI/n_phi;
      for (l = 0; l < 2; l++) {
	ps_Fl[2*l]   =  sqrt(twoJ[l])*cos(phi0);
	ps_Fl[2*l+1] = -sqrt(twoJ[l])*sin(phi0);
      }
      ps = (A0*ps_Fl).cst(); Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
      if (lastpos != globval.Cell_nLoc)
	cout << endl << "Particle lost at: " << lastpos
	     << "(" << globval.Cell_nLoc << ")" << endl;
      outf << scientific << setprecision(5)
      	   << setw(4) << n1 << setw(13) << Cell[loc].BeamPos << endl;
    }
    outf << endl;
  }

  outf.close();
}


void prt_long(const int n, const double A[][2])
{
  long int        lastpos;
  int             j, k, l, n1;
  double          twoJ[2], phi[2], phi0;
  ss_vect<double> ps, ps_Fl;
  ss_vect<tps>    A0;
  ofstream        outf;

  const int loc = globval.Cell_nLoc, n_phi = 100;

  const string file_name = "twoJ_phi.out";

  file_wr(outf, file_name.c_str());

  cout << endl << setw(4) << loc << " " << setw(10) << Cell[loc].Elem.PName
       << fixed << setprecision(5) << setw(9) << Cell[loc].S << endl;

  cout << endl;
  n1 = 0; A0.identity(); putlinmat(4, Cell[0].A, A0);
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(3)
	 << "[" << setw(11) << A[j][X_] << ", " << setw(11) << A[j][Y_] << "]"
	 << endl;
    ps.zero(); ps[x_] = A[j][X_]; ps[y_] = A[j][Y_];
    get_twoJ_phi(0, ps, twoJ, phi);
    for (k = 0; k <= n_phi; k++) {
      n1++;
      phi0 = k*2e0*M_PI/n_phi;
      for (l = 0; l < 2; l++) {
	ps_Fl[2*l]   =  sqrt(twoJ[l])*cos(phi0);
	ps_Fl[2*l+1] = -sqrt(twoJ[l])*sin(phi0);
      }
      ps = (A0*ps_Fl).cst();
      Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
      outf << scientific << setprecision(5)
      	   << setw(4) << n1 << setw(13) << Cell[loc].BeamPos << endl;
    }
    outf << endl;
  }

  outf.close();
}


void prt_ct(const int n, const double delta[])
{
  long int        lastpos;
  int             j, k;
  double          ct[n][globval.Cell_nLoc+1];
  ss_vect<double> ps;
  ofstream        outf;

  const string file_name = "ct.out";

  cout << endl;
  for (j = 0; j < n; j++) {
    cout << scientific << setprecision(2) << setw(10) << 1e2*delta[j] << endl;
    ps.zero(); ps[delta_] = delta[j];
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    for (k = 0; k <= globval.Cell_nLoc; k++) {
      ct[j][k] = Cell[k].BeamPos[ct_];
    }
  }

  file_wr(outf, file_name.c_str());
  for (k = 0; k <= globval.Cell_nLoc; k++) {
    outf << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << fixed << setprecision(5) << setw(9) << Cell[k].S
	 << setprecision(1) << setw(5) << get_code(Cell[k]);
   for (j = 0; j < n; j++)
      outf << scientific << setprecision(5)
	   << setw(13) << ct[j][k];
    outf << endl;
  }
  outf.close();
}


void check_cav_model()
{
  long int        lastpos;
  tps             p_s;
  ss_vect<double> ps;
  ss_vect<tps>    Id, map;

  // ARC3.
  const double p0 = 750e6;

  globval.Cavity_on = true;

  // globval.Energy contains p_0 [GeV/c].
  globval.Energy = 1e-9*p0;
  globval.gamma0 = sqrt(sqr(m_e)+sqr(1e9*globval.Energy))/m_e;
  globval.beta0  = sqrt(1e0-1e0/sqr(globval.gamma0));
  printf("\np0 = %12.5e, 1-beta0 = %12.5e, gamma0 = %12.5e\n",
	 1e9*globval.Energy, 1e0-globval.beta0, globval.gamma0);

  map.identity();
  Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  prt_lin_map(3, map);

  if (true)
    p_s =
      sqrt(1e0+2e0*map[delta_]/globval.beta0+sqr(map[delta_])-sqr(map[px_])
	   -sqr(map[py_]));
  else
    p_s = sqrt(sqr(1e9*globval.Energy/p0)-sqr(map[px_]) -sqr(map[py_]));

  Id.identity();
  Id[px_] /= p_s; Id[py_] /= p_s;
  // Transform from [p_x, p_y] to [x', y'].
  // map = map*Id;

  // Transform from ct to s.
  // map[ct_] *= globval.beta0*p_s;

  // prt_lin_map(3, map);
}


void get_lin_map(double beta[])
{
  long int     lastpos;
  int          k;
  ss_vect<tps> map;

  map.identity(); Cell_Pass(0, globval.Cell_nLoc, map, lastpos);
  prt_lin_map(3, map);

  // For a periodic cell.
  for (k = 0; k < 2; k++)
    beta[k] = sqrt(-map[2*k][2*k+1]/map[2*k+1][2*k]);
  cout << fixed << setprecision(5)
       << "beta = [" << beta[X_] << ", " << beta[Y_] << "]" << endl;
}


int main(int argc, char *argv[])
{
  double beta[2];

  // Problem with STL vector class.
  std::vector<std::string> dummy[2];
  dummy[X_].push_back("1"); dummy[Y_].push_back("2");

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false;  globval.bpm        = 0;

  string home_dir = "";

  const double A[][2] =
    {{1e-3, 1e-3}, {2.5e-3, 2.5e-3}, {5e-3, 5e-3}, {7.5e-3, 7.5e-3}};

  const double delta[] = {-5e-2, -2.5e-2, -1e-2, 1e-2, 2.5e-2, 5e-2};

  if (argc != 3) {
    cout << "*** bad command line." << endl;
    exit(1);
  }

  if (argv[2][0] == 'y')
    rdmfile(argv[1]);
  else {
    Read_Lattice((home_dir + argv[1]).c_str());
    if (false) no_sxt();
  }

  chk_bend();

  get_lin_map(beta);

  prtmfile("flat_file.dat");

  if (false) opt_nl_disp();

  if (false) scan_delta(20, -5e-2);

  switch (1) {
  case 1:
    get_optics("arc3");
    break;
  case 2:
    get_optics("chicane2arc3combiner2");
    break;
  case 3:
    get_optics("disp_wave");
    break;
  case 4:
    get_optics("linac1");
    break;
  case 5:
    get_optics("arc1");
    break;
  }

  // Messes up linear optics calculation.
  // if (false) get_disp(0.28);
  if (false) get_disp(0.6);

  if (false) {
    Ring_GetTwiss(true, 0.0); printglob();
  }

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (true) {
    prt_twoJ(4, A);
    prt_curly_H(6, delta);
    prt_eta(6, delta);
    prt_dnu(6, delta);
    prt_trsvrs(4, A);
    prt_long(4, A);
    prt_ct(6, delta);
  }
}
