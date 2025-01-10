#define NO 6

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


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


void prt_x0(const double delta)
{
  long int        lastpos;
  int             k;
  ss_vect<double> ps;
  FILE            *outf;

  const std::string file_name = "x0.out";

  outf = file_write(file_name.c_str());

  ps.zero(); ps[delta_] = delta;
  Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
  for (k = 0; k < lastpos; k++) {
    fprintf(outf, "%4d %15s %9.5f %4.1f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
	   k, Cell[k].Elem.PName, Cell[k].S, get_code(Cell[k]),
	   Cell[k].BeamPos[x_], Cell[k].BeamPos[px_],
	   Cell[k].BeamPos[py_], Cell[k].BeamPos[py_],
	   Cell[k].BeamPos[delta_], Cell[k].BeamPos[ct_]);
  }
  fclose(outf);
}


void scan_delta(const int n, const double delta)
{
  long int        lastpos;
  int             k;
  double          d;
  ss_vect<double> ps;
  std::ofstream        outf;

  const std::string file_name = "scan_delta.out";

  file_wr(outf, file_name.c_str());

  for (k = 0; k < n; k++) {
    d = (double)k/double(n-1)*delta;
    ps.zero(); ps[delta_] = d; Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    outf << std::scientific << std::setprecision(6)
	 << std::setw(14) << d << std::setw(14) << ps << std::endl;
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
    std::cout << std::endl << std::scientific << std::setprecision(3)
	 << std::setw(5) << n_iter << std::setw(11) << f << " T566 = " << T566
	 << " x_max = " << x_max << " px_max = " << px_max<< std::endl;
    for (k = 1; k <= n_prm; k++)
      std::cout << std::scientific << std::setprecision(5) << std::setw(13) << prms[k];
    std::cout << std::endl;

    f_ref = f;
  }
  
  return f;
}


void opt_nl_disp(void)
{
  int    i, j, iter;
  double *prms, **xi, fret, an;

  const double ftol  = 1e-10;

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

  std::cout << std::endl;
  n_iter = 0; n_prt = 1; f_ref = 1e30;
  dpowell(prms, xi, n_prm, ftol, &iter, &fret, f_opt);
  n_prt = 1; f_opt(prms);

  free_dvector(prms, 1, n_prm); free_dmatrix(xi, 1, n_prm, 1, n_prm);
}


template<typename T>
ss_vect<T> x_px2x_xp(const ss_vect<T> map)
{
  // Transform from [x, px, y, py, ct, delta] to [x, x', y, y', ct, delta]
  int        k;
  T          p_s;
  ss_vect<T> Id, Id_scl, map1;

  Id.identity(); Id_scl.identity(); map1 = map;
  p_s = get_p_s(Id);
  for (k = 0; k < 2; k++)
    Id_scl[2*k+1] *= p_s;
  map1 = map*Id_scl;
  for (k = 0; k < 2; k++)
    map1[2*k+1] /= p_s;

  return map1;
}


void analyze_nl_dyn(const double delta)
{
  int             k;
  tps             h;
  ss_vect<double> ps;
  ss_vect<tps>    Id, Id_scl, lin_map;

  Id.identity();
  Id_scl.identity(); Id_scl[delta_] *= delta;

  danot_(no_tps-1); get_map(false); danot_(no_tps);
  danot_(1); lin_map = map; danot_(no_tps);
  h = LieFact(map*Inv(lin_map));

  prt_lin_map(3, map);
  std::cout << std::endl << std::scientific << std::setprecision(6)
       << "R_56:    " << std::setw(14) << h_ijklm(map[ct_], 0, 0, 0, 0, 1) << std::endl;
  std::cout << std::scientific << std::setprecision(6)
       << "T_566:   " << std::setw(14) << h_ijklm(map[ct_], 0, 0, 0, 0, 2) << std::endl;

  std::cout << std::endl << std::fixed << std::setprecision(1)
       << "x0(" << -1e2*delta << "%) [m]:" << std::endl;
  for (k = 1; k <= no_tps-1; k++) {
    std::cout << std::scientific << std::setprecision(6)
	 << "  x0(delta^" << k << ") = "
	 << std::setw(14) << h_ijklm(map[x_]*Id_scl, 0, 0, 0, 0, k) << std::endl;
  }

  std::cout << std::endl << std::fixed << std::setprecision(1)
       << "px0(" << -1e2*delta << "%) [rad]:" << std::endl;
  for (k = 1; k <= no_tps-1; k++) {
    std::cout << std::scientific << std::setprecision(6)
	 << "px0(delta^" << k << ") = "
	 << std::setw(14) << h_ijklm(map[px_]*Id_scl, 0, 0, 0, 0, k) << std::endl;
  }
  // cout << scientific << setprecision(6) << setw(13) << h*Id_scl << endl;
  danot_(no_tps-1);
  // cout << scientific << setprecision(6) << setw(13) << map << endl;
  // cout << scientific << setprecision(6) << setw(13) << x_px2x_xp(map) << endl;
}


void get_optics()
{
  int     k;
  Vector2 alpha, beta, eta, etap;

  for (k = 0; k <= 1; k++) {
    alpha[k] = 0e0; beta[k] = 0e0; eta[k] = 0e0; etap[k] = 0e0;
  }

  // MATCHARC3.
  alpha[X_] = -2.384; alpha[Y_] = -0.7207;
  beta[X_] = 18.87; beta[Y_] = 1.931;

  // LINAC1.
  // beta[X_] = 4e0; beta[Y_] = 4e0;
  // ARC3.
  // beta[X_] = 5.9942; beta[Y_] = 1.8373;
 
  ttwiss(alpha, beta, eta, etap, 0e0);

  printf("\nalpha = [%13.6e, %13.6e]\n",
	 Cell[globval.Cell_nLoc].Alpha[X_],
	 Cell[globval.Cell_nLoc].Alpha[Y_]);
  printf("beta  = [%13.6e, %13.6e]\n",
	 Cell[globval.Cell_nLoc].Beta[X_],
	 Cell[globval.Cell_nLoc].Beta[Y_]);
  printf("nu    = [%13.6e, %13.6e]\n",
	 Cell[globval.Cell_nLoc].Nu[X_],
	 Cell[globval.Cell_nLoc].Nu[Y_]);
  printf("eta   = [%13.6e, %13.6e]\n",
	 Cell[globval.Cell_nLoc].Eta[X_],
	 Cell[globval.Cell_nLoc].Eta[Y_]);
  printf("etap  = [%13.6e, %13.6e]\n",
	 Cell[globval.Cell_nLoc].Etap[X_],
	 Cell[globval.Cell_nLoc].Etap[Y_]);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
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

  if (false) no_sxt();

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


int main(int argc, char *argv[])
{

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false;  globval.bpm         = 0;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  std::string home_dir = "";

  if (false)
    Read_Lattice((home_dir + argv[1]).c_str());
  else
    rdmfile(argv[1]);

  if (false) chk_bend();

  // opt_nl_disp();

  prt_x0(-5e-2);
  scan_delta(20, -5e-2);

  analyze_nl_dyn(5e-2);

  prtmfile("flat_file.dat");
}
