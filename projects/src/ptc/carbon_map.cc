#define NO 3

#include "tracy_lib.h"

int no_tps   = NO,
    ndpt_tps = 5;


const bool NSLS_II = false;


ss_vect<tps> get_fix_point(const int Fnum)
{
  long int     lastpos, loc;
  int          j;
  FieldMapType *FM;
  
  const int  n_iter = 3;

  
  loc = Elem_GetPos(Fnum, 1); FM = Cell[loc].Elem.FM;

  FM->Ld = 0.0; FM->L1 = 0.0;
  if (!NSLS_II) {
    // DIAMOND
    FM->cut = 0; FM->x0 = -30e-3;
  } else {
    // NSLS-II
    FM->cut = 50; FM->x0 = 55e-3;
  }

  trace = true;
  for (j = 1; j <= n_iter; j++) {
    map.identity(); Cell_Pass(loc, loc, map, lastpos);

    FM->phi -= map[px_].cst();
    FM->Lm = map[ct_].cst();
    FM->Ld += 2.0*map[x_].cst()/FM->phi;
    FM->L1 += Cell[loc].Elem.PL - FM->Lm;
  }

  map.identity(); Cell_Pass(loc, loc, map, lastpos);

  cout << endl;
  cout << scientific << setprecision(3)
       << "get_fix_pont:" << setw(11) << map.cst();
  cout << fixed << setprecision(5)
       << "phi [deg] = " << setw(7) << FM->phi*180.0/M_PI
       << ", L [m] = " << setw(7) << Cell[loc].Elem.PL << endl;
  cout << fixed << setprecision(5)
       << "Lr [m] = " << setw(7) << FM->Lr
       << ", Lm [m] = " << setw(7) << FM->Lm
       << ", Ld [m] = " << setw(7) << FM->Ld
       << ", L1 [m] = " << setw(7) << FM->L1 << endl;

  for (j = 1; j <= GetnKid(Fnum); j++) {
    loc = Elem_GetPos(Fnum, j);
    Cell[loc].Elem.FM->cut = FM->cut;
    Cell[loc].Elem.FM->x0 = FM->x0;
    Cell[loc].Elem.FM->phi = FM->phi;
    Cell[loc].Elem.FM->Lr = FM->Lr; Cell[loc].Elem.FM->Lm = FM->Lm;
    Cell[loc].Elem.FM->Ld = FM->Ld; Cell[loc].Elem.FM->L1 = FM->L1;
  }

  return map;
}


ss_vect<tps> rd_egt(const string &file_name)
{
  // The coordinates are: [x, x', y, y', ct, delta].
  string       line, symbol;
  int          i, j, k;
  double       val;
  ss_vect<tps> map;
  stringstream str;
  ifstream     inf;

  const bool prt = false;

  inf.open(file_name.c_str());

  // 0th order.
  if (prt) printf("\n");
  map.zero();
  getline(inf, line);
  str.clear(); str.str(""); str << line;
  str >> symbol;
  if (prt) printf("%s   ", symbol.c_str());
  for (k = 0; k < 6; k++) {
    str >> val;
    map[k] += val;
    if (prt) printf("%11.3e", val);
  }
  if (prt) printf("\n");

  // 1st order.
  for (j = 0; j < 6; j++) {
    getline(inf, line);
    str.clear(); str.str(""); str << line;
    str >> symbol;
    if (prt) printf("%s  ", symbol.c_str());
    for (k = 0; k < 6; k++) {
      str >> val;
      map[j] += val*tps(0e0, k+1);
      if (prt) printf("%11.3e", val);
    }
    if (prt) printf("\n");
  }

  // 2d order.
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      getline(inf, line);
      str.clear(); str.str(""); str << line;
      str >> symbol;
      if (prt) printf("%s ", symbol.c_str());
      for (k = 0; k <= j; k++) {
	str >> val;
	if (prt) printf("%11.3e", val);
	map[i] += val*tps(0e0, j+1)*tps(0e0, k+1);
	if (j != k) map[i] += val*tps(0e0, k+1)*tps(0e0, j+1);
      }
      if (prt) printf("\n");
    }
  }

  inf.close();

  return map;
}


template<typename T>
inline T get_p_z(const ss_vect<T> &x)
{
  T  p_s, p_s2;

  if (true)
    // Small angle axproximation.
    p_s = 1e0+x[delta_];
  else {
    p_s2 = sqr(1e0+x[delta_]) - sqr(x[px_]) - sqr(x[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
//      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}


void bend_fringe(const double L, const double phi, const double e1,
		 ss_vect<tps> &ps)
{
  double       rho, coeff;
  tps          u, pz1, pz2, pz3;
  ss_vect<tps> ps1;

  rho = L*180e0/(phi*M_PI);

  p_rot(e1, ps); 

  coeff = -1e0/(2e0*rho);
  ps1 = ps;
  pz1 = get_p_z(ps); pz2 = sqr(pz1); pz3 = cube(pz1);
  u = 1e0 + 4e0*coeff*ps1[px_]*ps1[y_]*ps1[py_]/pz3;
  ps[y_] = 2e0*ps1[y_]/(1e0+sqrt(u));
  ps[x_] = ps1[x_] - coeff*sqr(ps[y_])*(pz2+sqr(ps1[px_]))/pz3;
  ps[py_] = ps1[py_] + 2e0*coeff*ps1[px_]*ps[y_]/pz1;
  ps[ct_] = ps1[ct_] - coeff*ps1[px_]*sqr(ps[y_])*(1e0+ps1[delta_])/pz3;
}


ss_vect<tps> get_sbend(const double L, const double phi)
{
  double       rho;
  ss_vect<tps> Id, map;

  Id.identity();

  rho = L*180e0/(phi*M_PI);

  map.identity();
  map[x_]  +=
    rho*sin(phi*M_PI/180e0)*Id[px_] + rho*(1e0-cos(phi*M_PI/180e0))*Id[delta_];
  map[px_] += 2e0*tan(phi/2e0*M_PI/180e0)*Id[delta_];
  map[y_]  += -L*tan(phi*M_PI/(2e0*180e0))/rho*Id[y_] + L*Id[py_];
  map[py_] +=
    tan(phi*M_PI/(2e0*180e0))*(L*tan(phi*M_PI/(2e0*180e0))-2e0*rho)
    /sqr(rho)*Id[y_] - L*tan(phi*M_PI/(2e0*180e0))/rho*Id[py_];

  return map;
}


ss_vect<tps> get_map(const double T[5][5][5])
{
  int          i, j, k;
  ss_vect<tps> Id, map;

  Id.identity(); map.identity();
  for (i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      for (k = 1; k <= 4; k++) {
	map[i-1] += T[i][j][k]*Id[j-1]*Id[k-1];
	if (j != k) map[i-1] += T[i][j][k]*Id[k-1]*Id[j-1];
      }

  return map;
}


ss_vect<tps> get_e1_map(const double L, const double phi, const double e1,
			const bool entrance)
{
  int    i, j, k;
  double rho, sgn, T[5][5][5];

  for (i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++)
      for (k = 1; k <= 4; k++)
	T[i][j][k] = 0e0;

  rho = L*180e0/(phi*M_PI);

  sgn = (entrance)? 1e0 : -1e0;

  T[1][1][1] = -sgn*sqr(tan(e1*M_PI/180e0))/(2e0*rho);
  T[1][3][3] = sgn*1e0/(sqr(cos(e1*M_PI/180e0))*2e0*rho);
  T[2][1][1] = -sgn*cube(tan(e1*M_PI/180e0))/(2e0*sqr(rho));
  T[2][1][2] = -T[1][1][1];
  T[2][3][3] = -T[2][1][1];
  T[2][3][4] = T[1][1][1];
  T[3][1][3] = -T[1][1][1];
  T[4][1][3] = sgn*cube(tan(e1*M_PI/180e0))/(2e0*sqr(rho));
  T[4][1][4] = T[1][1][1];
  T[4][2][3] = -sgn*1/(sqr(cos(e1*M_PI/180e0))*2e0*rho);

  return get_map(T);
}


tps get_e1_h(const double L, const double phi, const double e1,
	     const bool entrance)
{
  double       h, sgn;
  ss_vect<tps> Id;

  Id.identity();

  h = phi*M_PI/(L*180e0);

  sgn = (entrance)? 1e0 : -1e0;

  return
    sgn*(-sqr(h)*cube(tan(e1*M_PI/180e0))/6e0*cube(Id[x_])
	 + h*sqr(tan(e1*M_PI/180e0))/2e0*sqr(Id[x_])*Id[px_]
	 + sqr(h)*cube(tan(e1*M_PI/180e0))/2e0*Id[x_]*sqr(Id[y_])
	 - h*sqr(tan(e1*M_PI/180e0))*Id[x_]*Id[y_]*Id[py_]
	 - h/(2e0*sqr(cos(e1*M_PI/180e0)))*Id[px_]*sqr(Id[y_]));
}


ss_vect<tps> dip_Brown_Helm(const double L, const double phi,
			    const double e1, const double e2)
{
  tps          h, h_e1, h_e2;
  ss_vect<tps> map, map_e1, map_e2, R;

  map = get_sbend(L, phi);
  map_e1 = get_e1_map(L, phi, phi/2e0, true);
  map_e2 = get_e1_map(L, phi, phi/2e0, false);
  h_e1 = get_e1_h(L, phi, phi/2e0, true);
  h_e2 = get_e1_h(L, phi, phi/2e0, false);

  prt_lin_map(3, map);

  daeps_(1e-10); h_e1 = 1e0*h_e1;
  daeps_(1e-10); h_e2 = 1e0*h_e2;
  cout << scientific << setprecision(3) << setw(11) << h_e1 << "\n";
  cout << scientific << setprecision(3) << setw(11) << h_e2 << "\n";

  h = LieFact_DF(map_e1, R); daeps_(1e-10); h_e1 = 1e0*h_e1;
  cout << scientific << setprecision(3) << setw(11) << h << "\n";
  h = LieFact_DF(map_e2, R); daeps_(1e-10); h = 1e0*h;
  cout << scientific << setprecision(3) << setw(11) << h << "\n";

  // return map_e2*map*map_e1;
  return map_e1*map*map_e2;
}


int main(int argc, char *argv[])
{
  bool         tweak;
  long int     lastpos;
  double       dx;
  Matrix       M;
  tps          h;
  ss_vect<tps> map, R, Id;
  ofstream     outf;

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;

  // 1: DIAMOND, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 1; sympl = false;

  // disable from TPSALib and LieLib log messages
  idprset(-1);

  Read_Lattice(argv[1]);

  // no_sxt();

//  Ring_GetTwiss(true, 0.0); printglob();

  if (true) {
    trace = false;

    globval.H_exact    = true;

    map.identity();
    // Tweak to remain within field map range at entrance.
    tweak = false;
    if (tweak) {
      dx = -1.4e-3; map[x_] += dx;
    }
    Cell_Pass(Elem_GetPos(ElemIndex("bb"), 1),
	      Elem_GetPos(ElemIndex("bb"), 1), map, lastpos);
    if (tweak) map[x_] -= dx;

    h = LieFact_DF(map, R);

    getlinmat(6, map, M);
    prt_lin_map(3, map);

    cout << "\n" << scientific << setprecision(3)
	 << setw(11) << map.cst() << "\n";
    cout << scientific << setprecision(3)
	  << endl << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;

    Id.identity(); Id[delta_] = 0e0; Id[ct_] = 0e0;

    // file_wr(outf, "h.dat");
    // Remove numeric noise.
    daeps_(1e-10); R = 1e0*R; h = 1e0*h; map = 1e0*map;
    prt_lin_map(3, R);
    cout << scientific << setprecision(3)
	 << setw(11) << h << "\n";
    // outf.close();
    exit(0);
  }

  if (false) {
    map = rd_egt("dip_egt.dat");
    h = LieFact_DF(map, R);

    daeps_(1e-10); R = 1e0*R; h = 1e0*h;
    prt_lin_map(3, R);
    cout << scientific << setprecision(3) << setw(11) << h << "\n";

    exit(0);
  }

  if (false) {
    map = get_fix_point(ElemIndex("bb"));

    prt_lin_map(3, map);

    getlinmat(6, map, M);
    cout << endl;
    cout << scientific << setprecision(3)
	 << "1-Det: " << setw(9) << 1-DetMat(6, M) << endl;

    exit(0);
  }

  if (true) {
    map = dip_Brown_Helm(0.933, 7.5, 7.5/2.0, 7.5/2.0);

    h = LieFact_DF(map, R);

    daeps_(1e-10); R = 1e0*R; h = 1e0*h;
    prt_lin_map(3, R);
    cout << scientific << setprecision(3) << setw(11) << h << "\n";

    exit(0);
  }

  if (false) {
    map.identity();
    bend_fringe(0.933, 7.5, 7.5/2.0, map);

    h = LieFact_DF(map, R);
    // Remove numeric noise.
    daeps_(1e-10); R = 1e0*R; h = 1e0*h;
    prt_lin_map(3, R);
    cout << scientific << setprecision(3)
	 << setw(11) << h << "\n";

    exit(0);
  }

  if (true) {
    trace = true;

    Ring_GetTwiss(true, 0.0); printglob();

    prt_lat("linlat.out", globval.bpm, true);  
  }
}
