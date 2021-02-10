#define NO 1

#include "tracy_lib.h"

int no_tps = NO;

const int    n_turn    = 2064;
const double delta_max = 2e-2;


void prt_kids(const int Fnum, ElemFamType elemf[])
{
  int k;

  for (k = 0; k < elemf[Fnum-1].nKid; k++)
    printf(" %2d", elemf[Fnum-1].KidList[k]);
  printf("\n");
}


void prt_fam(ElemFamType elemf[])
{
  int k;

  printf("\nFamilies:\n");
  for (k = 0; k < globval.Elem_nFam; k++)
    elemf[k].ElemF->print();
}


void prt_lat(ElemType *elems[])
{
  int k;

  printf("\nLattice:\n");
  for (k = 0; k <= globval.Cell_nLoc; k++)
    elems[k]->print();
}


void tst_lat(LatticeType &lat)
{
  int          k;
  ss_vect<tps> map;

  printf("\ntst_lat\n");
  map.identity();
  for (k = 0; k <= 10; k++)
    lat.elems[k]->Elem_Pass(map);
  prt_lin_map(3, map);
}


template<typename T>
void Elem_Pass_Lin(LatticeType &lat, ss_vect<T> ps)
{
  long int  k;
  MpoleType *Mp;

  for (k = 0; k <= globval.Cell_nLoc; k++) {
    if (lat.elems[k]->Pkind == Mpole) { 
      Mp = dynamic_cast<MpoleType*>(lat.elems[k]);
      if ((Mp->Pthick == thick) && (Mp->Porder <= Quad)) {
	ps =
	  is_double<ss_vect<T> >::ps(Mp->M_lin*ps);
	
	if (globval.emittance && !globval.Cavity_on
	    && (lat.elems[k]->PL != 0e0) && (Mp->Pirho != 0e0))
	  get_dI_eta_5(k, lat.elems);
      }
    } else
      lat.elems[k]->Elem_Pass(ps);
  }
}


void get_eps_x(LatticeType &lat, double &eps_x, double &sigma_delta,
	       double &U_0, double J[], double tau[], double I[],
	       const bool prt)
{
  int          k;
  ss_vect<tps> A;

  const double
    C_q_scl = 1e18*C_q/sqr(m_e),
    E_0     = 1e9*globval.Energy,
    C       = lat.elems[globval.Cell_nLoc]->S,
    T_0     = C/c0;

  globval.Cavity_on = false; globval.emittance = false;
  lat.Ring_GetTwiss(false, 0.0);
  A = putlinmat(6, globval.Ascr); A += globval.CODvect;
  globval.emittance = true;
  Elem_Pass_Lin(lat, A);
  lat.get_I(I, false);

  U_0 = 1e9*C_gamma*pow(globval.Energy, 4)*I[2]/(2e0*M_PI);
  eps_x = C_q_scl*sqr(globval.Energy)*I[5]/(I[2]-I[4]);
  sigma_delta = sqrt(C_q_scl*sqr(globval.Energy)*I[3]/(2e0*I[2]+I[4]));
  J[X_] = 1e0 - I[4]/I[2]; J[Z_] = 2e0 + I[4]/I[2]; J[Y_] = 4e0 - J[X_] - J[Z_];

  for (k = 0; k < 3; k++)
    tau[k] = 4e0*M_PI*T_0/(C_gamma*cube(1e-9*E_0)*J[k]*I[2]);

  if (prt) {
    printf("\n  I[1..5]:");
    for (k = 1; k <= 5; k++)
      printf(" %10.3e", I[k]);
    printf("\n");

    printf("\n  U_0   [keV]    = %5.1f\n", 1e-3*U_0);
    printf("  eps_x [nm.rad] = %6.4f\n", 1e9*eps_x);
    printf("  sigma_delta    = %9.3e\n", sigma_delta);
    printf("  J              = [%5.3f, %5.3f, %5.3f]\n", J[X_], J[Y_], J[Z_]);
    printf("  tau   [msec]   = [%e, %e, %e]\n",
	   1e3*tau[X_], 1e3*tau[Y_], 1e3*tau[Z_]);
  }
}


void get_lat(const char *file_name)
{
  LatticeType lat;
  double      eps_x, sigma_delta, U_0, J[3], tau[3], I[6];
  FILE        *inf, *outf;

  const string str = file_name;

  inf  = file_read((str + ".lat").c_str());
  outf = file_write((str + ".lax").c_str());

  Lattice_Read(inf, outf, lat.elemf);
  lat.Lat_Init();

  printf("\n%s\n", lat.elems[0]->PName);

  if (false) prt_fam(lat.elemf);
  if (false) prt_lat(lat.elems);

  lat.Ring_GetTwiss(true, 0e0); printglob(lat.elems[0]);
  if (globval.mat_meth)
    get_eps_x(lat, eps_x, sigma_delta, U_0, J, tau, I, true);

  if (false) tst_lat(lat);
}


int main(int argc, char *argv[])
{
  string file_name;
  double eps_x, sigma_delta, U_0, J[3], tau[3], I[6];

  globval.H_exact    = false; globval.quad_fringe = false;
  globval.Cavity_on  = false; globval.radiation   = false;
  globval.emittance  = false; globval.IBS         = false;
  globval.pathlength = false; globval.bpm         = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = !false;

  if (!true) {
    if (true)
      Read_Lattice(argv[1]);
    else
      rdmfile(argv[1]);
    lat.Ring_GetTwiss(true, 0e0); printglob(lat.elems[0]);
    // if (globval.mat_meth)
    //   lat.get_eps_x(eps_x, sigma_delta, U_0, J, tau, I, true);
  } else
    get_lat(argv[1]);

  exit(0);

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);
  prtmfile("flat_file.dat");
  exit(0);

  if (true) GetEmittance(ElemIndex("cav"), true);

  if (true) {
    globval.Cavity_on = true;
    get_dynap(delta_max, 25, n_turn, false);
  }
}
