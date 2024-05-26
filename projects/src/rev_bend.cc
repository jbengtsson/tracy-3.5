#include <assert.h>
#include <map>

#define NO 1

#include "tracy_lib.h"


int no_tps = NO;


double get_f_RF(const int Fnum)
{
  const int loc = Elem_GetPos(Fnum, 1);
  return Cell[loc].Elem.C->f_RF;
}


void set_f_RF(const int Fnum, const double f_RF)
{
  const int loc = Elem_GetPos(Fnum, 1);
  Cell[loc].Elem.C->f_RF = f_RF;
}


void scan_f_RF(const string &Fam_name, const double df_RF_max, const int n_step)
{
  const int
    Fnum = ElemIndex(Fam_name.c_str());
  const double
    f_RF_0    = get_f_RF(Fnum),
    f_RF_step = df_RF_max/n_step;

  double
    df_RF = 0e0;

  globval.Cavity_on = globval.radiation = globval.pathlength = true;

  printf("\n   f_Rf   eps_x    delta\n");
  printf("  [kHz]  [pm.rad]   [%%]\n");
  for (int k = 0; k <= n_step; k++) {
    GetEmittance(ElemIndex("cav"), true, false);
    printf("  %5.3f   %5.1f   %5.2f\n",
	   1e-3*df_RF, 1e12*globval.eps[X_], 1e2*globval.CODvect[delta_]);
    df_RF += f_RF_step;
    set_f_RF(Fnum, f_RF_0+df_RF);
  }

  prt_cod("cod.out", 0, true);
}


double get_elem_phi(const double Fnum, const double Knum)
{
  int    loc;
  double phi;

  loc = Elem_GetPos(Fnum, Knum);
  phi = Cell[loc].Elem.M->Pirho*Cell[loc].Elem.PL*180e0/M_PI;
  return phi;
}


void prt_beampos1(const string file_name, const long int lastpos)
{
  ofstream outf;

  file_wr(outf, file_name.c_str());

  outf << "#  k   Name                x           p_x           y"
       << "           p_y         delta         ct\n";
  outf << "#                         [m]         [rad]         [m]"
       << "         [rad]          []          [m]\n";
  for (int k = 0; k <= lastpos; k++)
    outf << scientific << setprecision(5)
	 << setw(4) << k << " " << setw(10) << Cell[k].Elem.PName
	 << setw(13) << Cell[k].BeamPos
	 << "\n";
  if (lastpos != globval.Cell_nLoc)
    outf << "  particle lost\n";

  outf.close();
}


double set_rev_bend(CellType &Cell)
{
  const double
    eps = 1e-5;

  using KeyType   = std::string;
  using ValueType = double;
  typedef std::map<KeyType, ValueType> DictType;

  string   key;
  double   phi, phi_0, rho, b_1;
  DictType dict;

  dict["b1_0           "] = 1.094181;
  dict["b1_1           "] = 0.151199;
  dict["b1_2           "] = 0.151101;
  dict["b1_3           "] = 0.101861;
  dict["b1_4           "] = 0.001569;
  dict["b1_5           "] = 0.000089;

  dict["b2u_6          "] = 0.001070;
  dict["b2u_5          "] = 0.050729;
  dict["b2u_4          "] = 0.074672;
  dict["b2u_3          "] = 0.076248;
  dict["b2u_2          "] = 0.114983;
  dict["b2u_1          "] = 0.152049;
  dict["b2_0           "] = 0.621695;
  dict["b2d_1          "] = 0.152220;
  dict["b2d_2          "] = 0.152122;
  dict["b2d_3          "] = 0.102549;
  dict["b2d_4          "] = 0.001579;
  dict["b2d_5          "] = 0.000090;

  key = Cell.Elem.PName;
  phi_0 = (dict.count(key) == 1)?
    dict[key] :
    0e0;
  phi = get_elem_phi(Cell.Fnum, Cell.Knum) - phi_0;
  rho = (fabs(phi) > eps)? Cell.Elem.PL*180e0/(phi*M_PI) : 0e0;
  b_1 = (fabs(rho) != 0e0)? 1e0/rho : 0e0;
  Cell.Elem.M->PBpar[Dip+HOMmax] = b_1;
  Cell.Elem.M->Pirho = phi_0*M_PI/(180e0*Cell.Elem.PL);
  Mpole_SetPB(Cell.Fnum, Cell.Knum, Dip);

  printf("  %10s %8.5f %10.3e %10.3e\n", key.c_str(), phi_0, phi, b_1);

  return phi;
}


void set_rev_bends(void)
{
  double
    dphi_tot = 0e0,
    phi_tot  = 0e0;

  printf("\n   Name              phi      dphi       b_1\n");
  for (int k = 0; k <= globval.Cell_nLoc; k++) { 
    if ((Cell[k].Elem.Pkind == Mpole) && (Cell[k].Elem.M->Pirho != 0e0)) {
      dphi_tot += set_rev_bend(Cell[k]);
      phi_tot += get_elem_phi(Cell[k].Fnum, Cell[k].Knum);
    }
  }
  printf("\nphi_tot = %5.3f Sigma = %9.3e\n", phi_tot, dphi_tot);
}


long int get_orbit(void)
{
    long int        lastpos;
    ss_vect<double> ps;

    ps.zero();
    Cell_Pass(0, globval.Cell_nLoc, ps, lastpos);
    if (lastpos == globval.Cell_nLoc)
      cout << scientific << setprecision(3) << "\n" << setw(11) << ps << "\n";
    else
      printf("\nparticle lost at element %ld (%ld)\n",
	     lastpos, globval.Cell_nLoc);
    return lastpos;
}


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;

  globval.mat_meth       = false;
}


int main(int argc, char *argv[])
{
  long int lastpos;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  Ring_GetTwiss(true, 0e-3);
  printglob();

  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!false) {
    GetEmittance(ElemIndex("cav"), true, true);
    scan_f_RF("cav", 0.5e3, 10);
  }

  if (false) {
    GetEmittance(ElemIndex("cav"), true, true);

    set_rev_bends();

    lastpos = get_orbit();

    prt_beampos1("cod.out", lastpos);

    assert(false);

    Ring_GetTwiss(true, 0e-3);
    printglob();

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);

    GetEmittance(ElemIndex("cav"), true, true);
  }
}
