#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void compute_fel_dnu(const string &marker)
{
  const std::vector<string>
    oct_name = {"oct_1", "oct_2", "oct_3", "oct_4", "oct_5",
		"oct_5", "oct_4", "oct_3", "oct_2", "oct_1"};
  const std::vector<double>
    oct_kid  = {   1,        1,       1,       1,      2,
		   2,        1,       1,       1,      1};


  std::vector<int>                   loc;
  std::vector< std::vector<double> > dnu;
  std::vector<double>                pair;

  for (auto k = 0; k < oct_name.size(); k++)
    loc.push_back(Elem_GetPos(ElemIndex(oct_name[k].c_str()), oct_kid[k]));
  for (auto k = 0; k < oct_name.size(); k++) {
    pair.clear();
    auto cell = Cell[loc[k]];
    pair.push_back(cell.Nu[X_]);
    pair.push_back(cell.Nu[Y_]);
    dnu.push_back(pair);
  }

  printf("\n   loc  name       dnu_x     dnu_y\n");
  printf("   %3d  %.8s  %6.3f    %6.3f\n",
	 loc[0], Cell[loc[0]].Elem.PName, 0e0, 0e0);
  for (auto k = 1; k < oct_name.size(); k++) {
    printf("   %3d  %.8s  %6.3f    %6.3f\n",
	   loc[k], Cell[loc[k]].Elem.PName, dnu[k][X_]-dnu[k-1][X_],
	   dnu[k][Y_]-dnu[k-1][Y_]);
  }
}


void get_kick_map
(const string &file_name, const int nx, const int ny, const long int loc1,
 const long int loc2, const double Ax, const double Ay)
{
  const double Brho = globval.Energy*1e9/c0;

  long int        lastpos;
  ss_vect<double> ps0, ps1, dps_map;
  ofstream        outf;

  cout << "\nget_kick_map:\n  " << setw(8) << Cell[loc1].Elem.PName << " -> "
       << setw(8) << Cell[loc2].Elem.PName << "\n";

  file_wr(outf, file_name.c_str());

  dps_map[x_] = 2e0*Ax/(nx-1e0);
  dps_map[y_] = 2e0*Ay/(ny-1e0);

  outf << "# Author:" << "\n";
  outf << "# Title" << "\n";
  outf << "# Cell Length [m]" << "\n";
  outf << fixed << setprecision(5) << Cell[loc2].S-Cell[loc1].S
       << "\n";
  outf << "# Number of Horizontal Points" << "\n";
  outf << nx << "\n";
  outf << "# Number of Vertical Points" << "\n";
  outf << ny << "\n";

  outf << "# Horizontal 2nd Order Kick [T2m2]" << "\n";
  outf << "START" << "\n";

  for (auto i1 = 0; i1 < nx; i1++)
    outf << scientific << setprecision(5) << setw(13)
	 << -Ax+i1*dps_map[x_];
  outf << "\n";
  cout << "  scanning horizontal plane:\n    ";
  ps0.zero();
  for (auto i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];
    cout << ".";
    outf << scientific << setprecision(5)
	 << setw(13) << ps0[y_];
    for (auto i2 = 0; i2 < nx; i2++) {
    ps0[x_] = -Ax + i2*dps_map[x_];
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      if (lastpos == loc2) {
	ps1 -= ps0;
	outf << scientific << setprecision(5)
	     << setw(13) << sqr(Brho)*ps1[px_];
      } else
	outf << scientific << setprecision(5) << setw(13) << NAN;
    }
    outf << "\n";
  }
  cout << "\n";

  outf << "# Vertical 2nd Order Kick [T2m2]" << "\n";
  outf << "START" << "\n";

  for (auto i1 = 0; i1 < nx; i1++)
    outf << scientific << setprecision(5) << setw(13)
	 << -Ax+i1*dps_map[x_];
  outf << "\n";

  cout << "  scanning vertical plane:\n    ";
  ps0.zero();
  for (auto i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];
    cout << ".";
    outf << scientific << setprecision(5) << setw(13) << ps0[y_];

    for (auto i2 = 0; i2 < nx; i2++) {
      ps0[x_] = -Ax + i2*dps_map[x_];
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      if (lastpos == loc2) {
	ps1 -= ps0;
	outf << scientific << setprecision(5)
	     << setw(13) << sqr(Brho)*ps1[py_];
      } else
	outf << scientific << setprecision(5) << setw(13) << NAN;
    }
    outf << "\n";
  }
  cout << "\n";

  outf.close();
}


void get_map_2D
(const string &file_name, const int nx, const int ny, const long int loc1,
 const long int loc2, const double Ax, const double Ay)
{
  long int        lastpos;
  ss_vect<double> ps0, ps1, dps_map;
  ofstream        outf;

  cout << "\nget_map_2D:\n  " << setw(8) << Cell[loc1].Elem.PName << " -> "
       << setw(8) << Cell[loc2].Elem.PName << "\n";

  file_wr(outf, file_name.c_str());

  dps_map[x_] = 2e0*Ax/(nx-1e0);
  dps_map[y_] = 2e0*Ay/(ny-1e0);

  outf << scientific << setprecision(5)
       << "# nx = " << nx << "     Ax = " << Ax
       << "     dx = " << dps_map[x_] << "\n";
  outf << scientific << setprecision(5)
       << "# ny = " << ny << "     Ay = " << Ay
       << "     dy = " << dps_map[y_] << "\n";
  outf << "#" << "\n";

  // cout << "\n";
  ps0.zero();
  for (auto i1 = 0; i1 < nx; i1++) {
    ps0[x_] = -Ax + i1*dps_map[x_];
    for (auto i2 = 0; i2 < ny; i2++) {
      ps0[y_] = -Ay + i2*dps_map[y_];
//       cout << setw(3) << i1+1 << setw(3) << i2+1
// 	   << scientific << setprecision(5)
// 	   << setw(13) << ps0[x_] << setw(13) << ps0[y_] << "\n";
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      ps1 -= ps0;
      outf << scientific << setprecision(5)
	   << setw(13) << ps0[x_] << setw(13) << ps0[y_]
	   << setw(13) << ps1[px_] << setw(13) << ps1[py_] << "\n";
    }
    outf << "\n";
  }
  outf.close();
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
  globval.EPU            = true;
}


int main(int argc, char *argv[])
{

  globval.mat_meth = false;

  FieldMap_filetype = 6;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  no_sxt();

  Ring_GetTwiss(true, 0e0);
  printglob();

  prtmfile("flat_file.dat");
  prt_lat("linlat1.out", globval.bpm, true);
  prt_lat("linlat.out", globval.bpm, true, 10);

  if (!true)
    GetEmittance(ElemIndex("cav"), false, true);

  if (false)
    compute_fel_dnu("fel_mark");

  if (!false) {
    const int
      n[] = {51, 51}, 
      i_0 = 18;
    const double
      A[] = {20e-3, 4e-3};

    get_kick_map
      ("helical_und.dat", n[X_], n[Y_], i_0, i_0, A[X_], A[Y_]);

    get_map_2D
      ("helical_und_2D.dat", n[X_], n[Y_], i_0, i_0, A[X_], A[Y_]);
  }
}
