#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


class dnu_dksi {
private:

public:
  double **dnu_mat_inv, **dksi_mat_inv;

  void get_dnu_mat(const std::vector<double> &Fnum, double **dnu_mat_inv);
  void set_dnu(const double nu_x, const double nu_y);
  void get_dksi_mat(const std::vector<double> &Fnum, double **dksi_mat_inv);
  void set_dksi(const double ksi_x, const double ksi_y);
};


void get_dnu_mat(const std::vector<double> &Fnum, double **dnu_mat_inv)
{
}


void set_dnu(const double nu_x, const double nu_y)
{
}


void get_dksi_mat(const std::vector<double> &Fnum, double **dksi_mat_inv)
{
}


void set_dksi(const double ksi_x, const double ksi_y)
{
}


int main(int argc, char *argv[])
{


  // 1: DIAMOND, 2: NSLS-II, 3: Oleg I, 4: Oleg II.
  FieldMap_filetype = 4; sympl = !false;

  reverse_elem = !false;

  trace = false;

  globval.mat_meth = !false;

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  globval.H_exact    = false; globval.quad_fringe    = false;
  globval.Cavity_on  = false; globval.radiation      = false;
  globval.emittance  = false; globval.IBS            = false;
  globval.pathlength = false; globval.bpm            = 0;
  globval.Cart_Bend  = false; globval.dip_edge_fudge = true;


}
