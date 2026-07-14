// Shared correction primitives — see correction/corr_utils.h.

namespace corr {

void zero_mult(std::vector<double> bn_an[])
{
  bn_an[Sext].clear();
  for (auto k = 0; k <= globval.Cell_nLoc; k++) {
    if (Cell[k].Elem.Pkind == Mpole) {
      bn_an[HOMmax+Sext].push_back(Cell[k].Elem.M->PB[HOMmax+Sext]);
      Cell[k].Elem.M->PB[HOMmax+Sext] = 0e0;
    }
  }
  printf("\nparam_data_type::zero_mult: zeroed b_3 for %d multipoles.\n",
	 (int)bn_an[HOMmax+Sext].size());
}


void restore_mult(std::vector<double> bn_an[])
{
  auto k = 0;
  for (auto j = 0; j <= globval.Cell_nLoc; j++) {
    if (Cell[j].Elem.Pkind == Mpole) {
      Cell[j].Elem.M->PB[HOMmax+Sext] = bn_an[HOMmax+Sext][k];
      k++;
    }
  }
  printf("\nparam_data_type::restore_mult:restored b_3 for %d multiupoles.\n",
	 (int)bn_an[HOMmax+Sext].size());
}


void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		   double s_dnu[], const bare_optics &bare)
{
  int       k;
  long int  j, ind;
  double    dbeta, dnu;

  const int n_sext = bare.n_sext;

  Ring_GetTwiss(false, 0.0);

  for (k = 0; k <= 1; k++) {
    m_dbeta[k] = 0.0; s_dbeta[k] = 0.0; m_dnu[k] = 0.0; s_dnu[k] = 0.0;
  }

  for (j = 0; j < n_sext; j++) {
    ind = bare.sexts[j];
    for (k = 0; k <= 1; k++) {
      dbeta = (Cell[ind].Beta[k]-bare.betas0_[j][k])/bare.betas0_[j][k];
      m_dbeta[k] += dbeta; s_dbeta[k] += sqr(dbeta);
      dnu = Cell[ind].Nu[k] - bare.nus0_[j][k];
      m_dnu[k] += dnu; s_dnu[k] += sqr(dnu);
    }
  }

  for (k = 0; k <= 1; k++) {
    m_dbeta[k] /= n_sext; m_dnu[k] /= n_sext;
    s_dbeta[k] = sqrt((s_dbeta[k]-n_sext*sqr(m_dbeta[k]))/(n_sext-1));
    s_dnu[k] = sqrt((s_dnu[k]-n_sext*sqr(m_dnu[k]))/(n_sext-1));
  }
}

}  // namespace corr
