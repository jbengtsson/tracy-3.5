// Closed-orbit-distortion correction driver — see correction/orbit_corr.h.

namespace corr {

void ini_COD_corr(const int n_bpm_Fam, const std::string bpm_names[],
		  const int n_hcorr_Fam, const std::string hcorr_names[],
		  const int n_vcorr_Fam, const std::string vcorr_names[],
		  const bool svd)
{
  int i, j, Fnum, n_bpm, n_hcorr, n_vcorr;

  n_bpm = 0;
  for (i = 0; i < n_bpm_Fam; i++)
    n_bpm += GetnKid(ElemIndex(bpm_names[i]));

  n_hcorr = 0;
  for (i = 0; i < n_hcorr_Fam; i++)
    n_hcorr += GetnKid(ElemIndex(hcorr_names[i]));

  n_vcorr = 0;
  for (i = 0; i < n_vcorr_Fam; i++)
    n_vcorr += GetnKid(ElemIndex(vcorr_names[i]));

  long int  bpms[n_bpm], hcorrs[n_hcorr], vcorrs[n_vcorr];

  n_bpm = 0;
  for (i = 0; i < n_bpm_Fam; i++) {
    Fnum = ElemIndex(bpm_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      bpms[n_bpm++] = Elem_GetPos(Fnum, j);
  }

  n_hcorr = 0;
  for (i = 0; i < n_hcorr_Fam; i++) {
    Fnum = ElemIndex(hcorr_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      hcorrs[n_hcorr++] = Elem_GetPos(Fnum, j);
  }

  n_vcorr = 0;
  for (i = 0; i < n_vcorr_Fam; i++) {
    Fnum = ElemIndex(vcorr_names[i]);
    for (j = 1; j <= GetnKid(Fnum); j++)
      vcorrs[n_vcorr++] = Elem_GetPos(Fnum, j);
  }

  std::cout << std::endl;
  std::cout << "ini_COD_corr: n_bpm = " << n_bpm << ", n_hcorr = " << n_hcorr
       << ", n_vcorr = " << n_vcorr << std::endl;

  gcmat(n_bpm, bpms, n_hcorr, hcorrs, 1, svd);
  gcmat(n_bpm, bpms, n_vcorr, vcorrs, 2, svd);

  if (true) {
    gtcmat(n_bpm, bpms, n_hcorr, hcorrs, 1, svd);
    gtcmat(n_bpm, bpms, n_vcorr, vcorrs, 2, svd);
  }
}


void orbit_corr::alloc(const std::vector<string> &bpm_Fam_names,
		       const std::vector<string> corr_Fam_names[])
{
  ::cod_ini(bpm_Fam_names, corr_Fam_names, orb_corr);
}


void orbit_corr::dealloc(void)
{
  int j;

  for (j = 0; j < 2; j++)
    orb_corr[j].dealloc();
}


void orbit_corr::prt_svdmat(void)
{
  orb_corr[X_].prt_svdmat();
  orb_corr[Y_].prt_svdmat();
}


bool orbit_corr::cod_corr(const orbit_cfg &cfg, const bare_optics &bare,
			  const int n_cell, const double scl,
			  const double h_maxkick, const double v_maxkick)
{
  bool                cod = false;
  long int            lastpos;
  double              m_dbeta[2], s_dbeta[2], m_dnu[2], s_dnu[2];
  ss_vect<double>     ps;
  std::vector<double> bn_an[2*HOMmax+1];  // local save buffer for zero/restore

  orb_corr[X_].clr_trims();
  orb_corr[Y_].clr_trims();

  zero_mult(bn_an);

  cod = getcod(0e0, lastpos);
  printf("\ncorr::cod_corr: %d\n", cod);

  if (!cod) {
    printf("  could not find closed orbit; threading beam\n");
    printf("  corr::cod_corr: n_cell = %d loc_Fam_name = \"%s\"\n",
	   n_cell, cfg.loc_Fam_name.c_str());

    orb_corr[X_].clr_trims(); orb_corr[Y_].clr_trims();
    thread_beam(n_cell, cfg.loc_Fam_name, cfg.bpm_Fam_names, cfg.corr_Fam_names,
		cfg.n_thread, scl);
    //prt_cod("codt.out", globval.bpm, true);
  }

  cod = ::cod_correct(cfg.n_orbit, scl, orb_corr);

  restore_mult(bn_an);

  get_dbeta_dnu(m_dbeta, s_dbeta, m_dnu, s_dnu, bare);
  printf("\ncod_corr: rms dbeta_x/beta_x = %4.2f%%"
	 ",   dbeta_y/beta_y = %4.2f%%\n",
	 1e2*s_dbeta[X_], 1e2*s_dbeta[Y_]);
  printf("          rms dnu_x          = %7.5f, dnu_y          = %7.5f\n",
	 s_dnu[X_], s_dnu[Y_]);

  prt_cod("cod.out", globval.bpm, true);

  return cod;
}


void orbit_corr::Orb_and_Trim_Stat(void)
{
  int     i, j;
  int     SextCounter = 0;
  int     bins[5]     = { 0, 0, 0, 0, 0 };
  double  bin         = 40.0e-6;              // bin size for stat
  double  tr;                                 // trim strength
  Vector2 Sext_max, Sext_sigma, TrimMax, orb;

  for (j = 0; j < 2; j++) {
   Sext_max[j] = Sext_sigma[j] = TrimMax[j] = 0e0;
  }
  SextCounter = 0;
  for (i = 0; i <= globval.Cell_nLoc; i++) {
    if ((Cell[i].Elem.Pkind == Mpole) && (Cell[i].Elem.M->n_design == Sext)) {
      SextCounter++;
      for (j = 0; j < 2; j++) {
	orb[j] = Cell[i].BeamPos[2*j];
	Sext_sigma[j] += sqr(orb[j]);
	if (fabs(orb[j]) > Sext_max[j]) Sext_max[j] = fabs(orb[j]);
      }
      j = (int) (sqrt(sqr(orb[X_])+sqr(orb[Y_]))/bin);
      if (j > 4) j = 4;
      if (j >= 0)
	bins[j]++;
      else
	printf("\nOrb_and_Trim_Stat: negative bin %d\n", j);
    } // sextupole handling
  } // looking throught the cells

  // Trim handling.
  for (j = 0; j < 2; j++)
    for (i = 0; i < (int)orb_corr[j].corrs.size(); i++) {
      if (j == 0)
	tr = Cell[orb_corr[j].corrs[i]].Elem.M->PBpar[HOMmax+Dip];
      else
	tr = Cell[orb_corr[j].corrs[i]].Elem.M->PBpar[HOMmax-Dip];
      TrimMax[j] = max(fabs(tr), TrimMax[j]);
    }


  for (j = 0; j < 2; j++)
    Sext_sigma[j] = sqrt(Sext_sigma[j]/SextCounter);
  printf("In sextupoles maximal horizontal orbit is:"
	 " %5.3f mm with sigma %5.3f mm\n",
	  1e3*Sext_max[X_], 1e3*Sext_sigma[X_]);
  printf("and maximal vertical orbit is:            "
	 " %5.3f mm with sigma %5.3f mm.\n",
	 1e3*Sext_max[Y_], 1e3*Sext_sigma[Y_]);

  for (i = 0; i < 4;  i++) {
    printf("There are %3d sextupoles with offset between  "
	   " %5.3f mm and %5.3f mm\n",
	   bins[i], i*bin*1e3, (i+1)*bin*1e3);
  }
  printf("There are %3d sextupoles with offset more than %5.3f mm \n",
	 bins[4], 4e3*bin);
  printf("Maximal hcorr is %5.3f mrad, maximal vcorr is %5.3f mrad\n",
	 1e3*TrimMax[X_], 1e3*TrimMax[Y_]);
}

}  // namespace corr
