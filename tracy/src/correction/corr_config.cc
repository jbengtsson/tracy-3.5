// Per-corrector configuration slices and bare-optics reference — see
// correction/corr_config.h. (param.dat parsing lives in param.cc.)

namespace corr {

void bare_optics::capture(void)
{
  // Store optics function values at the sextupoles.
  long int j, k;

  n_sext = 0;
  for (j = 0; j <= globval.Cell_nLoc; j++) {
    if ((Cell[j].Elem.Pkind == Mpole) && (Cell[j].Elem.M->n_design >= Sext)) {
      n_sext++; sexts[n_sext-1] = j;
      for (k = 0; k < 2; k++) {
	betas0_[n_sext-1][k] = Cell[j].Beta[k];
	nus0_[n_sext-1][k] = Cell[j].Nu[k];
      }
    }
  }
}

}  // namespace corr
