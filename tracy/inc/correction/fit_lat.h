#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the linear lattice to design tunes / chromaticities.
//
// The tune- and chromaticity-fit steps, extracted from err_and_corr_init. The
// fit-knob family names (the tune quadrupoles and the chroma sextupoles) are
// explicit arguments; err_and_corr_init passes the param_data_type tune_fam /
// chrom_fam fields (the tune_fams / chrom_fams param.dat keywords). Families are
// matched by element-name prefix over the global Cell[] lattice.
//
// State-free (globals + args only), so any app/corrector can fit the lattice.

namespace corr {

// Fit the ring tunes to (nu_x, nu_y) with the two quadrupole families whose
// element names start with fam_h / fam_v, then re-Twiss.
void fit_tune(const std::string &fam_h, const std::string &fam_v,
	      const double nu_x, const double nu_y);

// Fit the chromaticities to (chrom_x, chrom_y) with the two sextupole families
// whose element names start with fam_h / fam_v, then re-Twiss.
void fit_chrom(const std::string &fam_h, const std::string &fam_v,
	       const double chrom_x, const double chrom_y);

}  // namespace corr

#endif  // CORRECTION_FIT_LAT_H
