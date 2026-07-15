#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the linear lattice to design tunes / chromaticities.
//
// The tune- and chromaticity-fit steps, extracted from err_and_corr_init. The
// fit-knob family names (the tune quadrupoles and the chroma sextupoles) are now
// explicit arguments instead of literals buried in the loop; err_and_corr_init
// still passes them hardcoded. Families are matched by element-name prefix over
// the global Cell[] lattice, exactly as before.
//
// State-free (globals + args only), so any app/corrector can fit the lattice.
// A future config step promotes the family names to param.dat knobs.

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
