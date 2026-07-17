#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the lattice to design tunes / chromaticities.
//
// The fit-knob family names — the tune quadrupoles and the chroma sextupoles —
// are explicit arguments through the tune_fams / chrom_fams param.dat keywords.
// Families are matched by element-name PREFIX over the global Cell[] lattice,
// so a family name also matches any longer name starting with it.
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
