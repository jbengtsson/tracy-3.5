#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the lattice to design tunes / chromaticities.
//
// The fit knobs are element FAMILIES, named through the tune_fams / chrom_fams
// param.dat keywords and resolved by EXACT (case-insensitive) family name.
// Any number of families may be given.
//
// State-free (globals + args only), so any app/corrector can fit the lattice.

namespace corr {

// Fit the ring tunes to (nu_x, nu_y) with the quadrupole families in fams,
// then re-Twiss. Returns false, without fitting, if a name is not a family of
// the loaded lattice.
bool fit_tune(const std::vector<std::string> &fams, const double nu_x,
	      const double nu_y);

// Fit the chromaticities to (chrom_x, chrom_y) with the sextupole families in
// fams, then re-Twiss. Returns false as fit_tune does.
bool fit_chrom(const std::vector<std::string> &fams, const double chrom_x,
	       const double chrom_y);

}  // namespace corr

#endif  // CORRECTION_FIT_LAT_H
