#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the lattice to design tunes / chromaticities.
//
// The knobs are element families, named through the tune_fams / chrom_fams
// param.dat keywords and resolved by exact (case-insensitive) family name. Any
// number may be given: the fit is a 2 x N solve by SVD. db_2L / db_3L are
// whole-family integrated-strength steps, not per-magnet ones; targets are
// absolute.
//
// State-free (globals + args only), so any app/corrector can fit the lattice.

namespace corr {

// Fit the ring tunes to (nu_x, nu_y) with the quadrupole families in fams, then
// re-Twiss. Returns false on an unknown family, instability, a lost closed
// orbit, or eps not reached in imax steps; the first three restore the knobs,
// the last keeps the best iterate.
bool fit_tune(const std::vector<std::string> &fams, const double nu_x,
	      const double nu_y, const double db_2L = 1e-3,
	      const double eps = 1e-4, const int imax = 10);

// As fit_tune, with the sextupole families in fams. Radiation and the cavity
// are off for the duration and restored on every exit path.
bool fit_chrom(const std::vector<std::string> &fams, const double chrom_x,
	       const double chrom_y, const double db_3L = 1e0,
	       const double eps = 1e-4, const int imax = 10);

}  // namespace corr

#endif  // CORRECTION_FIT_LAT_H
