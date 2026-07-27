#ifndef CORRECTION_FIT_LAT_H
#define CORRECTION_FIT_LAT_H

// Fit the lattice to design tunes / chromaticities.
//
// The fit knobs are element FAMILIES, named through the tune_fams / chrom_fams
// param.dat keywords and resolved by EXACT (case-insensitive) family name. Any
// number of families may be given: the fit is a 2 x N least-squares solve by SVD
// pseudo-inverse, so a redundant or near-degenerate knob set costs a dropped
// singular direction rather than a failed solve.
//
// The free parameter per family is its TOTAL integrated strength db_nL, spread
// evenly over the family's members. So db_2L / db_3L below are whole-family
// steps, independent of how many magnets the family has or how many periods the
// lattice spans — NOT per-magnet steps.
//
// Targets are ABSOLUTE, and for the tunes they are compared against the raw,
// unwrapped globval.TotalTune.
//
// Both fits iterate to eps, or restore their knobs and return false. Diagnostics
// beyond the one before/after line are behind the global trace flag.
//
// State-free (globals + args only), so any app/corrector can fit the lattice.

namespace corr {

// Fit the ring tunes to (nu_x, nu_y) with the quadrupole families in fams, then
// re-Twiss. Returns false, with the knobs left as found, if a name is not a
// family of the loaded lattice, the ring goes unstable, the closed-orbit finder
// fails, or eps is not reached within imax steps.
bool fit_tune(const std::vector<std::string> &fams, const double nu_x,
	      const double nu_y, const double db_2L = 1e-3,
	      const double eps = 1e-4, const int imax = 10);

// Fit the chromaticities to (chrom_x, chrom_y) with the sextupole families in
// fams, then re-Twiss. Returns false as fit_tune does. Radiation and the cavity
// are off for the duration and restored on every exit path.
bool fit_chrom(const std::vector<std::string> &fams, const double chrom_x,
	       const double chrom_y, const double db_3L = 1e0,
	       const double eps = 1e-4, const int imax = 10);

}  // namespace corr

#endif  // CORRECTION_FIT_LAT_H
