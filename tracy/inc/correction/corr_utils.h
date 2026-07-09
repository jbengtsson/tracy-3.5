#ifndef CORRECTION_CORR_UTILS_H
#define CORRECTION_CORR_UTILS_H

// Shared correction primitives — no param_data_type / god-class state.
//
// These are used across correctors (orbit, coupling, ID) and are expected to be
// reused by a future LOCO module: LOCO fits *linear* optics, so it needs the
// sextupoles off (zero_mult/restore_mult) and reports residual beta-beat / tune
// shift (get_dbeta_dnu). Kept state-free (buffers/references passed in) so any
// caller can use them. See correction_refactor.md.

namespace corr {

// Save every sextupole b_3 into bn_an[HOMmax+Sext] and zero it in the lattice,
// then restore it. Linearizes the machine around a linear-optics operation.
// The same bn_an buffer must be handed to the matching restore_mult call.
void zero_mult(std::vector<double> bn_an[]);
void restore_mult(std::vector<double> bn_an[]);

// RMS beta-beat and tune shift at the sextupoles relative to the bare lattice
// reference (n_sext / sexts / betas0_ / nus0_, captured by get_bare()).
void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		   double s_dnu[], const int n_sext, const int sexts[],
		   const double betas0_[][2], const double nus0_[][2]);

}  // namespace corr

#endif  // CORRECTION_CORR_UTILS_H
