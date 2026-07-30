#ifndef CORRECTION_CORR_UTILS_H
#define CORRECTION_CORR_UTILS_H

// Shared correction primitives, used across correctors (orbit, coupling, ID)
// and by a LOCO module (not yet implemented).
//
// State-free — buffers and references are passed in — so any caller can use
// them without owning a corrector.

namespace corr {

// Save every sextupole b_3 into bn_an[HOMmax+Sext] and zero it in the lattice,
// then restore it. Linearizes the machine around a linear-optics operation.
// The same bn_an buffer must be handed to the matching restore_mult call.
void zero_mult(std::vector<double> bn_an[]);
void restore_mult(std::vector<double> bn_an[]);

// RMS beta-beat and tune shift at the sextupoles, relative to the bare-lattice
// reference captured by bare_optics::capture().
void get_dbeta_dnu(double m_dbeta[], double s_dbeta[], double m_dnu[],
		   double s_dnu[], const bare_optics &bare);

}  // namespace corr

#endif  // CORRECTION_CORR_UTILS_H
