#ifndef CORRECTION_LOCO_ORM_H
#define CORRECTION_LOCO_ORM_H

// Measured orbit response matrix (ORM): the closed-orbit response at every BPM
// to a kick at every corrector. LOCO fits the lattice model by driving the
// modelled ORM onto the measured one. The coupling corrector measures only the
// off-diagonal blocks, where a kick in one plane shows up in the other plane's
// BPM readings.

namespace corr {

// Measure one column of the ORM: symmetrically kick the corrector (fnum, knum)
// by +/-kick and finite-difference the resulting closed orbit read at the BPMs.
//
//   kick_type   +Dip for a horizontal kick, -Dip for a vertical one (the sign
//               convention SetdKLpar uses to pick the plane).
//   bpm_loc     lattice indices of the n_bpm BPMs, 0-based (bpm_loc[0..n_bpm-1]).
//   read_coord  which BeamPos component to read: x_ or y_.
//   resp        out, 1-based: resp[1..n_bpm] = d(orbit)/d(kick) [m/rad].
//
// The corrector is restored to its original setting on the way out. The closed
// orbit must exist at both kick settings; a failure aborts via chk_cod.
void measure_orm_column(const int fnum, const int knum, const int kick_type,
			const double kick, const int bpm_loc[], const int n_bpm,
			const int read_coord, double resp[]);

}  // namespace corr

#endif  // CORRECTION_LOCO_ORM_H
