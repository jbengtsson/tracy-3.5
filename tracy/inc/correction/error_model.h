#ifndef CORRECTION_ERROR_MODEL_H
#define CORRECTION_ERROR_MODEL_H

// Machine error model — loading engineering tolerances into the lattice.
//
// The standard file-based error loaders, extracted from param_data_type. Their
// only coupling to the god-class was the config file paths (ae/fe/ap_file),
// now explicit arguments; the rest is globals (Cell/globval, lsoc bpms_) and
// free functions. State-free, so any app/corrector can load errors directly.
// See error_models.md and correction_refactor.md.
//
// The girder "cormisal" model (GirderSetup/SetCorMis/CorMis_in) and Align_BPMs
// still live on param_data_type; they own the Girder/Lattice arrays and are a
// later slice.

namespace corr {

// Read alignment tolerances from ae_file and misalign the lattice. seed selects
// the RNG stream; Scale_it/Scale scale the rms values; new_rnd re-draws.
void LoadAlignTol(const std::string &ae_file, const bool Scale_it,
		  const double Scale, const bool new_rnd, const int seed);

// Read multipole field errors from fe_file and apply them.
void LoadFieldErr(const std::string &fe_file, const bool Scale_it,
		  const double Scale, const bool new_rnd);

// Read physical apertures from ap_file (scaled by scl_x/scl_y) and set them.
void LoadApers(const std::string &ap_file, const double scl_x,
	       const double scl_y);

// Read per-multipole misalignments from the fixed file "cormis.txt".
void ReadCorMis(const bool Scale_it, const double Scale);

}  // namespace corr

#endif  // CORRECTION_ERROR_MODEL_H
