#ifndef CORRECTION_GIRDER_MODEL_H
#define CORRECTION_GIRDER_MODEL_H

// Cormisal (girder-based) alignment error model, selected by n_meth == 1 (vs
// the file-based LoadAlignTol at n_meth == 0). GirderSetup builds a 3-level
// girder tree from the lattice; SetCorMis draws correlated random
// misalignments per girder and translates them to the elements on each girder;
// CorMis_in reads the rms amplitudes from "cormis.dat".
//
// WARNING — known bug: SetCorMis's girder-support interpolation divides by a
// girder's span (gsp[1]-gsp[0]); a zero-span girder -> division by zero -> NaN
// misalignments -> lost beam, even at zero error amplitude. GirderSetup can
// create zero-span level-3 girders from runs of consecutive zero-length magnets
// (e.g. adjacent zero-length correctors). See the inline TODO at the
// interpolation. The girder->element translation is also not physics-validated.

// Sizing limits. Global scope, not namespaced, because dynap.cc references
// iseednrmax unqualified.
const int igrmax = 2000, ilatmax = 10000, iseednrmax = 20;

namespace corr {

typedef struct girdertype {
  double gsp[2], gdx[2], gdy[2], gdt;
  long ilat[2], igir[2], gco[2], level;
} girdertype;

typedef struct latticetype {
  long igir;
  double smid;
} latticetype;

// Owns the girder-tree state shared by GirderSetup and SetCorMis.
struct girder_model {
  girdertype  Girder[igrmax];
  long        NGirderLevel[3];
  latticetype Lattice[ilatmax];

  void GirderSetup();
  void SetCorMis(double gxrms, double gyrms, double gtrms, double jxrms,
		 double jyrms, double exrms, double eyrms, double etrms,
		 double rancutx, double rancuty, double rancutt, long iseed);
};

// Stateless: read the girder/joint/element/BPM rms amplitudes and seeds from
// "cormis.dat" into the out-params (converting micron/micro-deg -> SI).
void CorMis_in(double *gdxrms, double *gdzrms, double *gdarms, double *jdxrms,
	       double *jdzrms, double *edxrms, double *edzrms, double *edarms,
	       double *bdxrms, double *bdzrms, double *bdarms, double *rancutx,
	       double *rancuty, double *rancutt, long *iseed, long *iseednr);

}  // namespace corr

#endif  // CORRECTION_GIRDER_MODEL_H
