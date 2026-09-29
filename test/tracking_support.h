// Assertions, error bounds and fixtures shared by the tracking regression
// cases.  Each tracking_*.cc file holds one group of cases and the
// run_*_tests() function that main() calls.
#ifndef TRACKING_SUPPORT_H
#define TRACKING_SUPPORT_H

#define NO 1

#include "tracy_lib.h"

namespace tracking_test {

const char* const coordinate_names[] = {"x", "px", "y", "py", "delta", "ct"};

// Bound for comparisons whose only error source is roundoff, scaled by
// max(1, |value|).  expect_value, expect_state and expect_matrix use it.
const double tolerance = 1e-14;

// Order checks compare a measured log2 ratio against an integer exponent.
const double order_tolerance = 0.15;

extern int failure_count;

// --- assertions ------------------------------------------------------------
//
// Every assertion rejects a nonfinite expected value, actual value, or bound
// before comparing: a comparison with NaN is false, so
// `|actual-expected| > allowed` alone would pass it.

void expect_near(const char* element, const char* quantity,
                 const double expected, const double actual,
                 const double allowed);
void expect_value(const char* element, const char* quantity,
                  const double expected, const double actual);
void expect_at_most(const char* element, const char* quantity,
                    const double bound, const double actual);
void expect_at_least(const char* element, const char* quantity,
                     const double bound, const double actual);
void expect_positive(const char* element, const char* quantity,
                     const double actual);
void expect_index(const char* test, const long expected, const long actual);
void expect_state(const char* element, const ss_vect<double>& expected,
                  const ss_vect<double>& actual);
void expect_matrix(const char* element, const Matrix& expected,
                   const Matrix& actual);
void expect_symplectic(const char* element, const Matrix& map);

// Set TRACKING_TEST_MEASURE=1 to print every residual that a bound was
// derived from.  Re-measure with this before changing any bound.
void report(const char* test, const char* quantity, const double value);

// --- residuals -------------------------------------------------------------

// Maximum for residual reductions.  std::max(largest, NaN) returns largest,
// which would drop a NaN residual before any assertion saw it; this keeps it.
double larger(const double largest, const double value);
double max_deviation(const ss_vect<double>& left,
                     const ss_vect<double>& right);

// In Tracy's (x, px, y, py, delta, ct) convention, Omega has canonical x-px
// and y-py blocks plus a delta-ct block in which delta plays the coordinate
// and ct the momentum: the drift gives Delta ct = -dH/ddelta, like
// Delta px = -dH/dx.  symplectic_product(M, i, j) is (M^T Omega M)[i][j].
double symplectic_product(const Matrix& map, const int row, const int column);
double symplectic_identity(const int row, const int column);

// --- configuration and oracles ---------------------------------------------

// Pins every global the tracking path branches on, to Read_Lattice's
// defaults for a ring (physlib.cc).  run_test calls it before each case, so a
// flag set in one case cannot leak into the next.
void configure_tracking();

ss_vect<double> make_state(const double x, const double px, const double y,
                           const double py, const double delta,
                           const double ct);

// Paraxial drift Hamiltonian: x += L px/(1+delta), likewise for y,
// and ct += L (px^2+py^2)/(2 (1+delta)^2).
void apply_drift_oracle(const double length, ss_vect<double>& state);

double radians(const double degrees);

// --- misalignment ----------------------------------------------------------
//
// The cases write a misalignment, and call the frame transforms directly,
// only through these three functions, so that a change in how Tracy stores a
// misalignment touches them alone.  The oracles read dS and dT, which remain
// the reporting view of the misalignment.

void set_misalignment(CellType& cell, const double dx, const double dy,
                      const double roll);

template<typename T>
void to_local(ss_vect<T>& state, const CellType& cell, const double c0,
              const double c1, const double s1)
{ GtoL(state, cell.dS, cell.dT, c0, c1, s1); }

template<typename T>
void to_global(ss_vect<T>& state, const CellType& cell, const double c0,
               const double c1, const double s1)
{ LtoG(state, cell.dS, cell.dT, c0, c1, s1); }

// --- fixtures --------------------------------------------------------------

CellType make_drift_cell(const char* name, const double length);
CellType make_thin_multipole(const char* name, const int order,
                             const double strength);
CellType make_thick_quadrupole(const char* name, const double length,
                               const double gradient, const int slices);
CellType make_bend(const char* name, const double length,
                   const double curvature, const double entrance_degrees,
                   const double exit_degrees, const int slices);
CellType make_cavity(const char* name);
// Cell[0..3]: START (zero length), D1, D2 and D3 (0.2 m each), apertures of
// +-1e-2 m except D2's upper limit in `plane` (X_ or Y_), 6e-3 m.  D3 lies
// downstream of the loss, so a pass that keeps tracking after it shows.
void set_up_aperture_line(const int plane);

// --- runners ---------------------------------------------------------------

void run_test(const char* name, void (*test)());

void run_element_tests();

} // namespace tracking_test

#endif
