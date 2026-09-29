// Single-element tracking against analytic oracles, direct 2D frame
// transforms, and a short fixed lattice with aperture loss.  Every case runs
// in the configuration pinned by configure_tracking().
#include "tracking_support.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace tracking_test {
namespace {

void linear_map(const long first_element, const long last_element,
                Matrix& map, const char* const name)
{
  ss_vect<tps> tps_map;
  long last_position = -1;

  tps_map.identity();
  Cell_Pass(first_element, last_element, tps_map, last_position);
  expect_index(name, last_element, last_position);
  getlinmat(ss_dim, tps_map, map);
}

void test_drift()
{
  const double length = 1.7;
  ss_vect<double> actual =
    make_state(2.1e-3, -3.2e-4, -1.4e-3, 2.5e-4, 1.2e-2, 8.0e-4);
  ss_vect<double> expected = actual;

  apply_drift_oracle(length, expected);

  Drift(length, actual);
  expect_state("drift", expected, actual);
}

void test_quadrupole()
{
  const double integrated_strength = 0.83;
  CellType cell = make_thin_multipole("quadrupole", Quad,
                                      integrated_strength);
  ss_vect<double> actual =
    make_state(2.3e-3, -4.1e-4, -1.7e-3, 3.6e-4, 7.0e-3, -2.0e-4);
  ss_vect<double> expected = actual;

  // From B_y+iB_x = b_2 (x+iy): Delta px = -b_2 x,
  // Delta py = b_2 y for an integrated normal quadrupole.  The kick does not
  // depend on delta: px and py are normalised to the reference momentum, so
  // the momentum dependence of the focusing enters through the 1/(1+delta) in
  // the drifts.
  expected[px_] -= integrated_strength*expected[x_];
  expected[py_] += integrated_strength*expected[y_];

  Mpole_Pass(cell, actual);
  expect_state("quadrupole", expected, actual);
  std::free(cell.Elem.M);
}

void test_quadrupole_linear_map()
{
  const double integrated_strength = 0.83;
  Matrix actual, expected;

  Cell[0] = make_thin_multipole("quadrupole map", Quad,
                                integrated_strength);
  UnitMat(ss_dim, expected);
  expected[px_][x_] = -integrated_strength;
  expected[py_][y_] = integrated_strength;

  // The TPSA path starts from the six-dimensional identity map.  Its
  // Jacobian must match the same independent thin-lens quadrupole matrix as
  // the double-particle kick test above.
  linear_map(0, 0, actual, "quadrupole TPSA map survival");
  expect_matrix("quadrupole TPSA map", expected, actual);
  expect_symplectic("quadrupole analytic map", expected);
  expect_symplectic("quadrupole TPSA map", actual);

  std::free(Cell[0].Elem.M);
}

void test_sextupole()
{
  const double integrated_strength = -12.5;
  CellType cell = make_thin_multipole("sextupole", Sext,
                                      integrated_strength);
  ss_vect<double> actual =
    make_state(3.0e-3, 1.1e-4, -2.0e-3, -2.7e-4, -9.0e-3, 4.0e-4);
  ss_vect<double> expected = actual;

  // From B_y+iB_x = b_3 (x+iy)^2: Delta px = -b_3 (x^2-y^2),
  // Delta py = 2 b_3 xy for an integrated normal sextupole.
  expected[px_] -= integrated_strength
                   *(expected[x_]*expected[x_]-expected[y_]*expected[y_]);
  expected[py_] += 2.0*integrated_strength*expected[x_]*expected[y_];

  Mpole_Pass(cell, actual);
  expect_state("sextupole", expected, actual);
  std::free(cell.Elem.M);
}

// --- thick quadrupole ------------------------------------------------------
//
// With no curvature the body integrates
//   H = (px^2+py^2)/(2 (1+delta)) + b2 (x^2-y^2)/2,
// the paraxial drift and the delta-independent quadrupole kick of the thin
// cases above.  Its exact solution, with w = sqrt(b2/(1+delta)) for b2 > 0:
//   x = x0 cos(ws) + px0/((1+delta) w) sin(ws), px = (1+delta) x',
//   y = y0 cosh(ws) + py0/((1+delta) w) sinh(ws), py = (1+delta) y',
// and ct' = (x'^2+y'^2)/2, integrated in closed form below.  The only error
// left is the integrator's, so these cases need no truncated oracle.

const double quadrupole_length = 0.3;
const double quadrupole_gradient = 2.3;

void apply_thick_quadrupole_oracle(const double length, const double b2,
                                   ss_vect<double>& state)
{
  const double p = 1.0+state[delta_];
  const double w = std::sqrt(b2/p), wL = w*length;
  const double a = state[x_], b = state[px_]/(p*w);
  const double c = state[y_], d = state[py_]/(p*w);
  const double cs = std::cos(wL), sn = std::sin(wL);
  const double ch = std::cosh(wL), sh = std::sinh(wL);

  state[x_] = a*cs+b*sn;
  state[px_] = p*w*(-a*sn+b*cs);
  state[y_] = c*ch+d*sh;
  state[py_] = p*w*(c*sh+d*ch);
  // Integrals of x'^2 and y'^2 over the length.
  const double horizontal =
    w*w*(a*a*(length/2.0-std::sin(2.0*wL)/(4.0*w))
         +b*b*(length/2.0+std::sin(2.0*wL)/(4.0*w))-a*b*sn*sn/w);
  const double vertical =
    w*w*(c*c*(std::sinh(2.0*wL)/(4.0*w)-length/2.0)
         +d*d*(std::sinh(2.0*wL)/(4.0*w)+length/2.0)+c*d*sh*sh/w);
  state[ct_] += (horizontal+vertical)/2.0;
}

ss_vect<double> quadrupole_initial_state(const double delta)
{
  return make_state(1.0e-3, -2.0e-4, -8.0e-4, 1.5e-4, delta, 3.0e-4);
}

// Distance between tracked and exact motion through the thick quadrupole.
double thick_quadrupole_residual(const double delta, const int slices)
{
  CellType cell = make_thick_quadrupole("thick quadrupole",
                                        quadrupole_length,
                                        quadrupole_gradient, slices);
  ss_vect<double> tracked = quadrupole_initial_state(delta);
  ss_vect<double> exact = tracked;

  apply_thick_quadrupole_oracle(quadrupole_length, quadrupole_gradient,
                                exact);
  Mpole_Pass(cell, tracked);
  std::free(cell.Elem.M);
  return max_deviation(tracked, exact);
}

// The integrator is fourth order (test_thick_quadrupole_integrator_order),
// so the error at 32 slices is the error at 8 slices times (8/32)^4, up to
// sixth-order terms.  The bound allows twice that; on 2026-09-29 the
// residual at 32 slices was 2.6e-12 against a bound of 5.2e-12.  An oracle
// or physics error does not shrink with the slice length, so it fails this
// bound even where it is small.
double thick_quadrupole_bound(const double delta)
{
  return 2.0*thick_quadrupole_residual(delta, 8)*std::pow(8.0/32.0, 4);
}

void test_thick_quadrupole()
{
  const double residual = thick_quadrupole_residual(0.0, 32);

  report("thick quadrupole", "distance from exact motion", residual);
  report("thick quadrupole", "bound", thick_quadrupole_bound(0.0));
  expect_at_most("thick quadrupole", "distance from exact motion",
                 thick_quadrupole_bound(0.0), residual);
}

void test_thick_quadrupole_off_momentum()
{
  const double deltas[] = {-1.0e-2, 1.0e-2};

  for (const double delta : deltas) {
    const char* const name = delta < 0.0 ? "thick quadrupole, delta < 0"
                                         : "thick quadrupole, delta > 0";
    const double residual = thick_quadrupole_residual(delta, 32);
    // Guard against a vacuous pass: the exact motion at delta must differ
    // from the on-momentum motion by far more than the bound, so a lost or
    // doubled 1/(1+delta) cannot pass.
    ss_vect<double> exact = quadrupole_initial_state(delta);
    ss_vect<double> on_momentum = exact;

    apply_thick_quadrupole_oracle(quadrupole_length, quadrupole_gradient,
                                  exact);
    on_momentum[delta_] = 0.0;
    apply_thick_quadrupole_oracle(quadrupole_length, quadrupole_gradient,
                                  on_momentum);
    on_momentum[delta_] = delta;

    report(name, "distance from exact motion", residual);
    expect_at_most(name, "distance from exact motion",
                   thick_quadrupole_bound(delta), residual);
    report(name, "chromatic effect", max_deviation(exact, on_momentum));
    expect_at_least(name, "chromatic effect",
                    1.0e3*thick_quadrupole_bound(delta),
                    max_deviation(exact, on_momentum));
  }
}

void test_thick_quadrupole_integrator_order()
{
  const double coarse = thick_quadrupole_residual(1.0e-2, 4);
  const double fine = thick_quadrupole_residual(1.0e-2, 8);

  // Against the exact solution, doubling the slice count must cut the error
  // by 2^4, as for the sector bend, here on the path without curvature.
  expect_positive("thick quadrupole", "residual", fine);
  report("thick quadrupole", "integrator order",
         std::log(coarse/fine)/std::log(2.0));
  expect_near("thick quadrupole", "integrator order", 4.0,
              std::log(coarse/fine)/std::log(2.0), order_tolerance);
}

// Sector-bend geometry shared by the body tests below.
const double bend_length = 0.7;
const double bend_curvature = 0.2;

// Initial direction in (x, px, y, py, delta), normalised so that the largest
// coordinate equals the amplitude scale A passed to the helpers below.  The
// delta component exercises the bend's dispersion.
const double bend_direction[5] = {2.0/3.0, -1.0/3.0, -1.0, 0.5, 0.8};

// The bend body is compared against its first-order map, so the residual is
// the leading neglected term, C*A^2, for an initial amplitude A.  C is set an
// order of magnitude above the residual measured on 2026-09-29: 2.798e-07 at
// A = 1e-3, delta included, i.e. C = 0.28.  Provisional: C is not derived from
// the bend's second-order map, and test_sector_bend_nonlinear_order fixes only
// the exponent.  If test_sector_bend_body fails, re-measure rather than
// raising C, and record the new value here.
const double bend_nonlinear_coefficient = 2.8;

CellType make_sector_bend(const char* name, const int slices)
{
  return make_bend(name, bend_length, bend_curvature, 0.0, 0.0, slices);
}

ss_vect<double> bend_initial_state(const double amplitude)
{
  return make_state(bend_direction[0]*amplitude, bend_direction[1]*amplitude,
                    bend_direction[2]*amplitude, bend_direction[3]*amplitude,
                    bend_direction[4]*amplitude, 0.0);
}

// Exact first-order sector-bend map with no edge angles, dispersion included:
// the delta column is (1-cos)/h for x, sin for px and (theta-sin)/h for ct.
void apply_sector_bend_linear_oracle(ss_vect<double>& state)
{
  const double angle = bend_length*bend_curvature;
  const double cosine = std::cos(angle);
  const double sine = std::sin(angle);
  const double x0 = state[x_];
  const double px0 = state[px_];
  const double delta0 = state[delta_];

  state[x_] = cosine*x0+sine/bend_curvature*px0
              +(1.0-cosine)/bend_curvature*delta0;
  state[px_] = -bend_curvature*sine*x0+cosine*px0+sine*delta0;
  state[y_] += bend_length*state[py_];
  state[ct_] += sine*x0+(1.0-cosine)/bend_curvature*px0
                +(angle-sine)/bend_curvature*delta0;
}

// Distance between tracked and linearised motion.  The oracle is truncated at
// first order, so this is dominated by the bend's physical nonlinearity, not
// by roundoff: it is expected to be O(A^2), and that scaling is asserted by
// test_sector_bend_nonlinear_order below.
double sector_bend_linear_residual(const double amplitude, const int slices)
{
  CellType cell = make_sector_bend("sector bend residual", slices);
  ss_vect<double> tracked = bend_initial_state(amplitude);
  ss_vect<double> linear = tracked;

  apply_sector_bend_linear_oracle(linear);
  Mpole_Pass(cell, tracked);
  std::free(cell.Elem.M);

  return max_deviation(tracked, linear);
}

// Distance between the same trajectory integrated with different step counts.
// This compares Tracy against itself on purpose: it measures the integrator's
// convergence order, independently of any physics oracle.
double sector_bend_slice_deviation(const double amplitude, const int slices,
                                   const int reference_slices)
{
  CellType coarse = make_sector_bend("sector bend coarse", slices);
  CellType fine = make_sector_bend("sector bend fine", reference_slices);
  ss_vect<double> coarse_state = bend_initial_state(amplitude);
  ss_vect<double> fine_state = coarse_state;

  Mpole_Pass(coarse, coarse_state);
  Mpole_Pass(fine, fine_state);
  std::free(coarse.Elem.M);
  std::free(fine.Elem.M);

  return max_deviation(coarse_state, fine_state);
}

void test_sector_bend_body()
{
  const double amplitude = 1.0e-3;
  const double residual = sector_bend_linear_residual(amplitude, 32);

  // Physical amplitude, so the linear oracle is compared against genuinely
  // nonlinear motion.
  const double bound = bend_nonlinear_coefficient*amplitude*amplitude;

  report("sector bend", "distance from linear motion", residual);
  expect_at_most("sector bend", "distance from linear motion", bound,
                 residual);
}

void test_sector_bend_nonlinear_order()
{
  const double amplitude = 1.0e-3;
  const double coarse = sector_bend_linear_residual(amplitude, 32);
  const double fine = sector_bend_linear_residual(amplitude/2.0, 32);

  // Halving the amplitude must quarter the departure from linear motion.
  // This constrains the leading nonlinear term itself, which no
  // single-amplitude comparison can do.
  expect_positive("sector bend", "nonlinear residual", fine);
  report("sector bend", "nonlinear order", std::log(coarse/fine)/std::log(2.0));
  expect_near("sector bend", "nonlinear order", 2.0,
              std::log(coarse/fine)/std::log(2.0), order_tolerance);
}

void test_sector_bend_slice_convergence()
{
  const double amplitude = 1.0e-3;
  const int reference_slices = 512;
  const double coarse = sector_bend_slice_deviation(amplitude, 4,
                                                    reference_slices);
  const double fine = sector_bend_slice_deviation(amplitude, 8,
                                                  reference_slices);

  // Tracy integrates thick elements with the fourth-order Forest-Ruth
  // composition (c_1, c_2, d_1, d_2 in t2elem.cc), so doubling the step count
  // must cut the error by about 2^4.  A regression in those coefficients
  // fails here even when every fixed-tolerance comparison still passes.
  expect_positive("sector bend", "slice deviation", fine);
  report("sector bend", "integrator order",
         std::log(coarse/fine)/std::log(2.0));
  expect_near("sector bend", "integrator order", 4.0,
              std::log(coarse/fine)/std::log(2.0), order_tolerance);
}

void test_bend_edges()
{
  const double curvature = 0.31;
  const double edge_angle_degrees = 7.5;
  ss_vect<double> actual =
    make_state(2.4e-3, -3.0e-4, -1.6e-3, 1.9e-4, 1.3e-2, 0.0);
  ss_vect<double> expected = actual;

  // Zero-gap hard edge: Delta px = h tan(phi) x, Delta py = -h tan(phi) y,
  // independent of delta in the default dip_edge_fudge mode.
  const double tangent = std::tan(radians(edge_angle_degrees));
  expected[px_] += curvature*tangent*expected[x_];
  expected[py_] -= curvature*tangent*expected[y_];

  EdgeFocus(curvature, edge_angle_degrees, 0.0, actual);
  expect_state("bend edge", expected, actual);
}

void test_corrector()
{
  const double horizontal_kick = 2.7e-4;
  const double vertical_kick = -1.9e-4;
  CellType cell = make_thin_multipole("corrector", Dip, horizontal_kick);
  ss_vect<double> actual =
    make_state(1.0e-3, -2.0e-4, -3.0e-3, 4.0e-4, 0.0, 0.0);
  ss_vect<double> expected = actual;

  cell.Elem.M->PB[HOMmax-Dip] = vertical_kick;

  // Integrated dipole fields give coordinate-independent canonical kicks.
  expected[px_] -= horizontal_kick;
  expected[py_] += vertical_kick;

  Mpole_Pass(cell, actual);
  expect_state("corrector", expected, actual);
  std::free(cell.Elem.M);
}

void test_cavity()
{
  CellType cell = make_cavity("cavity");
  ss_vect<double> actual =
    make_state(1.2e-3, 2.1e-4, -9.0e-4, -1.7e-4, 2.0e-3, 1.1e-2);
  ss_vect<double> expected = actual;

  // Two paraxial half drifts around the sinusoidal RF energy kick.
  globval.Cavity_on = true;
  apply_drift_oracle(cell.Elem.PL/2.0, expected);
  expected[delta_] -= cell.Elem.C->V_RF/(globval.Energy*1.0e9)
                      *std::sin(2.0*M_PI*cell.Elem.C->f_RF/c0*expected[ct_]
                                -cell.Elem.C->phi_RF);
  apply_drift_oracle(cell.Elem.PL/2.0, expected);

  Cav_Pass(cell, actual);
  expect_state("cavity", expected, actual);
  std::free(cell.Elem.C);
}

// --- 2D misalignment -------------------------------------------------------
//
// Direct GtoL/LtoG coordinate checks, a geometric round trip, displaced-
// quadrupole feed-down, rolled-quadrupole coupling, and a combined case, so
// that a frame-transform failure is separated from a field-kick failure.

void test_global_to_local_transform()
{
  const Vector2 displacement = {4.0e-4, -3.0e-4};
  const double angle = 0.23;
  const Vector2 rotation = {std::cos(angle), std::sin(angle)};
  CellType frame = {};
  const ss_vect<double> initial =
    make_state(2.1e-3, -7.0e-4, -1.4e-3, 9.0e-4, 8.0e-3, -6.0e-4);
  ss_vect<double> actual = initial, expected = initial;

  // Independent geometric transform: translate to the magnetic axis, then
  // rotate position and canonical momentum into the magnet frame.
  const double shifted_x = initial[x_]-displacement[X_];
  const double shifted_y = initial[y_]-displacement[Y_];
  expected[x_] = rotation[X_]*shifted_x+rotation[Y_]*shifted_y;
  expected[y_] = -rotation[Y_]*shifted_x+rotation[X_]*shifted_y;
  expected[px_] = rotation[X_]*initial[px_]+rotation[Y_]*initial[py_];
  expected[py_] = -rotation[Y_]*initial[px_]+rotation[X_]*initial[py_];

  set_misalignment(frame, displacement[X_], displacement[Y_], angle);
  to_local(actual, frame, 0.0, 0.0, 0.0);
  expect_state("GtoL off-center particle", expected, actual);
}

void test_local_to_global_transform()
{
  const Vector2 displacement = {-2.5e-4, 3.5e-4};
  const double angle = -0.19;
  const Vector2 rotation = {std::cos(angle), std::sin(angle)};
  CellType frame = {};
  const ss_vect<double> initial =
    make_state(-1.8e-3, 6.0e-4, 1.2e-3, -5.0e-4, -7.0e-3, 9.0e-4);
  ss_vect<double> actual = initial, expected = initial;

  // Independent inverse rotation followed by translation back to global axes.
  expected[x_] = rotation[X_]*initial[x_]-rotation[Y_]*initial[y_]
                 +displacement[X_];
  expected[y_] = rotation[Y_]*initial[x_]+rotation[X_]*initial[y_]
                 +displacement[Y_];
  expected[px_] = rotation[X_]*initial[px_]-rotation[Y_]*initial[py_];
  expected[py_] = rotation[Y_]*initial[px_]+rotation[X_]*initial[py_];

  set_misalignment(frame, displacement[X_], displacement[Y_], angle);
  to_global(actual, frame, 0.0, 0.0, 0.0);
  expect_state("LtoG off-center particle", expected, actual);
}

void test_geometric_transform_round_trip()
{
  CellType frame = {};
  const ss_vect<double> initial =
    make_state(3.0e-3, -1.1e-3, 2.0e-3, 8.0e-4, 1.4e-2, -1.0e-3);
  ss_vect<double> actual = initial;

  set_misalignment(frame, 5.0e-4, -7.0e-4, 0.41);
  to_local(actual, frame, 0.0, 0.0, 0.0);
  to_global(actual, frame, 0.0, 0.0, 0.0);
  expect_state("GtoL/LtoG geometric round trip", initial, actual);
}

void test_offset_quadrupole_feed_down()
{
  const double integrated_strength = 0.82;
  CellType cell = make_thin_multipole("offset quadrupole", Quad,
                                      integrated_strength);
  ss_vect<double> actual = make_state(0.0, 0.0, 0.0, 0.0, 5.0e-3, 2.0e-4);
  ss_vect<double> expected = actual;

  set_misalignment(cell, 4.0e-4, -3.0e-4, 0.0);

  // A particle on the design orbit is off axis in the displaced magnet.
  // Its local coordinates (-dx, -dy) generate dipole feed-down kicks.
  expected[px_] += integrated_strength*cell.dS[X_];
  expected[py_] -= integrated_strength*cell.dS[Y_];

  Mpole_Pass(cell, actual);
  expect_state("offset quadrupole feed-down", expected, actual);
  std::free(cell.Elem.M);
}

void test_rolled_quadrupole_coupling()
{
  const double integrated_strength = -0.67;
  const double roll = 0.31;
  const double cosine = std::cos(roll);
  const double sine = std::sin(roll);
  CellType cell = make_thin_multipole("rolled quadrupole", Quad,
                                      integrated_strength);
  const ss_vect<double> initial =
    make_state(1.7e-3, 0.0, -9.0e-4, 0.0, -3.0e-3, 0.0);
  ss_vect<double> actual = initial, expected = initial;

  set_misalignment(cell, 0.0, 0.0, roll);

  const double local_x = cosine*initial[x_]+sine*initial[y_];
  const double local_y = -sine*initial[x_]+cosine*initial[y_];
  const double local_px_kick = -integrated_strength*local_x;
  const double local_py_kick = integrated_strength*local_y;
  expected[px_] = cosine*local_px_kick-sine*local_py_kick;
  expected[py_] = sine*local_px_kick+cosine*local_py_kick;

  Mpole_Pass(cell, actual);
  expect_state("rolled quadrupole coupling", expected, actual);
  std::free(cell.Elem.M);
}

void test_misaligned_quadrupole()
{
  const double integrated_strength = 0.74;
  CellType cell = make_thin_multipole("misaligned quadrupole", Quad,
                                      integrated_strength);
  ss_vect<double> actual =
    make_state(1.9e-3, -2.2e-4, -1.3e-3, 3.1e-4, 4.0e-3, -7.0e-4);
  ss_vect<double> expected = actual, local;

  set_misalignment(cell, 3.2e-4, -2.6e-4, 0.17);

  // Independently compose translation, passive rotation, local quadrupole
  // kick, inverse rotation, and inverse translation.
  local = expected;
  const double shifted_x = expected[x_]-cell.dS[X_];
  const double shifted_y = expected[y_]-cell.dS[Y_];
  local[x_] = cell.dT[X_]*shifted_x+cell.dT[Y_]*shifted_y;
  local[y_] = -cell.dT[Y_]*shifted_x+cell.dT[X_]*shifted_y;
  local[px_] = cell.dT[X_]*expected[px_]+cell.dT[Y_]*expected[py_];
  local[py_] = -cell.dT[Y_]*expected[px_]+cell.dT[X_]*expected[py_];
  local[px_] -= integrated_strength*local[x_];
  local[py_] += integrated_strength*local[y_];
  expected[x_] = cell.dT[X_]*local[x_]-cell.dT[Y_]*local[y_]+cell.dS[X_];
  expected[y_] = cell.dT[Y_]*local[x_]+cell.dT[X_]*local[y_]+cell.dS[Y_];
  expected[px_] = cell.dT[X_]*local[px_]-cell.dT[Y_]*local[py_];
  expected[py_] = cell.dT[Y_]*local[px_]+cell.dT[X_]*local[py_];

  Mpole_Pass(cell, actual);
  expect_state("misaligned quadrupole", expected, actual);
  std::free(cell.Elem.M);
}

// --- fixed lattice ---------------------------------------------------------

void test_fixed_lattice()
{
  ss_vect<double> actual =
    make_state(1.4e-3, -3.0e-4, -8.0e-4, 2.2e-4, 6.0e-3, 3.0e-4);
  ss_vect<double> expected = actual;
  long last_position = -1;

  Cell[0] = make_drift_cell("D1", 0.3);
  Cell[1] = make_thin_multipole("QF", Quad, 0.6);
  Cell[2] = make_drift_cell("D2", 0.2);
  Cell[3] = make_thin_multipole("HC", Dip, 2.0e-4);

  apply_drift_oracle(0.3, expected);
  expected[px_] -= 0.6*expected[x_];
  expected[py_] += 0.6*expected[y_];
  apply_drift_oracle(0.2, expected);
  expected[px_] -= 2.0e-4;

  Cell_Pass(0, 3, actual, last_position);
  expect_index("fixed lattice survival", 3, last_position);
  expect_state("fixed lattice", expected, actual);

  std::free(Cell[1].Elem.M);
  std::free(Cell[3].Elem.M);
}

void test_fixed_lattice_transfer_map()
{
  const double drift_1 = 0.3;
  const double drift_2 = 0.2;
  const double quadrupole = 0.6;
  const double probe = 1.0e-6;
  const double expected[4][4] = {
    {1.0-quadrupole*drift_2,
     drift_1+drift_2-quadrupole*drift_1*drift_2, 0.0, 0.0},
    {-quadrupole, 1.0-quadrupole*drift_1, 0.0, 0.0},
    {0.0, 0.0, 1.0+quadrupole*drift_2,
     drift_1+drift_2+quadrupole*drift_1*drift_2},
    {0.0, 0.0, quadrupole, 1.0+quadrupole*drift_1}
  };
  char quantity[32];

  Cell[0] = make_drift_cell("D1", drift_1);
  Cell[1] = make_thin_multipole("QF", Quad, quadrupole);
  Cell[2] = make_drift_cell("D2", drift_2);

  // Symmetric double tracking extracts the linear map without relying on
  // Tracy's TPSA map construction. The oracle is the analytic D-Q-D product.
  for (int input = 0; input < 4; ++input) {
    ss_vect<double> plus, minus;
    long plus_last = -1, minus_last = -1;
    plus.zero();
    minus.zero();
    plus[input] = probe;
    minus[input] = -probe;
    Cell_Pass(0, 2, plus, plus_last);
    Cell_Pass(0, 2, minus, minus_last);
    expect_index("fixed lattice map positive probe", 2, plus_last);
    expect_index("fixed lattice map negative probe", 2, minus_last);
    for (int output = 0; output < 4; ++output) {
      std::snprintf(quantity, sizeof(quantity), "M[%s,%s]",
                    coordinate_names[output], coordinate_names[input]);
      expect_value("fixed lattice", quantity, expected[output][input],
                   (plus[output]-minus[output])/(2.0*probe));
    }
  }

  std::free(Cell[1].Elem.M);
}

void test_loss_location()
{
  ss_vect<double> actual = make_state(0.0, 2.0e-2, 0.0, 0.0, 0.0, 0.0);
  ss_vect<double> expected = actual;
  long last_position = -1;

  set_up_aperture_line(X_);
  apply_drift_oracle(0.2, expected);
  apply_drift_oracle(0.2, expected);

  // x reaches 8e-3 at D2, beyond its 6e-3 aperture.  Tracking must stop
  // there, before D3, and leave the state as it was at the loss.
  globval.Aperture_on = true;
  status.lossplane = 0;
  Cell_Pass(0, 3, actual, last_position);
  expect_index("fixed lattice loss", 2, last_position);
  expect_index("fixed lattice loss plane", 1, status.lossplane);
  expect_state("fixed lattice loss state", expected, actual);
}

} // namespace

void run_element_tests()
{
  run_test("drift", test_drift);
  run_test("quadrupole", test_quadrupole);
  run_test("quadrupole TPSA Jacobian and symplecticity",
           test_quadrupole_linear_map);
  run_test("sextupole", test_sextupole);
  run_test("thick quadrupole", test_thick_quadrupole);
  run_test("thick quadrupole off momentum",
           test_thick_quadrupole_off_momentum);
  run_test("thick quadrupole integrator order",
           test_thick_quadrupole_integrator_order);
  run_test("sector bend body", test_sector_bend_body);
  run_test("sector bend nonlinear order", test_sector_bend_nonlinear_order);
  run_test("sector bend slice convergence",
           test_sector_bend_slice_convergence);
  run_test("bend edge focusing", test_bend_edges);
  run_test("two-plane corrector", test_corrector);
  run_test("RF cavity", test_cavity);
  run_test("GtoL off-center particle", test_global_to_local_transform);
  run_test("LtoG off-center particle", test_local_to_global_transform);
  run_test("GtoL/LtoG geometric inverse", test_geometric_transform_round_trip);
  run_test("offset quadrupole feed-down", test_offset_quadrupole_feed_down);
  run_test("rolled quadrupole coupling", test_rolled_quadrupole_coupling);
  run_test("combined offset/roll quadrupole", test_misaligned_quadrupole);
  run_test("fixed lattice 6D state", test_fixed_lattice);
  run_test("fixed lattice transfer map", test_fixed_lattice_transfer_map);
  run_test("fixed lattice loss location", test_loss_location);
}

} // namespace tracking_test
