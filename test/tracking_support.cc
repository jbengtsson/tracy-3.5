#include "tracking_support.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

// Tracy applications select the active TPSA order at link time, so libtracy
// expects the application to define this symbol.
int no_tps = NO;

namespace tracking_test {

int failure_count = 0;

namespace {

bool both_finite(const double first, const double second)
{
  return std::isfinite(first) && std::isfinite(second);
}

bool measuring()
{
  static const bool enabled = std::getenv("TRACKING_TEST_MEASURE") != 0;
  return enabled;
}

// The conjugate of each coordinate and the sign of its Omega entry.
const int symplectic_conjugate[ss_dim] = {px_, x_, py_, y_, ct_, delta_};
const double symplectic_sign[ss_dim] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};

} // namespace

// --- assertions ------------------------------------------------------------

void expect_near(const char* element, const char* quantity,
                 const double expected, const double actual,
                 const double allowed)
{
  if (!both_finite(expected, actual) || !std::isfinite(allowed)
      || std::fabs(actual-expected) > allowed) {
    std::fprintf(stderr,
                 "%s %s: expected %.17g, actual %.17g, tolerance %.3g\n",
                 element, quantity, expected, actual, allowed);
    ++failure_count;
  }
}

void expect_value(const char* element, const char* quantity,
                  const double expected, const double actual)
{
  const double scale = std::max(1.0, std::max(std::fabs(expected),
                                              std::fabs(actual)));
  expect_near(element, quantity, expected, actual, tolerance*scale);
}

void expect_at_most(const char* element, const char* quantity,
                    const double bound, const double actual)
{
  if (!both_finite(bound, actual) || actual > bound) {
    std::fprintf(stderr, "%s %s: %.17g exceeds bound %.3g\n", element,
                 quantity, actual, bound);
    ++failure_count;
  }
}

void expect_at_least(const char* element, const char* quantity,
                     const double bound, const double actual)
{
  if (!both_finite(bound, actual) || actual < bound) {
    std::fprintf(stderr, "%s %s: %.17g is below bound %.3g\n", element,
                 quantity, actual, bound);
    ++failure_count;
  }
}

void expect_positive(const char* element, const char* quantity,
                     const double actual)
{
  // Guards the order checks: a residual that has collapsed to the roundoff
  // floor would make the measured exponent meaningless rather than wrong.
  if (!std::isfinite(actual) || actual <= 0.0) {
    std::fprintf(stderr, "%s %s: %.17g is not a usable positive residual\n",
                 element, quantity, actual);
    ++failure_count;
  }
}

void expect_index(const char* test, const long expected, const long actual)
{
  if (actual != expected) {
    std::fprintf(stderr, "%s location: expected %ld, actual %ld\n",
                 test, expected, actual);
    ++failure_count;
  }
}

void expect_state(const char* element, const ss_vect<double>& expected,
                  const ss_vect<double>& actual)
{
  for (int coordinate = 0; coordinate < ss_dim; ++coordinate)
    expect_value(element, coordinate_names[coordinate], expected[coordinate],
                 actual[coordinate]);
}

void expect_matrix(const char* element, const Matrix& expected,
                   const Matrix& actual)
{
  char quantity[32];

  for (int row = 0; row < ss_dim; ++row)
    for (int column = 0; column < ss_dim; ++column) {
      std::snprintf(quantity, sizeof(quantity), "M[%s,%s]",
                    coordinate_names[row], coordinate_names[column]);
      expect_value(element, quantity, expected[row][column],
                   actual[row][column]);
    }
}

void expect_symplectic(const char* element, const Matrix& map)
{
  char quantity[48];

  // Check M^T Omega M = Omega directly so this invariant stays independent of
  // Tracy matrix utilities.
  for (int row = 0; row < ss_dim; ++row)
    for (int column = 0; column < ss_dim; ++column) {
      std::snprintf(quantity, sizeof(quantity), "M^T Omega M[%s,%s]",
                    coordinate_names[row], coordinate_names[column]);
      expect_value(element, quantity, symplectic_identity(row, column),
                   symplectic_product(map, row, column));
    }
}

void report(const char* test, const char* quantity, const double value)
{
  if (measuring())
    std::printf("  measured %s %s: %.6e\n", test, quantity, value);
}

// --- residuals -------------------------------------------------------------

double larger(const double largest, const double value)
{
  if (std::isnan(largest) || std::isnan(value))
    return std::numeric_limits<double>::quiet_NaN();
  return std::max(largest, value);
}

double max_deviation(const ss_vect<double>& left, const ss_vect<double>& right)
{
  double largest = 0.0;

  for (int coordinate = 0; coordinate < ss_dim; ++coordinate)
    largest = larger(largest, std::fabs(left[coordinate]-right[coordinate]));
  return largest;
}

double symplectic_product(const Matrix& map, const int row, const int column)
{
  double product = 0.0;

  for (int intermediate = 0; intermediate < ss_dim; ++intermediate)
    product += symplectic_sign[intermediate]*map[intermediate][row]
               *map[symplectic_conjugate[intermediate]][column];
  return product;
}

double symplectic_identity(const int row, const int column)
{
  return column == symplectic_conjugate[row] ? symplectic_sign[row] : 0.0;
}

// --- configuration and oracles ---------------------------------------------

void configure_tracking()
{
  globval.Cavity_on = false;
  globval.radiation = false;
  globval.emittance = false;
  globval.quad_fringe = false;
  globval.H_exact = false;
  globval.Cart_Bend = false;
  globval.pathlength = false;
  globval.mat_meth = false;
  globval.Aperture_on = false;
  globval.dip_edge_fudge = true;
  globval.Energy = 3.0;
  SI_init();
}

ss_vect<double> make_state(const double x, const double px, const double y,
                           const double py, const double delta,
                           const double ct)
{
  ss_vect<double> state;

  state[x_] = x;
  state[px_] = px;
  state[y_] = y;
  state[py_] = py;
  state[delta_] = delta;
  state[ct_] = ct;
  return state;
}

void apply_drift_oracle(const double length, ss_vect<double>& state)
{
  const double relative_momentum = 1.0+state[delta_];
  state[x_] += length*state[px_]/relative_momentum;
  state[y_] += length*state[py_]/relative_momentum;
  state[ct_] += length*(state[px_]*state[px_]+state[py_]*state[py_])
                /(2.0*relative_momentum*relative_momentum);
}

double radians(const double degrees) { return degrees*M_PI/180.0; }

// --- misalignment ----------------------------------------------------------

void set_misalignment(CellType& cell, const double dx, const double dy,
                      const double roll)
{
  cell.dS[X_] = dx;
  cell.dS[Y_] = dy;
  cell.dT[X_] = std::cos(roll);
  cell.dT[Y_] = std::sin(roll);
}

// --- fixtures --------------------------------------------------------------

CellType make_drift_cell(const char* name, const double length)
{
  CellType cell = {};

  std::strncpy(cell.Elem.PName, name, sizeof(cell.Elem.PName)-1);
  cell.Elem.Pkind = drift;
  cell.Elem.PL = length;
  set_misalignment(cell, 0.0, 0.0, 0.0);
  return cell;
}

CellType make_thin_multipole(const char* name, const int order,
                             const double strength)
{
  CellType cell = {};

  std::strncpy(cell.Elem.PName, name, sizeof(cell.Elem.PName)-1);
  cell.Elem.Pkind = Mpole;
  cell.Elem.PL = 0.0;
  set_misalignment(cell, 0.0, 0.0, 0.0);
  Mpole_Alloc(&cell.Elem);
  cell.Elem.M->Pthick = thin;
  cell.Elem.M->Porder = order;
  cell.Elem.M->n_design = order;
  cell.Elem.M->PB[HOMmax+order] = strength;
  return cell;
}

CellType make_thick_quadrupole(const char* name, const double length,
                               const double gradient, const int slices)
{
  CellType cell = make_thin_multipole(name, Quad, gradient);

  cell.Elem.PL = length;
  cell.Elem.M->Pthick = thick;
  cell.Elem.M->PN = slices;
  return cell;
}

CellType make_bend(const char* name, const double length,
                   const double curvature, const double entrance_degrees,
                   const double exit_degrees, const int slices)
{
  const double angle = length*curvature;
  CellType cell = make_thin_multipole(name, 0, 0.0);

  cell.Elem.PL = length;
  cell.Elem.M->Pthick = thick;
  cell.Elem.M->PN = slices;
  cell.Elem.M->Pirho = curvature;
  cell.Elem.M->PTx1 = entrance_degrees;
  cell.Elem.M->PTx2 = exit_degrees;
  cell.Elem.M->Pgap = 0.0;
  cell.Elem.M->Pc0 = std::sin(angle/2.0);
  cell.Elem.M->Pc1 = cell.Elem.M->Pc0;
  return cell;
}

CellType make_cavity(const char* name)
{
  CellType cell = {};

  std::strncpy(cell.Elem.PName, name, sizeof(cell.Elem.PName)-1);
  cell.Elem.Pkind = Cavity;
  cell.Elem.PL = 0.4;
  Cav_Alloc(&cell.Elem);
  cell.Elem.C->V_RF = 2.0e6;
  cell.Elem.C->f_RF = 500.0e6;
  cell.Elem.C->phi_RF = 0.35;
  return cell;
}

void set_up_aperture_line(const int plane)
{
  Cell[0] = make_drift_cell("START", 0.0);
  Cell[1] = make_drift_cell("D1", 0.2);
  Cell[2] = make_drift_cell("D2", 0.2);
  Cell[3] = make_drift_cell("D3", 0.2);
  for (int i = 0; i <= 3; ++i) {
    Cell[i].maxampl[X_][0] = -1.0e-2;
    Cell[i].maxampl[X_][1] = i == 2 && plane == X_ ? 6.0e-3 : 1.0e-2;
    Cell[i].maxampl[Y_][0] = -1.0e-2;
    Cell[i].maxampl[Y_][1] = i == 2 && plane == Y_ ? 6.0e-3 : 1.0e-2;
  }
}

// --- runners ---------------------------------------------------------------

void run_test(const char* name, void (*test)())
{
  const int failures_before = failure_count;
  configure_tracking();
  test();
  std::printf("%s: %s\n", failure_count == failures_before ? "PASS" : "FAIL",
              name);
}

} // namespace tracking_test
