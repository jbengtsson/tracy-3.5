#define NO 1

#include "tracy_lib.h"

int no_tps = NO;


void compute_fel_dnu(const string &marker)
{
  const std::vector<string>
    oct_name = {"oct_1", "oct_2", "oct_3", "oct_4", "oct_5",
		"oct_5", "oct_4", "oct_3", "oct_2", "oct_1"};
  const std::vector<double>
    oct_kid  = {   1,        1,       1,       1,      2,
		   2,        1,       1,       1,      1};


  std::vector<int>                   loc;
  std::vector< std::vector<double> > dnu;
  std::vector<double>                pair;

  for (auto k = 0; k < oct_name.size(); k++)
    loc.push_back(Elem_GetPos(ElemIndex(oct_name[k].c_str()), oct_kid[k]));
  for (auto k = 0; k < oct_name.size(); k++) {
    pair.clear();
    auto cell = Cell[loc[k]];
    pair.push_back(cell.Nu[X_]);
    pair.push_back(cell.Nu[Y_]);
    dnu.push_back(pair);
  }

  printf("\n   loc  name       dnu_x     dnu_y\n");
  printf("   %3d  %.8s  %6.3f    %6.3f\n",
	 loc[0], Cell[loc[0]].Elem.PName, 0e0, 0e0);
  for (auto k = 1; k < oct_name.size(); k++) {
    printf("   %3d  %.8s  %6.3f    %6.3f\n",
	   loc[k], Cell[loc[k]].Elem.PName, dnu[k][X_]-dnu[k-1][X_],
	   dnu[k][Y_]-dnu[k-1][Y_]);
  }
}


void get_kick_map
(const string &file_name, const int nx, const int ny, const long int loc1,
 const long int loc2, const double Ax, const double Ay)
{
  const double Brho = globval.Energy*1e9/c0;

  long int        lastpos;
  ss_vect<double> ps0, ps1, dps_map;
  ofstream        outf;

  cout << "\nget_kick_map:\n  " << setw(8) << Cell[loc1].Elem.PName << " -> "
       << setw(8) << Cell[loc2].Elem.PName << "\n";

  file_wr(outf, file_name.c_str());

  dps_map[x_] = 2e0*Ax/(nx-1e0);
  dps_map[y_] = 2e0*Ay/(ny-1e0);

  outf << "# Author:" << "\n";
  outf << "# Title" << "\n";
  outf << "# Cell Length [m]" << "\n";
  outf << fixed << setprecision(5) << Cell[loc2].S-Cell[loc1].S
       << "\n";
  outf << "# Number of Horizontal Points" << "\n";
  outf << nx << "\n";
  outf << "# Number of Vertical Points" << "\n";
  outf << ny << "\n";

  outf << "# Horizontal 2nd Order Kick [T2m2]" << "\n";
  outf << "START" << "\n";

  for (auto i1 = 0; i1 < nx; i1++)
    outf << scientific << setprecision(5) << setw(13)
	 << -Ax+i1*dps_map[x_];
  outf << "\n";
  cout << "  scanning horizontal plane:\n    ";
  ps0.zero();
  for (auto i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];
    cout << ".";
    outf << scientific << setprecision(5)
	 << setw(13) << ps0[y_];
    for (auto i2 = 0; i2 < nx; i2++) {
    ps0[x_] = -Ax + i2*dps_map[x_];
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      if (lastpos == loc2) {
	ps1 -= ps0;
	outf << scientific << setprecision(5)
	     << setw(13) << sqr(Brho)*ps1[px_];
      } else
	outf << scientific << setprecision(5) << setw(13) << NAN;
    }
    outf << "\n";
  }
  cout << "\n";

  outf << "# Vertical 2nd Order Kick [T2m2]" << "\n";
  outf << "START" << "\n";

  for (auto i1 = 0; i1 < nx; i1++)
    outf << scientific << setprecision(5) << setw(13)
	 << -Ax+i1*dps_map[x_];
  outf << "\n";

  cout << "  scanning vertical plane:\n    ";
  ps0.zero();
  for (auto i1 = 0; i1 < ny; i1++) {
    ps0[y_] = Ay - i1*dps_map[y_];
    cout << ".";
    outf << scientific << setprecision(5) << setw(13) << ps0[y_];

    for (auto i2 = 0; i2 < nx; i2++) {
      ps0[x_] = -Ax + i2*dps_map[x_];
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      if (lastpos == loc2) {
	ps1 -= ps0;
	outf << scientific << setprecision(5)
	     << setw(13) << sqr(Brho)*ps1[py_];
      } else
	outf << scientific << setprecision(5) << setw(13) << NAN;
    }
    outf << "\n";
  }
  cout << "\n";

  outf.close();
}


void get_map_2D
(const string &file_name, const int nx, const int ny, const long int loc1,
 const long int loc2, const double Ax, const double Ay)
{
  long int        lastpos;
  ss_vect<double> ps0, ps1, dps_map;
  ofstream        outf;

  cout << "\nget_map_2D:\n  " << setw(8) << Cell[loc1].Elem.PName << " -> "
       << setw(8) << Cell[loc2].Elem.PName << "\n";

  file_wr(outf, file_name.c_str());

  dps_map[x_] = 2e0*Ax/(nx-1e0);
  dps_map[y_] = 2e0*Ay/(ny-1e0);

  outf << "# Units are [rad].\n";
  outf << scientific << setprecision(5)
       << "# nx = " << nx << "     Ax = " << Ax
       << "     dx = " << dps_map[x_] << "\n";
  outf << scientific << setprecision(5)
       << "# ny = " << ny << "     Ay = " << Ay
       << "     dy = " << dps_map[y_] << "\n";
  outf << "#" << "\n";

  // cout << "\n";
  ps0.zero();
  for (auto i1 = 0; i1 < nx; i1++) {
    ps0[x_] = -Ax + i1*dps_map[x_];
    for (auto i2 = 0; i2 < ny; i2++) {
      ps0[y_] = -Ay + i2*dps_map[y_];
//       cout << setw(3) << i1+1 << setw(3) << i2+1
// 	   << scientific << setprecision(5)
// 	   << setw(13) << ps0[x_] << setw(13) << ps0[y_] << "\n";
      ps1 = ps0;
      Cell_Pass(loc1, loc2, ps1, lastpos);
      ps1 -= ps0;
      outf << scientific << setprecision(5)
	   << setw(13) << ps0[x_] << setw(13) << ps0[y_]
	   << setw(13) << ps1[px_] << setw(13) << ps1[py_] << "\n";
    }
    outf << "\n";
  }
  outf.close();
}


struct ConversionStats {
  std::size_t dataRows = 0;
  std::size_t separatorsInserted = 0;
  std::size_t ascendingXRows = 0;
  std::size_t descendingXRows = 0;
  std::size_t undeterminedXRows = 0;
};

namespace {

  enum class XDirection { Unknown, Increasing, Decreasing };

  bool isBlank(const std::string& line) {
    return std::all_of(line.begin(), line.end(), [](unsigned char ch) {
      return std::isspace(ch) != 0;
    });
  }

  bool isComment(const std::string& line) {
    const auto first =
      std::find_if_not(line.begin(), line.end(), [](unsigned char ch) {
      return std::isspace(ch) != 0;
    });
    return first != line.end() && *first == '#';
  }

  bool nearlyEqual(double a, double b,
		   double absoluteTolerance = 1.0e-15,
		   double relativeTolerance = 1.0e-12) {
    const double scale = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) <= absoluteTolerance + relativeTolerance * scale;
  }

  void countDirection(XDirection direction, ConversionStats& stats) {
    switch (direction) {
    case XDirection::Increasing:
      ++stats.ascendingXRows;
      break;
    case XDirection::Decreasing:
      ++stats.descendingXRows;
      break;
    case XDirection::Unknown:
      ++stats.undeterminedXRows;
      break;
    }
  }

}  // namespace

// Copies the input table to outputPath and inserts one blank line immediately
// before every data row whose y value differs from the preceding data row.
//
// The first two columns must be x and y. Remaining columns are copied unchanged.
// Lines whose first non-whitespace character is '#' are copied as comments.
// The x sweep direction is detected independently for every constant-y row, so
// increasing, decreasing, and alternating (serpentine) scans are supported.
ConversionStats insertBlankLinesAtYChanges
(const std::string& inputPath, const std::string& outputPath,
 double absoluteTolerance = 1.0e-15, double relativeTolerance = 1.0e-12) {

  std::ifstream input(inputPath);
  if (!input) {
    throw std::runtime_error("Cannot open input file: " + inputPath);
  }

  std::ofstream output(outputPath, std::ios::trunc);
  if (!output) {
    throw std::runtime_error("Cannot open output file: " + outputPath);
  }

  ConversionStats stats;
  std::string line;
  std::size_t physicalLine = 0;

  bool havePreviousData = false;
  bool outputEndsWithBlankLine = false;
  double previousX = 0.0;
  double previousY = 0.0;
  XDirection currentDirection = XDirection::Unknown;

  while (std::getline(input, line)) {
    ++physicalLine;

    // std::getline removes '\n' but leaves '\r' for CRLF input.
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }

    if (isBlank(line)) {
      output << '\n';
      outputEndsWithBlankLine = true;
      continue;
    }

    if (isComment(line)) {
      output << line << '\n';
      outputEndsWithBlankLine = false;
      continue;
    }

    std::istringstream parser(line);
    double x = 0.0;
    double y = 0.0;
    if (!(parser >> x >> y)) {
      throw std::runtime_error
	("Expected numeric x and y in the first two columns at line " +
	 std::to_string(physicalLine));
    }
    if (!std::isfinite(x) || !std::isfinite(y)) {
      throw std::runtime_error
	("Non-finite x or y value at line " + std::to_string(physicalLine));
    }

    const bool yChanged =
      havePreviousData &&
      !nearlyEqual(y, previousY, absoluteTolerance, relativeTolerance);

    if (yChanged) {
      countDirection(currentDirection, stats);
      currentDirection = XDirection::Unknown;

      if (!outputEndsWithBlankLine) {
	output << '\n';
	++stats.separatorsInserted;
      }
    } else if (havePreviousData) {
      const bool xChanged =
	!nearlyEqual(x, previousX, absoluteTolerance, relativeTolerance);

      if (xChanged) {
	const XDirection stepDirection =
	  x > previousX ? XDirection::Increasing : XDirection::Decreasing;

	if (currentDirection == XDirection::Unknown) {
	  currentDirection = stepDirection;
	} else if (currentDirection != stepDirection) {
	  throw std::runtime_error
	    ( "x changes direction while y is constant at line " +
	      std::to_string(physicalLine));
	}
      }
    }

    // Preserve the original data text and all columns exactly.
    output << line << '\n';
    outputEndsWithBlankLine = false;

    previousX = x;
    previousY = y;
    havePreviousData = true;
    ++stats.dataRows;
  }

  if (input.bad()) {
    throw std::runtime_error("I/O error while reading: " + inputPath);
  }
  if (!output) {
    throw std::runtime_error("I/O error while writing: " + outputPath);
  }

  if (havePreviousData) {
    countDirection(currentDirection, stats);
  }

  return stats;
}


namespace kickmap_converter {

  struct SourcePoint {
    double x;
    double y;
    double horizontalKick;  // xpFactor
    double verticalKick;    // ypFactor
  };

  struct KickMapMetadata {
    double undulatorLengthMetres = 0e0;
    std::string title = "claris_und_kickmap";
    std::string author = "Generated from an SDDS-style kick-map table";
  };

  namespace {

    std::string trim(const std::string& text) {
      const std::size_t first = text.find_first_not_of(" \t\r\n");
      if (first == std::string::npos) {
        return {};
      }

      const std::size_t last = text.find_last_not_of(" \t\r\n");
      return text.substr(first, last - first + 1);
    }

    bool coordinateEqual(double a,
			 double b,
			 double absoluteTolerance = 1.0e-15,
			 double relativeTolerance = 1.0e-12) {
      const double scale = std::max(std::abs(a), std::abs(b));
      return std::abs(a - b) <=
	absoluteTolerance + relativeTolerance * scale;
    }

    std::string formatScientific(double value) {
      if (value == 0.0) {
        value = 0.0;  // Avoid emitting "-0.00000".
      }

      std::ostringstream stream;
      stream << std::scientific << std::setprecision(5) << value;
      std::string result = stream.str();

      const std::size_t exponentPosition = result.find('e');
      if (exponentPosition == std::string::npos ||
	  exponentPosition + 2 >= result.size()) {
        throw std::runtime_error("Could not format a floating-point value.");
      }

      const char sign = result.at(exponentPosition + 1);
      std::string exponentDigits = result.substr(exponentPosition + 2);

      while (exponentDigits.size() < 3) {
        exponentDigits.insert(exponentDigits.begin(), '0');
      }

      return result.substr(0, exponentPosition + 1) +
	sign + exponentDigits;
    }

    std::vector<SourcePoint> readSourcePoints(const std::string& inputPath) {
      std::ifstream input(inputPath);
      if (!input) {
        throw std::runtime_error("Cannot open input file: " + inputPath);
      }

      std::vector<SourcePoint> points;
      std::string line;
      std::size_t lineNumber = 0;

      while (std::getline(input, line)) {
        ++lineNumber;
        const std::string stripped = trim(line);

        if (stripped.empty() || stripped.front() == '#' ||
            stripped.front() == '&') {
	  continue;
        }

        std::istringstream parser(stripped);
        SourcePoint point{};

        if (parser >> point.x >> point.y >>
	    point.horizontalKick >> point.verticalKick) {
	  std::string extraToken;
	  if (parser >> extraToken) {
	    throw std::runtime_error(
				     "Unexpected extra column at line " +
				     std::to_string(lineNumber));
	  }

	  if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
	      !std::isfinite(point.horizontalKick) ||
	      !std::isfinite(point.verticalKick)) {
	    throw std::runtime_error(
				     "Non-finite value at line " +
				     std::to_string(lineNumber));
	  }

	  points.push_back(point);
	  continue;
        }

        // SDDS ASCII files may contain a single row-count value.
        parser.clear();
        parser.str(stripped);
        double rowCount = 0.0;
        std::string extraToken;
        if ((parser >> rowCount) && !(parser >> extraToken)) {
	  continue;
        }

        throw std::runtime_error(
				 "Cannot parse line " + std::to_string(lineNumber) +
				 ": " + stripped);
      }

      if (points.empty()) {
        throw std::runtime_error("No four-column data rows were found.");
      }

      return points;
    }

    std::vector<double> uniqueCoordinates(
					  const std::vector<SourcePoint>& points,
					  bool useXCoordinate) {

      std::vector<double> coordinates;
      coordinates.reserve(points.size());

      for (const SourcePoint& point : points) {
        coordinates.push_back(useXCoordinate ? point.x : point.y);
      }

      std::sort(coordinates.begin(), coordinates.end());

      std::vector<double> unique;
      unique.reserve(coordinates.size());

      for (double value : coordinates) {
        if (unique.empty() ||
            !coordinateEqual(unique.back(), value)) {
	  unique.push_back(value);
        }
      }

      return unique;
    }

    std::size_t findCoordinateIndex(const std::vector<double>& coordinates,
				    double value) {
      const auto position =
        std::lower_bound(coordinates.begin(), coordinates.end(), value);

      if (position != coordinates.end() &&
	  coordinateEqual(*position, value)) {
        return static_cast<std::size_t>(
					std::distance(coordinates.begin(), position));
      }

      if (position != coordinates.begin()) {
        const auto previous = std::prev(position);
        if (coordinateEqual(*previous, value)) {
	  return static_cast<std::size_t>(
					  std::distance(coordinates.begin(), previous));
        }
      }

      throw std::runtime_error(
			       "A data coordinate could not be matched to the grid.");
    }

    struct RectangularGrid {
      std::vector<double> xCoordinates;  // Ascending.
      std::vector<double> yCoordinates;  // Ascending internally.
      std::vector<std::vector<double>> horizontal;
      std::vector<std::vector<double>> vertical;
    };

    RectangularGrid makeRectangularGrid(
					const std::vector<SourcePoint>& points) {

      RectangularGrid grid;
      grid.xCoordinates = uniqueCoordinates(points, true);
      grid.yCoordinates = uniqueCoordinates(points, false);

      const std::size_t nx = grid.xCoordinates.size();
      const std::size_t ny = grid.yCoordinates.size();

      if (nx == 0 || ny == 0 ||
	  nx > std::numeric_limits<std::size_t>::max() / ny ||
	  points.size() != nx * ny) {
        throw std::runtime_error(
				 "The input does not contain one complete rectangular grid.");
      }

      grid.horizontal.assign(
			     ny, std::vector<double>(nx, 0.0));
      grid.vertical.assign(
			   ny, std::vector<double>(nx, 0.0));

      std::vector<std::vector<bool>> occupied(
					      ny, std::vector<bool>(nx, false));

      for (const SourcePoint& point : points) {
        const std::size_t ix =
	  findCoordinateIndex(grid.xCoordinates, point.x);
        const std::size_t iy =
	  findCoordinateIndex(grid.yCoordinates, point.y);

        if (occupied[iy][ix]) {
	  throw std::runtime_error(
				   "Duplicate point detected at x=" +
				   std::to_string(point.x) + ", y=" +
				   std::to_string(point.y));
        }

        occupied[iy][ix] = true;
        grid.horizontal[iy][ix] = point.horizontalKick;
        grid.vertical[iy][ix] = point.verticalKick;
      }

      for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
	  if (!occupied[iy][ix]) {
	    throw std::runtime_error(
				     "Missing point in the rectangular grid.");
	  }
        }
      }

      return grid;
    }

    void writeCoordinateHeader(
			       std::ostream& output,
			       const std::vector<double>& xCoordinates) {

      output << std::string(15, ' ');
      for (double x : xCoordinates) {
        output << std::setw(15) << formatScientific(x);
      }
      output << '\n';
    }

    void writeDataset(
		      std::ostream& output,
		      const std::string& heading,
		      const RectangularGrid& grid,
		      const std::vector<std::vector<double>>& values) {

      output << "# " << heading << '\n';
      output << "START\n";
      writeCoordinateHeader(output, grid.xCoordinates);

      // The reference format lists y from largest to smallest.
      for (std::size_t reverseIndex = grid.yCoordinates.size();
	   reverseIndex > 0;
	   --reverseIndex) {
        const std::size_t iy = reverseIndex - 1;

        output << std::setw(15)
               << formatScientific(grid.yCoordinates[iy]);

        for (double value : values[iy]) {
	  output << std::setw(15)
		 << formatScientific(value);
        }
        output << '\n';
      }
    }

  }  // namespace

  void writeKickMapWithoutThirdDataset(
				       const std::string& inputPath,
				       const std::string& outputPath,
				       const KickMapMetadata& metadata = {}) {

    const std::vector<SourcePoint> points =
      readSourcePoints(inputPath);
    const RectangularGrid grid =
      makeRectangularGrid(points);

    std::ofstream output(outputPath);
    if (!output) {
      throw std::runtime_error(
			       "Cannot open output file: " + outputPath);
    }

    if (!metadata.author.empty()) {
      output << "# Author : " << metadata.author << '\n';
    }
    if (!metadata.title.empty()) {
      output << "# " << metadata.title << '\n';
    }

    output << "# Undulator Length [m]\n";
    output << formatScientific(
			       metadata.undulatorLengthMetres) << '\n';

    output << "# Number of Horizontal Points\n";
    output << grid.xCoordinates.size() << '\n';

    output << "# Number of Vertical Points\n";
    output << grid.yCoordinates.size() << '\n';

    writeDataset(
		 output,
		 "Horizontal 2nd Order Kick [T2m2]",
		 grid,
		 grid.horizontal);

    writeDataset(
		 output,
		 "Vertical 2nd Order Kick [T2m2]",
		 grid,
		 grid.vertical);

    // Deliberately stop here. No third heading, START, or dataset is written.

    if (!output) {
      throw std::runtime_error(
			       "An error occurred while writing: " + outputPath);
    }
  }

}  // namespace kickmap_converter


void set_state(void)
{
  globval.H_exact        = false;
  globval.quad_fringe    = false;
  globval.Cavity_on      = false;
  globval.radiation      = false;
  globval.emittance      = false;
  globval.IBS            = false;
  globval.pathlength     = false;
  globval.Aperture_on    = false;
  globval.Cart_Bend      = false;
  globval.dip_edge_fudge = true;
  globval.EPU            = true;
}


int main(int argc, char *argv[])
{

  globval.mat_meth = false;

  FieldMap_filetype = 6;

  if (false) {
    const ConversionStats stats =
      insertBlankLinesAtYChanges(argv[1], argv[2]);

    std::cout << "Data rows: " << stats.dataRows << '\n'
	      << "Blank lines inserted: " << stats.separatorsInserted
	      << '\n'
	      << "Rows with increasing x: " << stats.ascendingXRows << '\n'
	      << "Rows with decreasing x: " << stats.descendingXRows << '\n'
	      << "Rows with undetermined x direction: "
	      << stats.undeterminedXRows << '\n';
    exit(0);
  }


  if (false) {
    kickmap_converter::KickMapMetadata metadata;

    kickmap_converter::writeKickMapWithoutThirdDataset
      (argv[1], argv[2], metadata);
    return 0;
  }

  if (true)
    Read_Lattice(argv[1]);
  else
    rdmfile(argv[1]);

  set_state();

  no_sxt();

  prtmfile("flat_file.dat");

  if (!false) {
    Ring_GetTwiss(true, 0e0);
    printglob();

    prt_lat("linlat1.out", globval.bpm, true);
    prt_lat("linlat.out", globval.bpm, true, 10);
  }

  if (false)
    GetEmittance(ElemIndex("cav"), false, true);

  if (false)
    compute_fel_dnu("fel_mark");

  if (!false) {
    const int
      n[] = {73, 73}, 
      // i_0 = 19;
      i_0 = 107;
    std::vector<double> A;

    // epu57v2lvg16kickmap2pure.dat
    // A = {28e-3, 5e-3};
    A = {9e-3, 9e-3};
    get_kick_map
      ("helical_und_km.dat", n[X_], n[Y_], i_0, i_0, A[X_], A[Y_]);

    // A = {3e-3, 3e-3};
    get_map_2D
      ("helical_und_2D.dat", n[X_], n[Y_], i_0, i_0, A[X_], A[Y_]);
  }
}
