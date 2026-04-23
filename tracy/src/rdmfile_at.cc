
// Requires C++11 but no later.

/*
  Passmethods used at MAX IV:
  ('*' implemented)

    * AperturePass
    BendLinearPass
    BndMPoleSymplectic4Pass
    * BndMPoleSymplectic4RadPass
    CavityPass
    * CorrectorPass
    * DriftPass
    * EAperturePass
    * IdTablePass
    * IdTableRadPass
    * IdentityPass
    * RFCavityPass
    * StrMPoleSymplectic4Pass
    * StrMPoleSymplectic4RadPass
    ThinMPolePass
    QuadLinearPass
                                                                              */


#include <array>
#include <sstream>
#include <algorithm>

static const bool dbg = false;

void string_to_c_str(const std::string &str, partsName &c_str) {
  // Tracy-2 element names are not "\0" terminated C strings (Pascal legacy).
  // Keep symbol names in the same canonical format used by ElemIndex:
  // lowercase pad with spaces up to SymbolLength.
  if (str.size() > NameLength)
    throw std::runtime_error("ElemName too long for fixed buffer");

  memset(c_str, 0, sizeof(partsName));
  for (size_t i = 0; i < str.size(); i++)
    c_str[i] = (char)std::tolower((unsigned char)str[i]);
  for (size_t i = str.size(); i < SymbolLength; i++)
    c_str[i] = ' ';
}

struct Value {
  bool isNumber = false;
  double number = 0.0;
  std::string text;

  static Value Num(double x)
  { Value v; v.isNumber = true; v.number = x; return v; }
  static Value Str(std::string s)
  { Value v; v.isNumber = false; v.text = std::move(s); return v; }
};

struct Element {
  std::string name;
  int number = -1;
  std::string passMethod;
  // Element properties dictionary.
  std::unordered_map<std::string, std::vector<Value>> props;
  bool hasT1 = false;
  bool hasT2 = false;
  bool hasR1 = false;
  bool hasR2 = false;
  std::array<double, 6> T1;
  std::array<double, 6> T2;
  std::array<std::array<double, 6>, 6> R1;
  std::array<std::array<double, 6>, 6> R2;
};

bool try_get_length(const Element& e, double& outLength)
{
  auto it = e.props.find("Length");
  if (it == e.props.end() || it->second.empty())
    return false;

  const Value& v = it->second.front();
  if (!v.isNumber)
    throw std::runtime_error("Element '" + e.name + "': Length is not numeric");

  outLength = v.number;
  return true;
}

double get_length(const Element& e)
{
  double L = 0.0;
  if (!try_get_length(e, L))
    throw std::runtime_error("Element '" + e.name + "' missing Length");
  return L;
}

static inline void ltrim_inplace(std::string& s) {
  size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  s.erase(0, i);
}
static inline void rtrim_inplace(std::string& s) {
  size_t i = s.size();
  while (i > 0 && std::isspace(static_cast<unsigned char>(s[i - 1]))) --i;
  s.erase(i);
}
static inline std::string trim(std::string s)
{ ltrim_inplace(s); rtrim_inplace(s); return s; }

static inline bool
starts_with(const std::string& s, const std::string& prefix) {
  return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static std::vector<std::string>
split_and_trim(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == delim) { out.push_back(trim(cur)); cur.clear(); }
    else cur.push_back(c);
  }
  out.push_back(trim(cur));
  return out;
}

static std::vector<std::string> split_csv_like(const std::string& s) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == ',') { out.push_back(trim(cur)); cur.clear(); }
    else cur.push_back(c);
  }
  if (!cur.empty() || !out.empty()) out.push_back(trim(cur));
  return out;
}

static bool try_parse_double(const std::string& token, double& out) {
  std::string t = trim(token);
  if (t.empty()) return false;

  const char* begin = t.c_str();
  char* end = nullptr;
  errno = 0;
  double val = std::strtod(begin, &end);

  if (begin == end) return false;
  while (end && *end && std::isspace(static_cast<unsigned char>(*end))) ++end;
  if (end && *end != '\0') return false;
  if (errno == ERANGE) return false;

  out = val;
  return true;
}

static double parse_strict_double(const std::string& token) {
  double out = 0.0;
  if (!try_parse_double(token, out))
    throw std::runtime_error("Invalid numeric token: '" + trim(token) + "'");
  return out;
}

static std::vector<double> parse_number_list(const std::string& rhs) {
  std::vector<double> vals;
  for (const auto& tok : split_csv_like(rhs)) {
    if (!tok.empty())
      vals.push_back(parse_strict_double(tok));
  }
  return vals;
}

static std::array<double, 6> parse_vec6(const std::string& rhs,
                                        const std::string& key) {
  const auto vals = parse_number_list(rhs);
  if (vals.size() != 6)
    throw std::runtime_error(key + " expects 6 values, got "
                             + std::to_string(vals.size()));

  std::array<double, 6> out{};
  std::copy(vals.begin(), vals.end(), out.begin());
  return out;
}

static std::array<std::array<double, 6>, 6>
parse_mat6x6(const std::string& rhs, const std::string& key) {
  const auto rows = split_and_trim(rhs, ';');

  std::vector<std::string> nonempty_rows;
  for (const auto& row : rows) {
    if (!row.empty())
      nonempty_rows.push_back(row);
  }

  if (nonempty_rows.size() != 6)
    throw std::runtime_error(key + " expects 6 rows, got "
                             + std::to_string(nonempty_rows.size()));

  std::array<std::array<double, 6>, 6> out{};
  for (size_t i = 0; i < 6; ++i) {
    const auto cols = split_csv_like(nonempty_rows[i]);
    if (cols.size() != 6)
      throw std::runtime_error(key + " row " + std::to_string(i+1)
                               + " expects 6 values, got "
                               + std::to_string(cols.size()));
    for (size_t j = 0; j < 6; ++j)
      out[i][j] = parse_strict_double(cols[j]);
  }
  return out;
}

static std::string statement_key(const std::string& stmt) {
  const std::string trimmed = trim(stmt);
  if (trimmed.empty() || trimmed.back() != ';')
    throw std::runtime_error("Malformed statement: " + stmt);
  const std::string body = trim(trimmed.substr(0, trimmed.size()-1));
  const auto comma = body.find(',');
  if (comma == std::string::npos)
    throw std::runtime_error("Statement missing ',' in: " + stmt);
  const std::string key = trim(body.substr(0, comma));
  if (key.empty())
    throw std::runtime_error("Statement missing key in: " + stmt);
  return key;
}


static std::string statement_rhs(const std::string& stmt) {
  const std::string trimmed = trim(stmt);
  if (trimmed.empty() || trimmed.back() != ';')
    throw std::runtime_error("Malformed statement: " + stmt);
  const std::string body = trim(trimmed.substr(0, trimmed.size()-1));
  const auto comma = body.find(',');
  if (comma == std::string::npos)
    throw std::runtime_error("Statement missing ',' in: " + stmt);
  return trim(body.substr(comma + 1));
}

static void assign_special_property(Element& cur, const std::string& key,
                                    const std::string& rhs) {
  if (dbg)
    printf("\nassign_special_property: %s\n", cur.name.c_str());
  if (key == "T1") {
    cur.T1 = parse_vec6(rhs, key);
    cur.hasT1 = true;
  } else if (key == "T2") {
    cur.T2 = parse_vec6(rhs, key);
    cur.hasT2 = true;
  } else if (key == "R1") {
    cur.R1 = parse_mat6x6(rhs, key);
    cur.hasR1 = true;
  } else if (key == "R2") {
    cur.R2 = parse_mat6x6(rhs, key);
    cur.hasR2 = true;
  } else {
    throw std::runtime_error("Unknown special property '" + key + "'");
  }
}

static std::pair<std::string, std::vector<Value>>
parse_property_statement(std::string line) {
  line = trim(line);
  if (line.empty()) throw std::runtime_error("Empty property line");
  if (line.back() != ';')
    throw std::runtime_error("Property line missing ';' terminator: " + line);

  line.pop_back();
  line = trim(line);

  auto tokens = split_csv_like(line);
  if (tokens.empty() || tokens[0].empty())
    throw std::runtime_error("Property line missing key: " + line);

  std::string key = tokens[0];
  std::vector<Value> values;
  for (size_t i = 1; i < tokens.size(); ++i) {
    const std::string& tok = tokens[i];
    if (tok.empty()) continue;
    double d = 0.0;
    if (try_parse_double(tok, d)) values.push_back(Value::Num(d));
    else values.push_back(Value::Str(tok));
  }
  return {key, std::move(values)};
}


// Parse 2D table properties where ';' is a row separator within the value data
// (xtable, ytable, xkick, ykick, xkick1, ykick1, B2). Values are flattened.
static std::pair<std::string, std::vector<Value>>
parse_table_statement(std::string line) {
  line = trim(line);
  if (line.empty())
    throw std::runtime_error("Empty table property line");
  // Strip trailing ';' if present.
  if (line.back() == ';') {
    line.pop_back();
    line = trim(line);
  }
  // Key is everything before the first ','.
  auto comma_pos = line.find(',');
  if (comma_pos == std::string::npos)
    throw std::runtime_error("Table property missing ',' separator: " + line);
  std::string key = trim(line.substr(0, comma_pos));
  std::string rest = line.substr(comma_pos + 1);
  // Split by ';' into row segments, then each segment by ',' into values.
  std::vector<Value> values;
  std::string seg;
  auto flush_seg = [&](std::string s) {
    for (const auto& tok : split_csv_like(s)) {
      std::string t = trim(tok);
      if (t.empty()) continue;
      double d = 0.0;
      if (!try_parse_double(t, d))
        throw std::runtime_error("Non-numeric value in " + key + ": " + t);
      values.push_back(Value::Num(d));
    }
  };
  for (char c : rest) {
    if (c == ';') { flush_seg(seg); seg.clear(); }
    else seg.push_back(c);
  }
  if (!trim(seg).empty()) flush_seg(seg);
  return {key, std::move(values)};
}

static void print_value(const Value& v) {
  if (v.isNumber) std::cout << v.number;
  else std::cout << v.text;
}

static void print_elem(const Element &elem) {
  std::cout << "\nelement number " << globval.Cell_nLoc << "\n";
  std::cout << "  ElemName:   " << elem.name << "\n";
  std::cout << "  ElemNbr:    " << elem.number << "\n";
  std::cout << "  PassMethod: " << elem.passMethod << "\n";
  std::cout << "  Properties: " << elem.props.size() << "\n";

  for (const auto& kv : elem.props) {
    std::cout << "    " << kv.first << ": ";
    for (size_t i = 0; i < kv.second.size(); ++i) {
      if (i) std::cout << ", ";
      print_value(kv.second[i]);
    }
    std::cout << "\n";
  }
}


static const std::vector<Value>& require_prop(const Element& e,
                                              const std::string& key) {
  auto it = e.props.find(key);
  if (it == e.props.end())
    throw std::runtime_error
      ("Element '" + e.name + "' missing property '" + key + "'");
  return it->second;
}

static double require_number(const Element& e, const std::string& key,
                             size_t idx = 0) {
  const auto& prop = require_prop(e, key);
  if (idx >= prop.size())
    throw std::runtime_error("Element '" + e.name + "' property '" + key
                             + "' too short");
  const Value& v = prop[idx];
  if (!v.isNumber)
    throw std::runtime_error("Element '" + e.name + "' property '" + key
                             + "' must be numeric");
  return v.number;
}

// Post-processing pass: assign Fnum/Knum/ElemFam from element names.
// Elements sharing the same PName belong to the same family.
static void assign_elem_families()
{
  std::unordered_map<std::string, int> name_to_fnum;
  globval.Elem_nFam = 0;
  bool dbg = false;

  for (long i = 0; i <= globval.Cell_nLoc; i++)
  {
    std::string name(Cell[i].Elem.PName);
    auto result = name_to_fnum.emplace(name, (int)globval.Elem_nFam + 1);
    const bool inserted = result.second;
    const int fnum = result.first->second;

    if (inserted)
    {
      globval.Elem_nFam++;
      ElemFam[fnum - 1].nKid = 0;
      strcpy(ElemFam[fnum - 1].ElemF.PName, Cell[i].Elem.PName);
    }

    Cell[i].Fnum = fnum;
    Cell[i].Knum = 0;
    ElemFam[fnum - 1].nKid++;
    Cell[i].Knum = ElemFam[fnum - 1].nKid;
    ElemFam[fnum - 1].KidList[Cell[i].Knum - 1] = i;
    if (dbg)
    {
      printf("  ElemName = '%s'\n", Cell[i].Elem.PName);
      printf("  Fnum     = %4d\n", Cell[i].Fnum);
      printf("  Knum     = %4d\n", Cell[i].Knum);
      printf("  nKid     = %4d\n", ElemFam[fnum - 1].nKid);
    }
    if (Cell[i].Knum == 1)
      ElemFam[fnum - 1].ElemF = Cell[i].Elem;
  }
  if (dbg)
    printf("\nRead in %d elements with %d families.\n",
           globval.Cell_nLoc, globval.Elem_nFam);
}

static void create_elem(const Element &curr_elem)
{
  CellType &cell = Cell[globval.Cell_nLoc];
  elemtype &elem = cell.Elem;

  cell.Fnum = 0;
  cell.Knum = 0;

  cell.dS[X_] = 0e0;
  cell.dS[Y_] = 0e0;
  cell.dT[X_] = 1e0;
  cell.dT[Y_] = 0e0;

  string_to_c_str(curr_elem.name, elem.PName);

  // Allocate element.
  if ((curr_elem.passMethod == "IdentityPass")
      || (curr_elem.passMethod == "AperturePass")) {
    // Marker.
    elem.PL = 0e0;
    elem.Pkind = PartsKind(marker);
  } else if (curr_elem.passMethod == "DriftPass") {
    // Drift.
    elem.Pkind = PartsKind(drift);
    Drift_Alloc(&elem);
  } else if ((curr_elem.passMethod == "CorrectorPass") ||
	     (curr_elem.passMethod == "StrMPoleSymplectic4Pass") ||
	     (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    elem.Pkind = PartsKind(Mpole);
    // Multipole.
    Mpole_Alloc(&elem);
  } else if (curr_elem.passMethod == "RFCavityPass") {
    // RF Cavity.
    elem.Pkind = PartsKind(Cavity);
    Cav_Alloc(&elem);
  }
  else if (curr_elem.passMethod == "IdTablePass" ||
           curr_elem.passMethod == "IdTableRadPass")
  {
    // Insertion Device (kick map).
    elem.Pkind = PartsKind(Insertion);
    Insertion_Alloc(&elem);
  }
  else if (curr_elem.passMethod == "GWigSymplecticPass")
  {
    // GWigSymplecticPass    - analytic,
    // GWigSymplecticRadPass - analytic.
    throw std::runtime_error("create_elem: unsupported PassMethod '" +
                             curr_elem.passMethod + "'");
  } else {
    throw std::runtime_error("create_elem: unsupported PassMethod '" +
                             curr_elem.passMethod + "'");
  }


  // Set element properties.
  if (dbg) {
    printf("\ncreate_elem: %4ld %2d\n", globval.Cell_nLoc, elem.Pkind);
    printf("  %s\n", elem.PName);
  }
  if (elem.Pkind != marker) {
    elem.PL = require_number(curr_elem, "Length");
    if (dbg) printf("  L          = %9.3e\n", elem.PL);
  }
  if ((curr_elem.passMethod != "IdentityPass") &&
      (curr_elem.passMethod != "CorrectorPass")) {
    auto it = curr_elem.props.find("EApertures");
    if (it != curr_elem.props.end() && !it->second.empty()) {
      if (it->second.size() < 2)
        throw std::runtime_error
	  ("Element '" + curr_elem.name + "' property 'EApertures' too short");
      auto X_max = require_number(curr_elem, "EApertures", 0);
      auto Y_max = require_number(curr_elem, "EApertures", 1);
      if (dbg) printf("  EApertures = [%9.3e, %9.3e]\n", X_max, Y_max);
    }
  }
  if (curr_elem.passMethod != "AperturePass") {
    auto it = curr_elem.props.find("Limits");
    double limits[2][2];
    if (it != curr_elem.props.end() && !it->second.empty()) {
      if (it->second.size() < 4)
        throw std::runtime_error
	  ("Element '" + curr_elem.name + "' property 'Limits' too short");
      limits[0][0] = require_number(curr_elem, "Limits", 0);
      limits[0][1] = require_number(curr_elem, "Limits", 1);
      limits[1][0] = require_number(curr_elem, "Limits", 2);
      limits[1][1] = require_number(curr_elem, "Limits", 3);
      if (dbg)
	printf("  Limits = %9.3e %9.3e %9.3e %9.3e\n",
	       limits[0][0], limits[0][1], limits[1][0], limits[1][1]);
    }
  }
  if (curr_elem.passMethod == "CorrectorPass") {
    if (elem.PL == 0e0)
      elem.M->Pthick = pthicktype(thin);
    else
      elem.M->Pthick = pthicktype(thick);
  }
  if ((curr_elem.passMethod == "StrMPoleSymplectic4Pass") ||
      (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    if (elem.PL == 0e0)
      elem.M->Pthick = pthicktype(thin);
    else
      elem.M->Pthick = pthicktype(thick);
    if ((curr_elem.passMethod == "BndMPoleSymplectic4RadPass")
	&& (elem.M->Pthick == thick)){
          auto phi = require_number(curr_elem, "BendingAngle");
	  elem.M->Pirho = phi/elem.PL;
	  if (dbg)
	    printf("  phi        = %10.3e\n", phi*180e0/M_PI);
          auto phi_1 = require_number(curr_elem, "EntranceAngle");
	  elem.M->PTx1 = phi_1;
	  if (dbg)
	    printf("  phi_1      = %10.3e\n", phi_1*180e0/M_PI);
          auto phi_2 = require_number(curr_elem, "ExitAngle");
	  elem.M->PTx2 = phi_2;
	  if (dbg)
	    printf("  phi_2      = %10.3e\n", phi_2*180e0/M_PI);
    }
    auto n_int = (int)std::round(require_number(curr_elem, "NumIntSteps"));
    auto max_order = (int)std::round(require_number(curr_elem, "MaxOrder"));
    elem.M->PN = n_int;
    elem.M->Porder = max_order + 1;
    if (dbg) {
      printf("  n_int      = %d\n", elem.M->PN);
      printf("  max_order  = %d\n", elem.M->Porder);
    }
    auto it_an = curr_elem.props.find("PolynomA");
    auto it_bn = curr_elem.props.find("PolynomB");
    if (dbg) printf("   n       b_n         a_n\n");
    for (auto n = 1; n <= elem.M->Porder; n++) {
      elem.M->PB[HOMmax+n] = require_number(curr_elem, "PolynomB", n-1);
      elem.M->PB[HOMmax-n] = require_number(curr_elem, "PolynomA", n-1);
      elem.M->PBpar[HOMmax+n] = elem.M->PB[HOMmax+n];
      elem.M->PBpar[HOMmax-n] = elem.M->PB[HOMmax-n];
      if (dbg)
	printf("  %2d   %10.3e  %10.3e]\n",
	       n, elem.M->PB[HOMmax+n], elem.M->PB[HOMmax-n]);
    }
    // Set n_design based on element name prefix.
    // TODO: Idealy this should depend on Porder.
    //It may not be reliable as PolynomB length gets padded with zeroes.
    switch (curr_elem.name[0])
    {
    case 'D':
    case 'R':
      elem.M->n_design = 1;
      break;
    case 'Q':
      elem.M->n_design = 2;
      break;
    case 'S':
      elem.M->n_design = 3;
      break;
    case 'O':
      elem.M->n_design = 4;
      break;
    default:
      elem.M->n_design = 0;
      break;
    }
  }
  if (curr_elem.passMethod == "RFCavityPass") {
    // RF Cavity.
    auto V_RF = require_number(curr_elem, "Voltage");
    auto f_RF = require_number(curr_elem, "Frequency");
    auto E_0 = require_number(curr_elem, "Energy");
    // TODO: TimeLag needs to be translated to RF phase.
    auto TimeLag = require_number(curr_elem, "TimeLag");

    globval.Energy = 1e-9 * E_0;
    elem.C->V_RF = V_RF; // [V]
    elem.C->f_RF = f_RF; // [Hz]

    if (dbg) {
      printf("  V_RF       = %9.3e\n", V_RF);
      printf("  f_RF       = %9.3e\n", f_RF);
      printf("  E_0        = %9.3e\n", E_0);
      printf("  TimeLag    = %9.3e\n", TimeLag);
    }
  }
  if (curr_elem.passMethod == "IdTablePass" ||
      curr_elem.passMethod == "IdTableRadPass")
  {
    InsertionType *ID = elem.ID;

    const auto &xtab = curr_elem.props.find("xtable")->second;
    const auto &ytab = curr_elem.props.find("ytable")->second;
    const int nx = (int)xtab.size();
    const int nz = (int)ytab.size();

    if (nx > IDXMAX || nz > IDZMAX)
    {
      printf("create_elem: ID table too large:"
             " nx=%d (max %d), nz=%d (max %d)\n",
             nx, IDXMAX, nz, IDZMAX);
      exit(1);
    }

    ID->nx = nx;
    ID->nz = nz;
    for (int j = 0; j < nx; j++)
      ID->tabx[j] = xtab[j].number;
    // LinearInterpolation2 expects tabz in decreasing order; AT ytable is
    // increasing, so reverse it.
    for (int i = 0; i < nz; i++)
      ID->tabz[i] = ytab[nz - 1 - i].number;
    if (dbg)
      printf("  tabx: [%9.3e, %9.3e], tabz: [%9.3e, %9.3e]\n",
             ID->tabx[0], ID->tabx[nx-1], ID->tabz[nz-1], ID->tabz[0]);

    // Energy must be set before kick/B2 normalization.
    auto it_e = curr_elem.props.find("Energy");
    if (it_e != curr_elem.props.end())
      globval.Energy = 1e-9 * it_e->second.at(0).number;
    const double Brho = globval.Energy * 1e9 / c0;
    const double Brho2 = Brho * Brho;
    if (dbg) {
      printf("  Energy     = %9.3e [GeV]\n", globval.Energy);
      printf("  Brho       = %9.3e [T.m]\n", Brho);
    }

    // Second order kick maps (always present).
    // AT stores kicks divided by Brho^2; Tracy's Insertion_Pass rescales by
    // 1/Brho^2 at tracking time, so multiply back here.
    const auto &xk = curr_elem.props.find("xkick")->second;
    const auto &yk = curr_elem.props.find("ykick")->second;
    for (int i = 0; i < nz; i++)
      for (int j = 0; j < nx; j++)
      {
        ID->thetax[i][j] = xk[(nz - 1 - i) * nx + j].number * Brho2;
        ID->thetaz[i][j] = yk[(nz - 1 - i) * nx + j].number * Brho2;
      }
    ID->secondorder = true;
    if (dbg) {
      double xk_min = ID->thetax[0][0], xk_max = ID->thetax[0][0];
      double yk_min = ID->thetaz[0][0], yk_max = ID->thetaz[0][0];
      for (int i = 0; i < nz; i++)
        for (int j = 0; j < nx; j++) {
          xk_min = std::min(xk_min, ID->thetax[i][j]);
          xk_max = std::max(xk_max, ID->thetax[i][j]);
          yk_min = std::min(yk_min, ID->thetaz[i][j]);
          yk_max = std::max(yk_max, ID->thetaz[i][j]);
        }
      printf("  xkick2:    [%9.3e, %9.3e] [T.m]\n", xk_min, xk_max);
      printf("  ykick2:    [%9.3e, %9.3e] [T.m]\n", yk_min, yk_max);
    }

    // First order kick maps (optional).
    auto it_xk1 = curr_elem.props.find("xkick1");
    auto it_yk1 = curr_elem.props.find("ykick1");
    if (it_xk1 != curr_elem.props.end() && it_yk1 != curr_elem.props.end() &&
        !it_xk1->second.empty() && !it_yk1->second.empty())
    {
      const auto &xk1 = it_xk1->second;
      const auto &yk1 = it_yk1->second;
      for (int i = 0; i < nz; i++)
        for (int j = 0; j < nx; j++)
        {
          ID->thetax1[i][j] = xk1[(nz - 1 - i) * nx + j].number;
          ID->thetaz1[i][j] = yk1[(nz - 1 - i) * nx + j].number;
        }
      ID->firstorder = true;
      if (dbg) {
        double xk1_min = ID->thetax1[0][0], xk1_max = ID->thetax1[0][0];
        double yk1_min = ID->thetaz1[0][0], yk1_max = ID->thetaz1[0][0];
        for (int i = 0; i < nz; i++)
          for (int j = 0; j < nx; j++) {
            xk1_min = std::min(xk1_min, ID->thetax1[i][j]);
            xk1_max = std::max(xk1_max, ID->thetax1[i][j]);
            yk1_min = std::min(yk1_min, ID->thetaz1[i][j]);
            yk1_max = std::max(yk1_max, ID->thetaz1[i][j]);
          }
        printf("  xkick1:    [%9.3e, %9.3e] [rad]\n", xk1_min, xk1_max);
        printf("  ykick1:    [%9.3e, %9.3e] [rad]\n", yk1_min, yk1_max);
      }
    }
    else
      ID->firstorder = false;

    ID->Pmethod = Meth_First;
    ID->PN = curr_elem.props.find("Nslice")->second.at(0).number;
    ID->linear = true;
    ID->scaling = 1.0;
    if (dbg)
      printf("  Nslice     = %d\n", ID->PN);

    // B2 field map (optional, for radiation via IdTableRadPass).
    auto it_b2 = curr_elem.props.find("B2");
    if (it_b2 != curr_elem.props.end() && !it_b2->second.empty())
    {
      const auto &b2 = it_b2->second;
      for (int i = 0; i < nz; i++)
        for (int j = 0; j < nx; j++)
          ID->B2[i][j] = b2[(nz - 1 - i) * nx + j].number;
      ID->long_comp = true;
    }
    else
      ID->long_comp = false;

    if (dbg)
      printf("  ID: nx=%d, nz=%d, 1st=%d, 2nd=%d, B2=%d\n",
             nx, nz, ID->firstorder, ID->secondorder, ID->long_comp);
  }
  if ((curr_elem.passMethod == "DriftPass") || (curr_elem.passMethod == "CorrectorPass") || (curr_elem.passMethod == "StrMPoleSymplectic4Pass") || (curr_elem.passMethod == "BndMPoleSymplectic4RadPass"))
  {
    // Misalignment / entrance-exit transforms are parsed and validated;
    // wire them into Tracy-specific fields here if needed by the local API.
    (void)curr_elem.hasT1;
    (void)curr_elem.hasT2;
    (void)curr_elem.hasR1;
    (void)curr_elem.hasR2;
  }

  if (globval.Cell_nLoc == 0)
      cell.S = 0e0;
    else
      cell.S = Cell[globval.Cell_nLoc-1].S + elem.PL;
}

void rdmfile_at(const std::string& file_name) {
  std::ifstream in(file_name);
  if (!in) throw std::runtime_error("Failed to open file: " + file_name);

  Element cur{};
  bool hasCur = false;

  auto fail = [&](const std::string& msg, std::size_t lineNo) -> void {
    throw std::runtime_error(msg + " at line " + std::to_string(lineNo));
  };

  auto require_cur = [&](const char* what, std::size_t lineNo) -> void {
    if (!hasCur) fail(std::string(what) + " before ElemName", lineNo);
  };

  auto start_element = [&](std::string name) -> void {
    hasCur = true;
    cur = Element{};
    cur.name = std::move(name);
  };

  auto flush_cur = [&]() {
    if (!hasCur) return;

    if (cur.name.empty())
      throw std::runtime_error("Element missing ElemName");
    if (cur.number < 0)
      throw std::runtime_error
        ("Element '" + cur.name + "' missing/invalid ElemNbr");
    if (cur.passMethod.empty())
      throw std::runtime_error
        ("Element '" + cur.name + "' missing PassMethod");

    globval.Cell_nLoc++;
    if (false && dbg) print_elem(cur);
    create_elem(cur);

    cur = Element{};
    hasCur = false;
  };

  auto handle_statement =
    [&](const std::string& stmt, std::size_t lineNo) -> void {
    const std::string key = statement_key(stmt);

    if (key == "T1" || key == "T2" || key == "R1" || key == "R2") {
      require_cur(key.c_str(), lineNo);
      assign_special_property(cur, key, statement_rhs(stmt));
      return;
    }

    const bool is_table_key =
      (key == "xtable" || key == "ytable" ||
       key == "xkick"  || key == "ykick"  ||
       key == "xkick1" || key == "ykick1" || key == "B2");
    auto kv = is_table_key ? parse_table_statement(stmt)
                           : parse_property_statement(stmt);
    const auto& vals = kv.second;

    if (key == "ElemName") {
      if (vals.size() != 1 || vals[0].isNumber)
        fail("ElemName expects one text value", lineNo);
      flush_cur();
      start_element(vals[0].text);
      return;
    }

    if (key == "ElemNbr") {
      require_cur("ElemNbr", lineNo);
      if (vals.size() != 1 || !vals[0].isNumber)
        fail("ElemNbr expects one numeric value", lineNo);
      cur.number = (int)std::lround(vals[0].number);
      return;
    }

    if (key == "PassMethod") {
      require_cur("PassMethod", lineNo);
      if (vals.size() != 1 || vals[0].isNumber)
        fail("PassMethod expects one text value", lineNo);
      cur.passMethod = vals[0].text;
      return;
    }

    require_cur("Property", lineNo);

    auto& dst = cur.props[key];
    dst.insert(dst.end(), vals.begin(), vals.end());
  };

  std::string line;
  std::size_t lineNo = 0;
  std::string stmt;
  std::size_t stmtStartLine = 0;

  globval.Cell_nLoc = -1;

  while (std::getline(in, line)) {
    ++lineNo;

    const std::string t = trim(line);
    if (t.empty()) continue;

    if (stmt.empty())
      stmtStartLine = lineNo;
    else
      stmt += ' ';
    stmt += t;

    // Some AT flat files may emit empty optional kick-map lines without a
    // trailing ';' (e.g. "xkick1," / "ykick1,"). Treat these as complete
    // one-line statements so they don't swallow the next property key.
    if (stmt.find(';') == std::string::npos &&
        (starts_with(stmt, "xkick1,") || starts_with(stmt, "ykick1,"))) {
      handle_statement(stmt + ";", stmtStartLine);
      stmt.clear();
      stmtStartLine = 0;
      continue;
    }

    if (stmt.find(';') == std::string::npos)
      continue;

    handle_statement(stmt, stmtStartLine);
    stmt.clear();
    stmtStartLine = 0;
  }

  if (!stmt.empty())
    fail("Unterminated statement", stmtStartLine);

  flush_cur();


  globval.dPcommon = 1e-8;
  globval.CODeps = 1e-14;
  globval.CODimax = 40;

  SI_init();

  globval.mat_meth = false;

  printf("\nrdmfile_at: read %ld elements, C = %7.5f\n",
         globval.Cell_nLoc+1, Cell[globval.Cell_nLoc].S);

  assign_elem_families();

  // Compute harmonic number for all cavity elements now that C is known.
  {
    const double C_ring = Cell[globval.Cell_nLoc].S;
    for (long i = 0; i <= globval.Cell_nLoc; ++i)
    {
      if (Cell[i].Elem.Pkind == PartsKind(Cavity))
      {
        Cell[i].Elem.C->harm_num =
            (int)std::round(Cell[i].Elem.C->f_RF * C_ring / c0);
        if (dbg)
          printf("rdmfile_at: cavity element %ld harm_num=%d\n",
                 i, Cell[i].Elem.C->harm_num);
      }
    }
  }
}
