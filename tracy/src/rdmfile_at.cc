
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
    IdTablePass
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
  if (str.size() > NameLength)
    throw std::runtime_error("ElemName too long for fixed buffer");
  std::memset(c_str, 0, sizeof(c_str));
  std::memcpy(c_str, str.data(), str.size());
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
  } else if (curr_elem.passMethod == "GWigSymplecticPass") {
    // GWigSymplecticPass    - analytic,
    // GWigSymplecticRadPass - analytic,
    // IdTablePass           - kick map.
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
  }
  if (curr_elem.passMethod == "RFCavityPass") {
    // RF Cavity.
    auto V_RF = require_number(curr_elem, "Voltage");
    auto f_RF = require_number(curr_elem, "Frequency");
    auto E_0 = require_number(curr_elem, "Energy");
    
    globval.Energy = 1e-9*E_0;

    if (dbg) {
      printf("  V_RF       = %9.3e\n", V_RF);
      printf("  f_RF       = %9.3e\n", f_RF);
      printf("  E_0        = %9.3e\n", E_0);
    }
  }
  if ((curr_elem.passMethod == "DriftPass")
      || (curr_elem.passMethod == "CorrectorPass")
      || (curr_elem.passMethod == "StrMPoleSymplectic4Pass")
      || (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
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

    auto kv = parse_property_statement(stmt);
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
}

