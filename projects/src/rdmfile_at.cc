
// Requires C++11 but no later.

#if 1
#include <cerrno>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "tracy_lib.h"
#endif

static const bool dbg = false;

void string_to_c_str(const std::string &str, partsName &c_str) {
  // Tracy-2 element names are not "\0" terminated C strings (Pascal legacy).
  if (str.size() > NameLength)
    throw std::runtime_error("ElemName too long for fixed buffer");
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
  std::unordered_map<std::string, std::vector<Value>> props;
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

  out = val;
  return true;
}

static std::pair<std::string, std::vector<Value>>
parse_property_line(std::string line) {
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

static void print_elem(Element &elem) {
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

static void create_elem(Element &curr_elem)
{
  CellType *cell = &Cell[globval.Cell_nLoc];
  elemtype *elem = &cell->Elem;

  // While PName is fixed size array - i.e., Pascal legacy - to keep it tidy. 
  elem->PName[0] = '\0';

  cell->Fnum = 0;
  cell->Knum = 0;

  cell->dS[X_] = 0e0;
  cell->dS[Y_] = 0e0;
  cell->dT[X_] = 1e0;
  cell->dT[Y_] = 0e0;

  string_to_c_str(curr_elem.name, elem->PName);

  if ((curr_elem.passMethod == "IdentityPass")
      || (curr_elem.passMethod == "AperturePass")) {
    elem->PL = 0e0;
    elem->Pkind = PartsKind(marker);
  } else if (curr_elem.passMethod == "DriftPass") {
    elem->Pkind = PartsKind(drift);
    Drift_Alloc(elem);
  } else if ((curr_elem.passMethod == "CorrectorPass")
	     || (curr_elem.passMethod == "StrMPoleSymplectic4Pass")
	     || (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    elem->Pkind = PartsKind(Mpole);
    Mpole_Alloc(elem);
 } else if (curr_elem.passMethod == "RFCavityPass") {
    elem->Pkind = PartsKind(Cavity);
    Cav_Alloc(elem);
  } else {
    std::cout << "create_elem: *** undef. passMethod - "
	      << curr_elem.passMethod << "\n";
    exit(1);
  }

  if (dbg) {
    printf("\ncreate_elem: %4ld %2d\n", globval.Cell_nLoc, elem->Pkind);
    printf("  %s\n", elem->PName);
  }
  if (elem->Pkind != marker) {
    auto L = curr_elem.props.find("Length")->second.at(0).number;
    elem->PL = L;
    if (dbg) printf("  L          = %9.3e\n", elem->PL);
  }
  if ((curr_elem.passMethod != "IdentityPass")
      && (curr_elem.passMethod != "CorrectorPass")) {
    auto it = curr_elem.props.find("EApertures");
    if (it != curr_elem.props.end() && !it->second.empty()) {
      auto X_max = it->second.at(0).number;
      auto Y_max = it->second.at(1).number;
      if (dbg) printf("  EApertures = [%9.3e, %9.3e]\n", X_max, Y_max);
    }
  }
  if (curr_elem.passMethod != "AperturePass") {
    auto it = curr_elem.props.find("Limits");
    double limits[2][2];
    if (it != curr_elem.props.end() && !it->second.empty()) {
      limits[0][0] = it->second.at(0).number;
      limits[0][1] = it->second.at(1).number;
      limits[1][0] = it->second.at(2).number;
      limits[1][1] = it->second.at(3).number;
      if (dbg)
	printf("  Limits = %9.3e %9.3e %9.3e %9.3e\n",
	       limits[0][0], limits[0][1], limits[1][0], limits[1][1]);
    }
  }
  if ((curr_elem.passMethod == "StrMPoleSymplectic4Pass")
      || (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    if (elem->PL == 0e0)
      elem->M->Pthick = pthicktype(thin);
    else
      elem->M->Pthick = pthicktype(thick);
    auto n_int =
      (int)std::round(curr_elem.props.find("NumIntSteps")->second.at(0).number);
    auto max_order =
      (int)std::round(curr_elem.props.find("MaxOrder")->second.at(0).number);
    elem->M->PN = n_int;
    elem->M->Porder = max_order + 1;
    if (dbg) {
      printf("  n_int      = %d\n", elem->M->PN);
      printf("  max_order  = %d\n", elem->M->Porder);
    }
    auto it_an = curr_elem.props.find("PolynomA");
    auto it_bn = curr_elem.props.find("PolynomB");
    if (dbg) printf("   n       b_n         a_n\n");
    for (auto n = 1; n <= elem->M->Porder; n++) {
      elem->M->PB[HOMmax+n] = it_bn->second.at(n-1).number;
      elem->M->PB[HOMmax-n] = it_an->second.at(n-1).number;
      elem->M->PBpar[HOMmax+n] = elem->M->PB[HOMmax+n];
      elem->M->PBpar[HOMmax-n] = elem->M->PB[HOMmax-n];
      if (dbg)
	printf("  %2d   %10.3e  %10.3e]\n",
	       n, elem->M->PB[HOMmax+n], elem->M->PB[HOMmax-n]);
    }
  }
  if (curr_elem.passMethod == "RFCavityPass") {
    auto V_RF = curr_elem.props.find("Voltage")->second.at(0).number;
    auto f_RF = curr_elem.props.find("Frequency")->second.at(0).number;
    auto E_0 = curr_elem.props.find("Energy")->second.at(0).number;
    
    globval.Energy = 1e-9*E_0;

    if (dbg) {
      printf("  V_RF       = %9.3e\n", V_RF);
      printf("  f_RF       = %9.3e\n", f_RF);
      printf("  E_0        = %9.3e\n", E_0);
    }
  }

  if (globval.Cell_nLoc == 0)
      cell->S = 0e0;
    else
      cell->S = cell->S + elem->PL;
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
    // If Element::number defaults to 0 in your struct, you may want to
    // explicitly mark it unset instead:
    // cur.number = -1;
  };

  // Finalize + emit + clear the current element.
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

  // Helper: parse lines like "Key: value".
  auto parse_colon_field = [&](const std::string& t, const std::string& key)
    -> std::string {
    return trim(t.substr(key.size()));
  };

  std::string line;
  std::size_t lineNo = 0;

  globval.Cell_nLoc = -1;

  while (std::getline(in, line)) {
    ++lineNo;

    const std::string t = trim(line);
    if (t.empty()) continue;

    if (starts_with(t, "ElemName:")) {
      flush_cur();
      start_element(parse_colon_field(t, "ElemName:"));
      continue;
    }

    if (starts_with(t, "ElemNbr:")) {
      require_cur("ElemNbr", lineNo);
      const std::string rhs = parse_colon_field(t, "ElemNbr:");
      try {
        cur.number = std::stoi(rhs);
      } catch (const std::exception&) {
        fail("Invalid ElemNbr ('" + rhs + "')", lineNo);
      }
      continue;
    }

    if (starts_with(t, "PassMethod:")) {
      require_cur("PassMethod", lineNo);
      cur.passMethod = parse_colon_field(t, "PassMethod:");
      continue;
    }

    // Otherwise it must be a property line.
    require_cur("Property", lineNo);

    auto kv = parse_property_line(t);
    auto& dst = cur.props[kv.first];
    dst.insert(dst.end(), kv.second.begin(), kv.second.end());
  }

  // Flush the last element at EOF.
  flush_cur();

  globval.dPcommon = 1e-8;
  globval.CODeps = 1e-14;
  globval.CODimax = 40;

  SI_init();

  globval.mat_meth = !false;
}
