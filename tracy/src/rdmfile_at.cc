
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

// Parse property lines that represent tables where ';' is a value/row separator.
// Handles: xtable/ytable (1D) and xkick/ykick/xkick1/ykick1/B2 (2D, rows by ';', cols by ',').
// Values are flattened into a single vector.
static std::pair<std::string, std::vector<Value>>
parse_table_line(std::string line) {
  line = trim(line);
  if (line.empty())
    throw std::runtime_error("Empty property line");

  // Strip trailing ';' if present.
  if (line.back() == ';') {
    line.pop_back();
    line = trim(line);
  }

  // Key is everything before the first ','.
  auto comma_pos = line.find(',');
  if (comma_pos == std::string::npos)
    throw std::runtime_error("Property line missing key separator: " + line);

  std::string key = trim(line.substr(0, comma_pos));
  std::string rest = line.substr(comma_pos + 1);

  // Split by ';' into segments, then each segment by ',' into values.
  std::vector<Value> values;
  std::vector<std::string> segments;
  std::string cur;
  for (char c : rest) {
    if (c == ';') {
      segments.push_back(cur);
      cur.clear();
    } else
      cur.push_back(c);
  }
  if (!trim(cur).empty())
    segments.push_back(cur);

  for (const auto &seg : segments) {
    auto tokens = split_csv_like(seg);
    for (const auto &tok : tokens) {
      std::string t = trim(tok);
      if (t.empty())
        continue;
      double d = 0.0;
      if (!try_parse_double(t, d))
        throw std::runtime_error("Non-numeric value in " + key + ": " + t);
      values.push_back(Value::Num(d));
    }
  }

  return {key, std::move(values)};
}

static void print_value(const Value &v) {
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

// Post-processing pass: assign Fnum/Knum/ElemFam from element names, matching
// Elements sharing the same PName belong to the same family.
static void assign_elem_families()
{
  std::unordered_map<std::string, int> name_to_fnum;
  globval.Elem_nFam = 0;

  for (long i = 0; i <= globval.Cell_nLoc; i++) {
    CellType &cell = Cell[i];
    std::string name(cell.Elem.PName);

    auto result = name_to_fnum.emplace(name, (int)globval.Elem_nFam + 1);
    const bool inserted = result.second;
    const int  fnum     = result.first->second;

    if (inserted) {
      // First kid of this family: initialise the family prototype.
      globval.Elem_nFam++;
      memset(ElemFam[fnum-1].ElemF.PName, 0, sizeof(partsName));
      memcpy(ElemFam[fnum-1].ElemF.PName, cell.Elem.PName,
             sizeof(partsName));
      ElemFam[fnum-1].nKid = 0;
    }

    cell.Fnum = fnum;
    ElemFam[fnum-1].nKid++;
    cell.Knum = ElemFam[fnum-1].nKid;

    // KidList and ElemF.Pkind are only populated for i > 0, exclude ring-origin marker.
    if (i > 0) {
      ElemFam[fnum-1].KidList[cell.Knum - 1] = i;
      ElemFam[fnum-1].ElemF.Pkind = cell.Elem.Pkind;
    }
  }
}

static void create_elem(Element &curr_elem)
{
  CellType &cell = Cell[globval.Cell_nLoc];
  elemtype &elem = cell.Elem;

  // While PName is fixed size array - i.e., Pascal legacy - to keep it tidy. 
  elem.PName[0] = '\0';

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
  } else if (curr_elem.passMethod == "IdTablePass" ||
             curr_elem.passMethod == "IdTableRadPass") {
    // Insertion Device (kick map).
    elem.Pkind = PartsKind(Insertion);
    Insertion_Alloc(&elem);
  } else if (curr_elem.passMethod == "GWigSymplecticPass") {
    // GWigSymplecticPass    - analytic,
    // GWigSymplecticRadPass - analytic.
    std::cout << "create_elem: *** undef. passMethod not implemented - "
              << curr_elem.passMethod << "\n";
    exit(1);
  } else {
    std::cout << "create_elem: *** undef. passMethod - "
	      << curr_elem.passMethod << "\n";
    exit(1);
  }


  // Set element properties.
  if (dbg) {
    printf("\ncreate_elem: %4ld %2d\n", globval.Cell_nLoc, elem.Pkind);
    printf("  %s\n", elem.PName);
  }
  if (elem.Pkind != marker) {
    auto L = curr_elem.props.find("Length")->second.at(0).number;
    elem.PL = L;
    if (dbg) printf("  L          = %9.3e\n", elem.PL);
  }
  if ((curr_elem.passMethod != "IdentityPass") &&
      (curr_elem.passMethod != "CorrectorPass")) {
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
  if ((curr_elem.passMethod == "StrMPoleSymplectic4Pass") ||
      (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    if (elem.PL == 0e0)
      elem.M->Pthick = pthicktype(thin);
    else
      elem.M->Pthick = pthicktype(thick);
    if ((curr_elem.passMethod == "BndMPoleSymplectic4RadPass")
	&& (elem.M->Pthick == thick)){
      auto phi = curr_elem.props.find("BendingAngle")->second.at(0).number;
	  elem.M->Pirho = phi/elem.PL;
      if (dbg)
	    printf("  phi        = %10.3e\n", phi*180e0/M_PI);
          auto phi_1 =
	    curr_elem.props.find("EntranceAngle")->second.at(0).number;
      elem.M->PTx1 = phi_1;
      if (dbg)
	    printf("  phi_1      = %10.3e\n", phi_1*180e0/M_PI);
      auto phi_2 = curr_elem.props.find("ExitAngle")->second.at(0).number;
      elem.M->PTx2 = phi_2;
      if (dbg)
	    printf("  phi_2      = %10.3e\n", phi_2*180e0/M_PI);
    }
    auto n_int =
      (int)std::round(curr_elem.props.find("NumIntSteps")->second.at(0).number);
    auto max_order =
        (int)std::round(curr_elem.props.find("MaxOrder")->second.at(0).number);
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
      elem.M->PB[HOMmax+n] = it_bn->second.at(n-1).number;
      elem.M->PB[HOMmax-n] = it_an->second.at(n-1).number;
      elem.M->PBpar[HOMmax+n] = elem.M->PB[HOMmax+n];
      elem.M->PBpar[HOMmax-n] = elem.M->PB[HOMmax-n];
      if (dbg)
	printf("  %2d   %10.3e  %10.3e]\n",
	       n, elem.M->PB[HOMmax+n], elem.M->PB[HOMmax-n]);
    }
  }
  if (curr_elem.passMethod == "RFCavityPass") {
    // RF Cavity.
    auto V_RF = curr_elem.props.find("Voltage")->second.at(0).number;
    auto f_RF = curr_elem.props.find("Frequency")->second.at(0).number;
    auto E_0 = curr_elem.props.find("Energy")->second.at(0).number;

    globval.Energy = 1e-9*E_0;
    elem.C->V_RF   = V_RF;   // [V]
    elem.C->f_RF   = f_RF;   // [Hz]

    if (dbg) {
      printf("  V_RF       = %9.3e\n", V_RF);
      printf("  f_RF       = %9.3e\n", f_RF);
      printf("  E_0        = %9.3e\n", E_0);
    }
  }

  if (curr_elem.passMethod == "IdTablePass" ||
      curr_elem.passMethod == "IdTableRadPass") {
    InsertionType *ID = elem.ID;

    const auto &xtab = curr_elem.props.find("xtable")->second;
    const auto &ytab = curr_elem.props.find("ytable")->second;
    const int nx = (int)xtab.size();
    const int nz = (int)ytab.size();

    if (nx > IDXMAX || nz > IDZMAX) {
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

    // Energy must be set before kick/B2 normalization.
    auto it_e = curr_elem.props.find("Energy");
    if (it_e != curr_elem.props.end())
      globval.Energy = 1e-9 * it_e->second.at(0).number;
    const double Brho = globval.Energy * 1e9 / c0;

    // Second order kick maps (always present).
    // The AT flat file stores kicks already divided by Brho^2, but Tracy's
    // Insertion_Pass applies its own 1/Brho^2 scaling at tracking time.
    // Multiply by Brho^2 here to recover the raw (T^2 m^2) values.
    const double Brho2 = Brho * Brho;
    const auto &xk = curr_elem.props.find("xkick")->second;
    const auto &yk = curr_elem.props.find("ykick")->second;
    for (int i = 0; i < nz; i++)
      for (int j = 0; j < nx; j++) {
        ID->thetax[i][j] = xk[(nz - 1 - i) * nx + j].number * Brho2;
        ID->thetaz[i][j] = yk[(nz - 1 - i) * nx + j].number * Brho2;
      }
    ID->secondorder = true;

    // First order kick maps (optional).
    auto it_xk1 = curr_elem.props.find("xkick1");
    auto it_yk1 = curr_elem.props.find("ykick1");
    if (it_xk1 != curr_elem.props.end() && it_yk1 != curr_elem.props.end() &&
        !it_xk1->second.empty() && !it_yk1->second.empty()) {
      const auto &xk1 = it_xk1->second;
      const auto &yk1 = it_yk1->second;
      for (int i = 0; i < nz; i++)
        for (int j = 0; j < nx; j++) {
          ID->thetax1[i][j] = xk1[(nz - 1 - i) * nx + j].number;
          ID->thetaz1[i][j] = yk1[(nz - 1 - i) * nx + j].number;
        }
      ID->firstorder = true;
    } else
      ID->firstorder = false;

    ID->Pmethod = Meth_First;
    ID->PN = curr_elem.props.find("Nslice")->second.at(0).number;
    ID->linear = true;
    ID->scaling = 1.0;

    // B2 field map (optional, for radiation via IdTableRadPass).
    // Expected in the same units as Tracy's radiate_ID, i.e. already
    // normalized by L*Brho^2.  Store as-is.
    auto it_b2 = curr_elem.props.find("B2");
    if (it_b2 != curr_elem.props.end() && !it_b2->second.empty()) {
      const auto &b2 = it_b2->second;
      for (int i = 0; i < nz; i++)
        for (int j = 0; j < nx; j++) {
          ID->B2[i][j] = b2[(nz - 1 - i) * nx + j].number;
        }
      ID->long_comp = true;
    } else
      ID->long_comp = false;

    if (dbg)
      printf("  ID: nx=%d, nz=%d, 1st=%d, 2nd=%d\n", nx, nz, ID->firstorder,
             ID->secondorder);
  }

  if ((curr_elem.passMethod == "DriftPass")
  || (curr_elem.passMethod == "CorrectorPass")
  || (curr_elem.passMethod == "StrMPoleSymplectic4Pass")
  || (curr_elem.passMethod == "BndMPoleSymplectic4RadPass")) {
    // Misalignment.
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

    // Catch parameter lines comming from 2d properites in AT
    std::pair<std::string, std::vector<Value>> kv;
    if (starts_with(t, "xtable,") || starts_with(t, "ytable,") ||
        starts_with(t, "xkick,") || starts_with(t, "ykick,") ||
        starts_with(t, "xkick1,") || starts_with(t, "ykick1,") ||
        starts_with(t, "B2,"))
      kv = parse_table_line(t);
    else
      kv = parse_property_line(t);
    auto& dst = cur.props[kv.first];
    dst.insert(dst.end(), kv.second.begin(), kv.second.end());
  }

  // Flush the last element at EOF.
  flush_cur();

  in.close();

  globval.dPcommon = 1e-8;
  globval.CODeps = 1e-14;
  globval.CODimax = 40;

  SI_init();

  globval.mat_meth = false;

  printf("\nrdmfile_at: read %ld elements, C = %7.5f\n",
	 globval.Cell_nLoc+1, Cell[globval.Cell_nLoc].S);

  assign_elem_families();

  // Compute harmonic number for all cavity elements now that C is known.
  const double C_ring = Cell[globval.Cell_nLoc].S;
  for (long i = 0; i <= globval.Cell_nLoc; i++) {
    if (Cell[i].Elem.Pkind == PartsKind(Cavity) && Cell[i].Elem.C->f_RF != 0.0)
      Cell[i].Elem.C->harm_num =
        (int)std::round(Cell[i].Elem.C->f_RF * C_ring / c0);
  }
}
