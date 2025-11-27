
const bool dbg = !false;


std::string get_ElemName(const std::string &line)
{
  std::string        key, val;
  std::istringstream iss;

  iss.clear();
  iss.str(line);
  getline(iss, key, ':');
  if (key == "ElemName") {
    getline(iss, val);
    if (dbg) cout << setw(SymbolLength+1) << val;
  } else {
    std::cerr << "*** No ElemName parameter!\n";
    exit(1);
  }
  return val;
}


int get_ElemNbr(const std::string &line, const std::string &name)
{
  std::string        key, val;
  std::istringstream iss;
  int                loc;

  iss.clear();
  iss.str(line);
  getline(iss, key, ':');
  if (key == "ElemNbr") {
    getline(iss, val);
    loc = atoi(val.c_str());
    if ((loc < 1 || loc > Cell_nLocMax) || (loc < 1 || loc > Elem_nFamMax)) {
      std::cerr
	<< "***: Cell_nLocMax or Elem_nFamMax exceeded (" << loc << ")!\n";
      exit(1);
    }
    globval.Cell_nLoc = loc;

    strncpy(Cell[loc].Elem.PName, name.c_str(), NameLength);
    Cell[loc].Fnum = loc;
    Cell[loc].Knum = 1;

    auto Fnum = loc;
    strncpy(ElemFam[Fnum-1].ElemF.PName, name.c_str(), NameLength);
    ElemFam[Fnum-1].nKid = 1;
    ElemFam[Fnum-1].KidList[0] = loc;
    globval.Elem_nFam = loc;


    if (dbg) cout << setw(5) << globval.Cell_nLoc;
  } else {
    std::cerr << "*** No ElemNbr parameter!\n";
    exit(1);
  }

  return loc;
}


void get_PassMethod(const std::string &line, const int loc)
{
  std::string        key, val;
  std::istringstream iss;

  iss.clear();
  iss.str(line);
  getline(iss, key, ':');
  if (key == "PassMethod") {
    getline(iss, val);
    if (val.find("DriftPass") != std::string::npos) {
      Cell[loc].Elem.Pkind = drift;
      Drift_Alloc(&Cell[loc].Elem);
    } else if (val.find("IdentityPass") != std::string::npos ||
	       val.find("AperturePass") != std::string::npos) {
      Cell[loc].Elem.Pkind = marker;
    } else if (val.find("StrMPoleSymplectic4Pass") != std::string::npos ||
	       val.find("BndMPoleSymplectic4RadPass") != std::string::npos) {
      Cell[loc].Elem.Pkind = Mpole;
      Mpole_Alloc(&Cell[loc].Elem);
      Cell[loc].Elem.M->Pthick = thick;
    } else if (val.find("CorrectorPass") != std::string::npos) {
      Cell[loc].Elem.Pkind = Mpole;
      Mpole_Alloc(&Cell[loc].Elem);
      Cell[loc].Elem.M->Pthick = thin;
      Cell[loc].Elem.M->Porder = 1;
    } else if (val.find("RFCavityPass") != std::string::npos) {
      Cell[loc].Elem.Pkind = Cavity;
      Cav_Alloc(&Cell[loc].Elem);
    } else {
      Cell[loc].Elem.Pkind = undef;
      std::cerr << "*** undefined PassMethod: " << val << "!\n";
      exit(1);
    }

    if (dbg) cout << setw(2) << Cell[loc].Elem.Pkind;
  }
}


void get_Length(const std::string line, const int loc)
{
  std::string        key, val;
  std::istringstream iss;

  iss.clear();
  iss.str(line);
  getline(iss, key, ',');
  if (Cell[loc].Elem.Pkind != marker) {
    if (key.find("Length") != std::string::npos) {
      std::string val;
      getline(iss, val);
      Cell[loc].Elem.PL = atof(val.c_str());
      if (dbg)
	cout << scientific << setprecision(3)
	     << setw(10) << Cell[loc].Elem.PL << "\n";
    } else {
      std::cerr << "*** No Length parameter!\n";
      exit(1);
    }
  }
}


void rdmfile_new(const char *filename) {
  std::string   line, val;
  std::ifstream inf(filename);

  if (!inf.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    exit(1);
  }

  std::cout << "Reading flat file: " << filename << std::endl;

  while (std::getline(inf, line)) {
    if (line.empty()) continue;
    val = get_ElemName(line);

    std::getline(inf, line);
    auto loc = get_ElemNbr(line, val);

    std::getline(inf, line);
    get_PassMethod(line, loc);
 
    std::getline(inf, line);
    get_Length(line, loc);

    if (dbg) cout << "\n";
  }

    // else if (key.find("EApertures") != std::string::npos) {
    //   std::string x1, x2, y1, y2;
    //   getline(iss, x1, ',');
    //   getline(iss, x2, ',');
    //   getline(iss, y1, ',');
    //   getline(iss, y2);
    //   Cell[loc].maxampl[X_][0] = atof(x1.c_str());
    //   Cell[loc].maxampl[X_][1] = atof(x2.c_str());
    //   Cell[loc].maxampl[Y_][0] = atof(y1.c_str());
    //   Cell[loc].maxampl[Y_][1] = atof(y2.c_str());
    // }

    // else if (key.find("PolynomB") != std::string::npos) {
    //   int n = 0;
    //   std::string val;
    //   while (getline(iss, val, ',')) {
    //     double b = atof(val.c_str());
    //     if (Cell[loc].Elem.Pkind == Mpole && b != 0.0) {
    //       Cell[loc].Elem.M->PB[HOMmax + n] = b;
    //       Cell[loc].Elem.M->PBpar[HOMmax + n] = b;
    //       Cell[loc].Elem.M->Porder = std::max(n, Cell[loc].Elem.M->Porder);
    //     }
    //     n++;
    //   }
    // }

    // else if (key.find("NumIntSteps") != std::string::npos) {
    //   std::string val;
    //   getline(iss, val);
    //   if (Cell[loc].Elem.Pkind == Mpole)
    //     Cell[loc].Elem.M->PN = atoi(val.c_str());
    // }

    // else if (key.find("Voltage") != std::string::npos) {
    //   std::string val;
    //   getline(iss, val);
    //   if (Cell[loc].Elem.Pkind == Cavity)
    //     Cell[loc].Elem.C->V_RF = atof(val.c_str());
    // }

    // else if (key.find("Frequency") != std::string::npos) {
    //   std::string val;
    //   getline(iss, val);
    //   if (Cell[loc].Elem.Pkind == Cavity)
    //     Cell[loc].Elem.C->f_RF = atof(val.c_str());
    // }

    // else if (key.find("HarmonicNumber") != std::string::npos) {
    //   std::string val;
    //   getline(iss, val);
    //   if (Cell[loc].Elem.Pkind == Cavity)
    //     Cell[loc].Elem.C->harm_num = atoi(val.c_str());
    // }

    // else if (key.find("PhaseLag") != std::string::npos) {
    //   std::string val;
    //   getline(iss, val);
    //   if (Cell[loc].Elem.Pkind == Cavity)
    //     Cell[loc].Elem.C->phi_RF = atof(val.c_str());
    // }
 
  std::cout << "rdmfile_new: read " << globval.Cell_nLoc << " elements."
	    << "\n";
  inf.close();
}
