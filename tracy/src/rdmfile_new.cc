#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>

extern "C" {
  #include "tracy_lib.h"
}

void rdmfile_new(const char *filename) {
  std::ifstream inf(filename);
  if (!inf.is_open()) {
    std::cerr << "Cannot open file: " << filename << std::endl;
    exit(1);
  }

  std::string line;
  char        name[NameLength] = {0};
  int         fnum = 0, knum = 1;
  int          i = 0;

  std::cout << "reading machine file (new format): " << filename << std::endl;

  while (std::getline(inf, line)) {
    if (line.empty()) continue;
    cout << line << "\n";

    std::istringstream iss(line);
    std::string key;
    getline(iss, key, ':');

    if (key == "ElemName") {
      std::string val;
      getline(iss, val);
      if (sscanf(val.c_str(), "%29s", name) != 1) {
        std::cerr << "Error parsing ElemName." << std::endl;
        exit(1);
      }
      i++;
      fnum = i;

      if (i < 1 || i > Cell_nLocMax) {
        std::cerr
	  << "Fatal: Cell index out-of-bounds (i = " << i << "). Aborting.\n";
        exit(1);
      }

      strncpy(Cell[i].Elem.PName, name, NameLength);
      Cell[i].Fnum = fnum;
      Cell[i].Knum = knum;

      cout << "\n" << name << "\n";
      strncpy(ElemFam[fnum-1].ElemF.PName, name, NameLength);
      cout << "\nSo far, so good!\n";
      ElemFam[fnum-1].nKid = 1;
      ElemFam[fnum-1].KidList[0] = i;
      globval.Elem_nFam = std::max((long)fnum, globval.Elem_nFam);
    }

    else if (key == "PassMethod") {
      std::string val;
      getline(iss, val);
      if (val.find("DriftPass") != std::string::npos) {
        Cell[i].Elem.Pkind = drift;
        Drift_Alloc(&Cell[i].Elem);
      } else if (val.find("IdentityPass") != std::string::npos ||
                 val.find("AperturePass") != std::string::npos) {
        Cell[i].Elem.Pkind = marker;
      } else if (val.find("StrMPoleSymplectic4Pass") != std::string::npos ||
                 val.find("BndMPoleSymplectic4RadPass") != std::string::npos) {
        Cell[i].Elem.Pkind = Mpole;
        Mpole_Alloc(&Cell[i].Elem);
        Cell[i].Elem.M->Pthick = thick;
      } else if (val.find("CorrectorPass") != std::string::npos) {
        Cell[i].Elem.Pkind = Mpole;
        Mpole_Alloc(&Cell[i].Elem);
        Cell[i].Elem.M->Pthick = thin;
        Cell[i].Elem.M->Porder = 1;
      } else if (val.find("RFCavityPass") != std::string::npos) {
        Cell[i].Elem.Pkind = Cavity;
        Cav_Alloc(&Cell[i].Elem);
      } else {
        Cell[i].Elem.Pkind = undef;
        std::cerr << "Warning: unsupported PassMethod: " << val << std::endl;
      }
    }

    else if (key.find("Length") != std::string::npos) {
      std::string val;
      getline(iss, val);
      Cell[i].Elem.PL = atof(val.c_str());
    }

    else if (key.find("EApertures") != std::string::npos) {
      std::string x1, x2, y1, y2;
      getline(iss, x1, ',');
      getline(iss, x2, ',');
      getline(iss, y1, ',');
      getline(iss, y2);
      Cell[i].maxampl[X_][0] = atof(x1.c_str());
      Cell[i].maxampl[X_][1] = atof(x2.c_str());
      Cell[i].maxampl[Y_][0] = atof(y1.c_str());
      Cell[i].maxampl[Y_][1] = atof(y2.c_str());
    }

    else if (key.find("PolynomB") != std::string::npos) {
      int n = 0;
      std::string val;
      while (getline(iss, val, ',')) {
        double b = atof(val.c_str());
        if (Cell[i].Elem.Pkind == Mpole && b != 0.0) {
          Cell[i].Elem.M->PB[HOMmax + n] = b;
          Cell[i].Elem.M->PBpar[HOMmax + n] = b;
          Cell[i].Elem.M->Porder = std::max(n, Cell[i].Elem.M->Porder);
        }
        n++;
      }
    }

    else if (key.find("NumIntSteps") != std::string::npos) {
      std::string val;
      getline(iss, val);
      if (Cell[i].Elem.Pkind == Mpole)
        Cell[i].Elem.M->PN = atoi(val.c_str());
    }

    else if (key.find("Voltage") != std::string::npos) {
      std::string val;
      getline(iss, val);
      if (Cell[i].Elem.Pkind == Cavity)
        Cell[i].Elem.C->V_RF = atof(val.c_str());
    }

    else if (key.find("Frequency") != std::string::npos) {
      std::string val;
      getline(iss, val);
      if (Cell[i].Elem.Pkind == Cavity)
        Cell[i].Elem.C->f_RF = atof(val.c_str());
    }

    else if (key.find("HarmonicNumber") != std::string::npos) {
      std::string val;
      getline(iss, val);
      if (Cell[i].Elem.Pkind == Cavity)
        Cell[i].Elem.C->harm_num = atoi(val.c_str());
    }

    else if (key.find("PhaseLag") != std::string::npos) {
      std::string val;
      getline(iss, val);
      if (Cell[i].Elem.Pkind == Cavity)
        Cell[i].Elem.C->phi_RF = atof(val.c_str());
    }
  }

  globval.Cell_nLoc = i;
  std::cout << "rdmfile_new: read " << globval.Cell_nLoc << " elements." << std::endl;
  inf.close();
}
