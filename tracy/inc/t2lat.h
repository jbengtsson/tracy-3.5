/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


// maximum number of LEGO blocks (Cell_nLoc)
#define Cell_nLocMax 20000

// maximum number of families for Elem_NFam
#define Elem_nFamMax 3000


class Lattice_Type {
 private:
 public:
  ElemFamType ElemFam[Elem_nFamMax];
  CellType    Cell[Cell_nLocMax+1];

  bool Lattice_Read(FILE **fi_, FILE **fo_);
};

long ElemIndex(const std::string &name1);
void Read_Lattice(const char *fic);
