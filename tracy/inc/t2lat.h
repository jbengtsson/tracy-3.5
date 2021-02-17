/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/

#ifndef T2LAT_H
#define T2LAT_H

#define Cell_nLocMax 20000  // maximum number of LEGO blocks (Cell_nLoc).
#define DBNameLen    39

typedef char   DBNameType[DBNameLen];

long int ElemIndex(const std::string &name1);

#endif
