/* Tracy-3

   J. Bengtsson, BNL 2007

*/

extern tps  sigma_;

bool GetCOD(long imax, double eps, double dP, long &lastpos);

template<typename T>
void Elem_Pass(const long i, ss_vect<T> &x);

template<typename T>
void Cell_Pass(const long i0, const long i1, ss_vect<T> &x, long &lastpos);

void Cell_Pass(const long i0, const long i1, tps &sigma, long &lastpos);

void Cell_SetdP(const double dP);


void Elem_Pass_M(const long i, Vector &xref, Matrix &x);

void Cell_Pass_M(long i0, long i1, Vector &xref, Matrix &mat, long &lastpos);


void Cell_Concat(double dP);

void Cell_fPass(ss_vect<double> &x, long &lastpos);

void Cell_fPass_M(ss_vect<double> &xref, Matrix &mat, long &lastpos);


void Cell_Init(void);
