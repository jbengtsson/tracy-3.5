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

void Cell_Init(void);
