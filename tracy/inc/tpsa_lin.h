 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

void daeps_(const double eps);

void danot_(const int no);

void daini_(int no, int nv, int fio);

void lieini(const int no, const int nv, const int nd2i);

void daall_(tps_buf &x, const int nd2, const char *daname,
	    const int no, const int nv);

void dadal_(tps_buf &x, const int ss_dim);

void davar_(tps_buf &x, const double r, const int i);

void dacon_(tps_buf &x, const double r);

void dapek_(const tps_buf &x, const int jj[], double &r);

void dapok_(tps_buf &x, const int jj[], const double r);

double getmat(const ss_vect<tps> &map, const int i, const int j);

void putmat(ss_vect<tps> &map, const int i, const int j, const double r);

void getlinmat(const int nv, const ss_vect<tps> &map, Matrix &mat);

void putlinmat(const int nv, const Matrix &mat, ss_vect<tps> &map);

void dacop_(const tps_buf &x, tps_buf &z);

void daadd_(const tps_buf &x, const tps_buf &y, tps_buf &z);

void dasub_(const tps_buf &x, const tps_buf &y, tps_buf &z);

void damul_(const tps_buf &x, const tps_buf &y, tps_buf &z);

void dadiv_(const tps_buf &x, const tps_buf &y, tps_buf &z);

void dacad_(const tps_buf &x, const double y, tps_buf &z);

void dacsu_(const tps_buf &x, const double y, tps_buf &z);

void dacmu_(const tps_buf &x, const double y, tps_buf &z);

void dasuc_(const tps_buf &x, const double y, tps_buf &z);

void dacdi_(const tps_buf &x, const double y, tps_buf &z);

void dadic_(const tps_buf &x, const double y, tps_buf &z);

void dapos_(const tps_buf &x, tps_buf &z);

void dasqr_(const tps_buf &x, tps_buf &z);

void dacma_(const tps_buf &x, const tps_buf &y, const double rb, tps_buf &z);

void dalin_(const tps_buf &x, const double ra, const tps_buf &y,
	    const double rb, tps_buf &z);

void dainv_(const tps_buf &x, tps_buf &z);

void dasqrt(const tps_buf &x, tps_buf &z);

void daexp(const tps_buf &x, tps_buf &z);

void dalog(const tps_buf &x, tps_buf &z);

void dasin(const tps_buf &x, tps_buf &z);

void dacos(const tps_buf &x, tps_buf &z);

void dasinh(const tps_buf &x, tps_buf &z);

void dacosh(const tps_buf &x, tps_buf &z);

void datan(const tps_buf &x, tps_buf &z);

void daarctan(const tps_buf &x, tps_buf &z);

void dafun_(const char *fun, const tps_buf &x, tps_buf &z);

void dacct_(const ss_vect<tps> &x, const int i,
	    const ss_vect<tps> &y, const int j, ss_vect<tps> &z, const int k);

void dainv_(const ss_vect<tps> &x, const int i, ss_vect<tps> &z, const int k);

void Rotmap(const int n, ss_vect<tps> &map, const Matrix &R);

void daabs_(const tps_buf &x, double &r);

void daabs2_(const tps_buf &x, double &r);
