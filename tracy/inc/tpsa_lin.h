 /* Author:	 Johan Bengtsson

   Definitions:  Polymorphic number class.              */

#ifndef TPSA_LIN_H
#define TPSA_LIN_H

void daeps_(const double eps);

void danot_(const long int no);

void daini_(long int no, long int nv, long int fio);

void lieini(const long int no, const long int nv, const long int nd2i);

void daall_(std::vector<double> &x, const long int nd2, const char *daname,
	    const long int no, const long int nv);

void dadal_(std::vector<double> &x, const long int ss_dim);

void davar_(std::vector<double> &x, const double r, const long int i);

void dacon_(std::vector<double> &x, const double r);

void dapek_(const std::vector<double> &x, const long int jj[], double &r);

void dapok_(std::vector<double> &x, const long int jj[], const double r);

double get_m_ij(const ss_vect<tps> &map, const long int i,
		    const long int j);

void put_m_ij(ss_vect<tps> &map, const long int i, const long int j,
	      const double r);

ss_vect<tps> stlmattomap(const std::vector< std::vector<double> > &stlmat);

std::vector< std::vector<double> > maptostlmat(const ss_vect<tps> &map);

arma::mat maptomat(const ss_vect<tps> &map);

ss_vect<tps> mattomap(const arma::mat &mat);

arma::mat stlmattomat(const std::vector< std::vector<double> > &stlmat);

std::vector< std::vector<double> > mattostlmat(const arma::mat &mat);

void dacop_(const std::vector<double> &x, std::vector<double> &z);

void daadd_(const std::vector<double> &x, const std::vector<double> &y,
	    std::vector<double> &z);

void dasub_(const std::vector<double> &x, const std::vector<double> &y,
	    std::vector<double> &z);

void damul_(const std::vector<double> &x, const std::vector<double> &y,
	    std::vector<double> &z);

void dadiv_(const std::vector<double> &x, const std::vector<double> &y,
	    std::vector<double> &z);

void dacad_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dacsu_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dacmu_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dasuc_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dacdi_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dadic_(const std::vector<double> &x, const double y,
	    std::vector<double> &z);

void dapos_(const std::vector<double> &x, std::vector<double> &z);

void dasqr_(const std::vector<double> &x, std::vector<double> &z);

void dacma_(const std::vector<double> &x, const std::vector<double> &y,
	    const double rb, std::vector<double> &z);

void dalin_(const std::vector<double> &x, const double ra,
	    const std::vector<double> &y,
	    const double rb, std::vector<double> &z);

void dainv_(const std::vector<double> &x, std::vector<double> &z);

void dasqrt(const std::vector<double> &x, std::vector<double> &z);

void daexp(const std::vector<double> &x, std::vector<double> &z);

void dalog(const std::vector<double> &x, std::vector<double> &z);

void dasin(const std::vector<double> &x, std::vector<double> &z);

void dacos(const std::vector<double> &x, std::vector<double> &z);

void dasinh(const std::vector<double> &x, std::vector<double> &z);

void dacosh(const std::vector<double> &x, std::vector<double> &z);

void datan(const std::vector<double> &x, std::vector<double> &z);

void daarctan(const std::vector<double> &x, std::vector<double> &z);

void dafun_(const char *fun, const std::vector<double> &x,
	    std::vector<double> &z);

void dacct_(const ss_vect<tps> &x, const long int i,
	    const ss_vect<tps> &y, const long int j, ss_vect<tps> &z,
	    const long int k);

void dainv_(const ss_vect<tps> &x, const long int i, ss_vect<tps> &z,
	    const long int k);

void Rotmap(const long int n, ss_vect<tps> &map, const arma::mat &R);

void daabs_(const std::vector<double> &x, double &r);

void daabs2_(const std::vector<double> &x, double &r);

#endif
