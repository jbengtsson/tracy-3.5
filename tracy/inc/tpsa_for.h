/* Author:	 Johan Bengtsson.

   Definitions:  Interface to Fortran TPSA- and Lie library.  */


extern "C" {
  // Interface to FORTRAN TPSA-library
  void daini_(const int &, const int &, const int &);
  void daall_(int &, const int &, const char [], 
	      const int &, const int &);
  void dadal_(int &, const int &);
  void dacon_(int &, const double &);
  void davar_(int &, const double &, const int &);
  void daadd_(const int &, const int &, int &);
  void dasub_(const int &, const int &, int &);
  void damul_(const int &, const int &, int &);
  void dadiv_(const int &, const int &, int &);
  void dacad_(const int &, const double &, int &);
  // dacsu <=> dacad(a, -r, b)
  //void dasuc_(const int &, const double &, int &);
  void dacmu_(const int &, const double &, int &);
  // dacdi <=> dacmu(a, 1/r, b)
  //void dadic_(const int &, const double &, int &);
  void dafun_(const char [], const int &, int &);
  void dader_(const int &, const int &, int &);
  void danot_(const int &);
  void daeps_(const double &);
  void dacop_(const int &, int &);
  void daabs_(const int &, double &);
  void daabs2_(const int &, double &);
  void dapek_(const int &, const int [], double &);
  void dapok_(int &, const int [], const double &);
  void dapri_(const int &, const int &);
  void darea_(const int &, const int &);
  void dapoi_(const int &, const int &, int &, const int &);
  void hash_(const int &, const int &, const int [], int &, int &);
  void dehash_(const int &, const int &, const int &, const int &, int []);
  void daimp_(const double [], const int [], const int [], int &);
  void dainv_(const int [], const int &, int [], const int &);
//  void dapin_(const int [], const int &, int [], const int &, const int []);
  void daexp_(const int &, double [], int [], int [], char []); 
  void dacct_(const int [], const int &, const int [], const int &,
	      int [], const int &);

  // Interface to FORTRAN Lie-lib
  void lieinit_(const int &, const int &, const int &,
		const int &, const int &, const int &);
  void idprset_(const int &);
  void exp1d_(const int &, const int &, int &,
	      const double &, const int &);
  void daflo_(const int [], const int &, int &);
  void expnd2_(const int &, const int [], int [],
	       const double &, const int &);
  void daflod_(const int [], const int [], int []);
  void etinv_(const int [], int []);
  void etpin_(const int [], int [], const int []);
  void etcct_(const int [], const int [], int []);
  bool mapnorm_(const int [], int &, int [], int [], int [],
		int &, const int &);
  bool mapnormf_(const int [], int [], int [], int [], int [],
		 int [], const int &, const int &);
  void gofix_(const int [], int [], int [], const int &);
  void dhdj_(const int &, int []);
  void ctor_(const int &, int &, int &);
  void rtoc_(const int &, const int &, int &);
  void fexpo_(const int &, const int [], int [], const int &, const int &,
	      const double &, const int &);
  void liefact_(const int [], int [], int &);
  void flofacg_(const int [], int [], const double &);
  void intd_(const int [], int &, double &);
  void difd_(const int &, int [], double &);
  void etpoi_(const int &, const int &, int &);
  void take_(const int &, const int &, int &);
  void taked_(const int [], const int &, int []);
  void gettura_(double [], double []);
  void etmtree_(const int [], int []);
  void etppush2_(const int [], const ss_vect<double> &, ss_vect<double> &);
}
