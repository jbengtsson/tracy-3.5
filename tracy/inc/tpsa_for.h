/* Author:	 Johan Bengtsson.

   Definitions:  Interface to Fortran TPSA- and Lie library.  */

#ifndef TPSA_FOR_H
#define TPSA_FOR_H

extern "C" {
  // Interface to FORTRAN TPSA-library
  void daini_(const long int &, const long int &, const long int &);
  void daall_(long int &, const long int &, const char [], const long int &,
	      const long int &);
  void dadal_(long int &, const long int &);
  void dacon_(long int &, const double &);
  void davar_(long int &, const double &, const long int &);
  void daadd_(const long int &, const long int &, long int &);
  void dasub_(const long int &, const long int &, long int &);
  void damul_(const long int &, const long int &, long int &);
  void dadiv_(const long int &, const long int &, long int &);
  void dacad_(const long int &, const double &, long int &);
  // dacsu <=> dacad(a, -r, b)
  /* void dasuc_(const long int &, const double &, long int &); */
  void dacmu_(const long int &, const double &, long int &);
  // dacdi <=> dacmu(a, 1/r, b)
  /* void dadic_(const long int &, const double &, long int &); */
  void dafun_(const char [], const long int &, long int &);
  void dader_(const long int &, const long int &, long int &);
  void danot_(const long int &);
  long int getno_(void);
  void daeps_(const double &);
  void dacop_(const long int &, long int &);
  void daabs_(const long int &, double &);
  void daabs2_(const long int &, double &);
  void dapek_(const long int &, const long int [], double &);
  void dapok_(long int &, const long int [], const double &);
  void dapri_(const long int &, const long int &);
  void darea_(const long int &, const long int &);
  void dapoi_(const long int &, const long int &, long int &, const long int &);
  void hash_(const long int &, const long int &, const long int [], long int &,
	     long int &);
  void dehash_(const long int &, const long int &, const long int &,
	       const long int &, long int []);
  void daimp_(const double [], const long int [], const long int [],
	      long int &);
  void dainv_(const long int [], const long int &, long int [],
	      const long int &);
  /* void dapin_(const long int [], const long int &, long int [], */
  /*  	      const long int &, const long int []); */
  void daexp_(const long int &, double [], long int [], long int [], char []); 
  void dacct_(const long int [], const long int &, const long int [],
	      const long int &, long int [], const long int &);

  // Long Interface to FORTRAN Lie-lib
  void lieinit_(const long int &, const long int &, const long int &,
		const long int &, const long int &, const long int &);
  void idprset_(const long int &);
  void exp1d_(const long int &, const long int &, long int &,
	      const double &, const long int &);
  void daflo_(const long int [], const long int &, long int &);
  void expnd2_(const long int &, const long int [], long int [],
	       const double &, const long int &);
  void daflod_(const long int [], const long int [], long int []);
  void etinv_(const long int [], long int []);
  void etpin_(const long int [], long int [], const long int []);
  void etcct_(const long int [], const long int [], long int []);
  bool mapnorm_(const long int [], long int &, long int [], long int [],
		long int [], long int &, const long int &);
  bool mapnormf_(const long int [], long int [], long int [], long int [],
		 long int [], long int [], const long int &, const long int &);
  void gofix_(const long int [], long int [], long int [], const long int &);
  void dhdj_(const long int &, long int []);
  void ctor_(const long int &, long int &, long int &);
  void ctoi_(const long int &, long int &);
  void cpart_(const long int &, long int &);
  void etctr_(const long int []);
  void etcjg_(long int []);
  void trx_(const long int &, long int &, long int []);
  void rtoc_(const long int &, const long int &, long int &);
  void fexpo_(const long int &, const long int [], long int [],
	      const long int &, const long int &, const double &,
	      const long int &);
  void liefact_(const long int [], long int [], long int &);
  void flofacg_(const long int [], long int [], const double &);
  void intd_(const long int [], long int &, double &);
  void difd_(const long int &, long int [], double &);
  void etpoi_(const long int &, const long int &, long int &);
  void take_(const long int &, const long int &, long int &);
  void taked_(const long int [], const long int &, long int []);
  void gettura_(double [], double []);
  void etmtree_(const long int [], long int []);
  void etppush2_(const long int [], const ss_vect<double> &, ss_vect<double> &);
}

#endif
