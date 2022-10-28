/* Author:        Johan Bengtsson

    Definitions:  Interface to Fortran library for Truncated Power
		  Series Algebra.

    Note, linear case is a special case, see e.g. daexp

*/

#ifndef TPSA_FOR_PM_H
#define TPSA_FOR_PM_H

// Name length for FORTRAN library is 10; 10+1 for C.
const int name_len_for = 10; 

extern int  bufsize;  // Note, max no of monomials is (no+nv)!/(nv!*no!)
 

long int fact(long int n);

long int nok(long int n, long int k);

double getmat(const ss_vect<tps> &map, const int i, const int j);

void putmat(ss_vect<tps> &map, const int i, const int j, const double r);

void getlinmat(const int nv, const ss_vect<tps> &map, Matrix &mat);

ss_vect<tps> putlinmat(const int nv, const Matrix &mat);

void idprset(const int level);

tps atan2(const tps &b,const tps &a);

#endif
