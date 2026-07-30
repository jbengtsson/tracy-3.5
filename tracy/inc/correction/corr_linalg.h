#ifndef CORRECTION_CORR_LINALG_H
#define CORRECTION_CORR_LINALG_H

// SVD substrate for the response-matrix correctors. A response matrix M
// (m observations x n knobs) is decomposed, singular values below a cut are
// zeroed to suppress ill-conditioned knob combinations, and the knob vector is
// recovered by back-substitution.
//
// Index convention is Numerical-Recipes 1-based throughout: matrices are
// M[1..m][1..n], vectors w[1..n], b[1..m], x[1..n], as returned by
// dmatrix()/dvector(). Callers own the storage.

namespace corr {

// Copy M into the scratch U, SVD-decompose it (M = U diag(w) V^T), then zero
// every singular value below s_cut so svd_backsub drops that direction.
// When prt is set, echo the cut and the singular-value spectrum, n_prt per line.
void svd_decomp_cut(double **M, const int m, const int n, double **U, double *w,
		    double **V, const double s_cut, const bool prt,
		    const int n_prt = 8);

// Solve M x = b in the least-squares / minimum-norm sense from a decomposition
// produced by svd_decomp_cut. Directions whose w was zeroed are excluded.
void svd_backsub(double **U, double *w, double **V, const int m, const int n,
		 double b[], double x[]);

}  // namespace corr

#endif  // CORRECTION_CORR_LINALG_H
