#ifndef SINGULARVALUEDECOMPOSITION_H
#define SINGULARVALUEDECOMPOSITION_H

#include <vector>

#include "Matrix.h"

/** Singular Value Decomposition.
  * <P>
  * For an m-by-n matrix A with m >= n, the singular value decomposition is
  * an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
  * an n-by-n orthogonal matrix V so that A = U*S*V'.
  * <P>
  * The singular values, sigma[k] = S[k][k], are ordered so that
  * sigma[0] >= sigma[1] >= ... >= sigma[n-1].
  * <P>
  * The singular value decomposition always exists, so the constructor will
  * never fail.  The matrix condition number and the effective numerical
  * rank can be computed from this decomposition.
  */
class SingularValueDecomposition
{
public:
	/** Construct the singular value decomposition
      * @param Arg    Rectangular matrix
      */
	explicit SingularValueDecomposition(const Matrix& Arg);
	
	/** Return the left singular vectors
	  * @return     U
	  */
	Matrix getU();
	/** Return the right singular vectors
	  * @return     V
	  */
	Matrix getV();
	/** Return the one-dimensional array of singular values
	  * @return     diagonal of S.
	  */
	std::vector<double> getSingularValues() const;
	/** Return the diagonal matrix of singular values
	  * @return     S
	  */
	Matrix getS();
	/** Two norm
	  * @return     max(S)
	  */
	double norm2() const;
	/** Two norm condition number
	  * @return     max(S)/min(S)
	  */
	double cond() const;
	/** Effective numerical matrix rank
	  * @return     Number of nonnegligible singular values.
	  */
	int rank() const;
	
private:
	/** Arrays for internal storage of U.
	  * @serial internal storage of U.
	  */
	std::vector<std::vector<double>> U;
	/** Arrays for internal storage of V.
	  * @serial internal storage of V.
	  */
	std::vector<std::vector<double>> V;
	
	/** Array for internal storage of singular values.
	  * @serial internal storage of singular values.
	  */
	std::vector<double> s;
	
	/** Row dimensions.
	  * @serial row dimension.
	  */
	int m;
	/** Column dimensions.
	  * @serial column dimension.
	  */
	int n;
};

#endif /* SINGULARVALUEDECOMPOSITION_H */