#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

#include "IllegalArgumentException.h"

/** C++ port of:
  *    Jama = Java Matrix class.
  * <P>
  *    The Java Matrix Class provides the fundamental operations of numerical
  *    linear algebra.  Various constructors create Matrices from two dimensional
  *    arrays of double precision floating point numbers.  Various "gets" and
  *    "sets" provide access to submatrices and matrix elements.  Several methods
  *    implement basic matrix arithmetic, including matrix addition and
  *    multiplication, matrix norms, and element-by-element array operations.
  *    Methods for reading and printing matrices are also included.  All the
  *    operations in this version of the Matrix Class involve real matrices.
  *    Complex matrices may be handled in a future version.
  * <P>
  *    Five fundamental matrix decompositions, which consist of pairs or triples
  *    of matrices, permutation vectors, and the like, produce results in five
  *    decomposition classes.  These decompositions are accessed by the Matrix
  *    class to compute solutions of simultaneous linear equations, determinants,
  *    inverses and other matrix functions.  The five decompositions are:
  * <P><UL>
  *    <LI>Cholesky Decomposition of symmetric, positive definite matrices.
  *    <LI>LU Decomposition of rectangular matrices.
  *    <LI>QR Decomposition of rectangular matrices.
  *    <LI>Singular Value Decomposition of rectangular matrices.
  *    <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
  * </UL>
  * <DL>
  * <DT><B>Example of use:</B></DT>
  * <P>
  * <DD>Solve a linear system A x = b and compute the residual norm, ||b - A x||.
  * <P><PRE>
  *       double[][] vals = {{1.,2.,3},{4.,5.,6.},{7.,8.,10.}};
  *       Matrix A = new Matrix(vals);
  *       Matrix b = Matrix.random(3,1);
  *       Matrix x = A.solve(b);
  *       Matrix r = A.times(x).minus(b);
  *       double rnorm = r.normInf();
  * </PRE></DD>
  * </DL>
  * @author The MathWorks, Inc. and the National Institute of Standards and Technology.
  * @version 5 August 1998
  */
class Matrix
{
public:
	/** Construct an m-by-n matrix of zeros.
	  * @param m    Number of rows.
	  * @param n    Number of colums.
	  */
	Matrix(int m = 0, int n = 0);
	/** Construct an m-by-n constant matrix.
	  * @param m    Number of rows.
	  * @param n    Number of colums.
	  * @param s    Fill the matrix with this scalar value.
	  */
	Matrix(int m, int n, double s);
	/** Construct a matrix from a 2-D array.
	  * @param A    Two-dimensional array of doubles.
	  * @exception  IllegalArgumentException All rows must have the same length
	  * @see        #constructWithCopy
	  */
	Matrix(const std::vector<std::vector<double>>& A);
	/** Construct a matrix quickly without checking arguments.
	  * @param A    Two-dimensional array of doubles.
	  * @param m    Number of rows.
	  * @param n    Number of colums.
	  */
	Matrix(const std::vector<std::vector<double>>& A, int m, int n);
	/** Construct a matrix from a one-dimensional packed array
	  * @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
	  * @param m    Number of rows.
	  * @exception  IllegalArgumentException Array length must be a multiple of m.
	  */
	Matrix(const std::vector<double>& vals, int m);
	
	double operator()(int i, int j) const;
	
	/** Construct a matrix from a copy of a 2-D array.
	  * @param A    Two-dimensional array of doubles.
	  * @exception  IllegalArgumentException All rows must have the same length
	  */
	static Matrix constructWithCopy(const std::vector<std::vector<double>>& A);
	/** Make a deep copy of a matrix */
	Matrix copy();
	/** Clone the Matrix object. */
	Matrix clone();
	/** Access the internal two-dimensional array.
	  * @return     Pointer to the two-dimensional array of matrix elements.
	  */
	std::vector<std::vector<double>>& getArray();
	/** Copy the internal two-dimensional array.
	  * @return     Two-dimensional array copy of matrix elements.
	  */
	std::vector<std::vector<double>> getArrayCopy() const;
	/** Make a one-dimensional column packed copy of the internal array.
	  * @return     Matrix elements packed in a one-dimensional array by columns.
	  */
	std::vector<double> getColumnPackedCopy() const;
	/** Make a one-dimensional row packed copy of the internal array.
	  * @return     Matrix elements packed in a one-dimensional array by rows.
	  */
	std::vector<double> getRowPackedCopy() const;
	/** Get row dimension.
	  * @return     m, the number of rows.
	  */
	int getRowDimension() const;
	/** Get column dimension.
	  * @return     n, the number of columns.
	  */
	int getColumnDimension() const;
	/** Get a single element.
	  * @param i    Row index.
	  * @param j    Column index.
	  * @return     A(i,j)
	  * @exception  ArrayIndexOutOfBoundsException
	  */
	double get(int i, int j) const;
	/** Get a submatrix.
	  * @param i0   Initial row index
	  * @param i1   Final row index
	  * @param j0   Initial column index
	  * @param j1   Final column index
	  * @return     A(i0:i1,j0:j1)
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	Matrix getMatrix(int i0, int i1, int j0, int j1);
	/** Get a submatrix.
	  * @param r    Array of row indices.
	  * @param c    Array of column indices.
	  * @return     A(r(:),c(:))
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	Matrix getMatrix(const std::vector<int>& r, const std::vector<int>& c);
	/** Get a submatrix.
	  * @param i0   Initial row index
	  * @param i1   Final row index
	  * @param c    Array of column indices.
	  * @return     A(i0:i1,c(:))
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	Matrix getMatrix(int i0, int i1, const std::vector<int>& c);
	/** Get a submatrix.
	  * @param r    Array of row indices.
	  * @param j0   Initial column index
	  * @param j1   Final column index
	  * @return     A(r(:),j0:j1)
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	Matrix getMatrix(const std::vector<int>& r, int j0, int j1);
	/** Set a single element.
	  * @param i    Row index.
	  * @param j    Column index.
	  * @param s    A(i,j).
	  * @exception  ArrayIndexOutOfBoundsException
	  */
	void set(int i, int j, double s);
	/** Set a submatrix.
	  * @param i0   Initial row index
	  * @param i1   Final row index
	  * @param j0   Initial column index
	  * @param j1   Final column index
	  * @param X    A(i0:i1,j0:j1)
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	void setMatrix(int i0, int i1, int j0, int j1, const Matrix& X);
	/** Set a submatrix.
	  * @param r    Array of row indices.
	  * @param c    Array of column indices.
	  * @param X    A(r(:),c(:))
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	void setMatrix(const std::vector<int>& r, const std::vector<int>& c, const Matrix& X);
	/** Set a submatrix.
	  * @param r    Array of row indices.
	  * @param j0   Initial column index
	  * @param j1   Final column index
	  * @param X    A(r(:),j0:j1)
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	void setMatrix(const std::vector<int>& r, int j0, int j1, const Matrix& X);
	/** Set a submatrix.
	  * @param i0   Initial row index
	  * @param i1   Final row index
	  * @param c    Array of column indices.
	  * @param X    A(i0:i1,c(:))
	  * @exception  ArrayIndexOutOfBoundsException Submatrix indices
	  */
	void setMatrix(int i0, int i1, const std::vector<int>& c, const Matrix& X);
	/** Matrix transpose.
	  * @return    A'
	  */
	Matrix transpose();
	/** One norm
	  * @return    maximum column sum.
	  */
	double norm1() const;
	/** Two norm
	  * @return    maximum singular value.
	  */
	double norm2() const;
	/** Infinity norm
	  * @return    maximum row sum.
	  */
	double normInf() const;
	/** Frobenius norm
	  * @return    sqrt of sum of squares of all elements.
	  */
	double normF() const;
	/**  Unary minus
	  * @return    -A
	  */
	Matrix uminus();
	/** C = A + B
	  * @param B    another matrix
	  * @return     A + B
	  */
	Matrix plus(const Matrix& B) const;
	/** A = A + B
	  * @param B    another matrix
	  * @return     A + B
	  */
	Matrix& plusEquals(const Matrix& B);
	/** C = A - B
	  * @param B    another matrix
	  * @return     A - B
	  */
	Matrix minus(const Matrix& B) const;
	/** A = A - B
	  * @param B    another matrix
	  * @return     A - B
	  */
	Matrix& minusEquals(const Matrix& B);
	/** Element-by-element multiplication, C = A.*B
	  * @param B    another matrix
	  * @return     A.*B
	  */
	Matrix arrayTimes(const Matrix& B) const;
	/** Element-by-element multiplication in place, A = A.*B
	  * @param B    another matrix
	  * @return     A.*B
	  */
	Matrix& arrayTimesEquals(const Matrix& B);
	/** Element-by-element right division, C = A./B
	  * @param B    another matrix
	  * @return     A./B
	  */
	Matrix arrayRightDivide(const Matrix& B) const;
	/** Element-by-element right division in place, A = A./B
	  * @param B    another matrix
	  * @return     A./B
	  */
	Matrix& arrayRightDivideEquals(const Matrix& B);
	/** Element-by-element left division, C = A.\B
	  * @param B    another matrix
	  * @return     A.\B
	  */
	Matrix arrayLeftDivide(const Matrix& B) const;
	/** Element-by-element left division in place, A = A.\B
	  * @param B    another matrix
	  * @return     A.\B
	  */
	Matrix& arrayLeftDivideEquals(const Matrix& B);
	/** Multiply a matrix by a scalar, C = s*A
	  * @param s    scalar
	  * @return     s*A
	  */
	Matrix times(double s) const;
	/** Multiply a matrix by a scalar in place, A = s*A
	  * @param s    scalar
	  * @return     replace A by s*A
	  */
	Matrix& timesEquals(double s);
	/** Linear algebraic matrix multiplication, A * B
	  * @param B    another matrix
	  * @return     Matrix product, A * B
	  * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	  */
	Matrix times(const Matrix& B) const;
	/**
	  * Linear algebraic matrix multiplication, A * B
	  * B being a triangular matrix
	  * <b>Note:</b>
	  * Actually the matrix should be a <b>column orienten, upper triangular
	  * matrix</b> but use the <b>row oriented, lower triangular matrix</b>
	  * instead (transposed), because this is faster due to the easyer array
	  * access.
	  *
	  * @param B    another matrix
	  * @return     Matrix product, A * B
	  *
	  * @exception  IllegalArgumentException Matrix inner dimensions must agree.
	  */
	Matrix timesTriangular(Matrix& B) const;
	/**
	  * X.diffEquals() calculates differences between adjacent columns of this
	  * matrix. Consequently the size of the matrix is reduced by one. The result
	  * is stored in this matrix object again.
	  */
	void diffEquals();
	/**
	  * X.logEquals() calculates the natural logarithem of each element of the
	  * matrix. The result is stored in this matrix object again.
	  */
	void logEquals();
	/**
	  * X.powEquals() calculates the power of each element of the matrix. The
	  * result is stored in this matrix object again.
	  */
	void powEquals(double exp);
	/**
	  * X.powEquals() calculates the power of each element of the matrix.
	  *
	  * @return Matrix
	  */
	Matrix pow(double exp) const;
	/**
	  * X.thrunkAtLowerBoundariy(). All values smaller than the given one are set
	  * to this lower boundary.
	  */
	void thrunkAtLowerBoundary(double value);
	/** Matrix trace.
	  * @return     sum of the diagonal elements.
	  */
	double trace() const;
	/** Generate matrix with random elements
	  * @param m    Number of rows.
	  * @param n    Number of colums.
	  * @return     An m-by-n matrix with uniformly distributed random elements.
	  */
	static Matrix random(int m, int n);
	/** Generate identity matrix
	  * @param m    Number of rows.
	  * @param n    Number of colums.
	  * @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
	  */
	static Matrix identity(int m, int n);
	/** Print the matrix to stdout.   Line the elements up in columns
	  * with a Fortran-like 'Fw.d' style format.
	  * @param w    Column width.
	  * @param d    Number of digits after the decimal.
	  */
	void print(int w, int d) const;
	static Matrix readCSV(const char* file);
	/**
	  * Returns the mean values along the specified dimension.
	  *
	  * @param dim
	  *            If 1, then the mean of each column is returned in a row
	  *            vector. If 2, then the mean of each row is returned in a
	  *            column vector.
	  * @return A vector containing the mean values along the specified
	  *         dimension.
	  */
	Matrix mean(int dim) const;
	/**
	  * Calculate the full covariance matrix.
	  * @return the covariance matrix
	  */
	Matrix cov();
	/**
	  * Returns the sum of the component of the matrix.
	  * @return the sum
	  */
	double sum() const;
	/**
	  * returns a new Matrix object, where each value is set to the absolute value
	  * @return a new Matrix with all values being positive
	  */
	Matrix abs() const;
	/**
	  * Writes the Matrix to an ascii-textfile that can be read by Matlab.
	  * Usage in Matlab: load('filename', '-ascii');
	  * @param filename the name of the ascii file to create, e.g. "C:\\temp\\matrix.ascii"
	  * @throws IllegalArgumentException if there is a problem with the filename
	  */
	void writeAscii(const char* filename) const;
	
private:
	friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
	friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
	friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
	friend Matrix& operator+=(Matrix& lhs, const Matrix& rhs); 
	friend Matrix& operator-=(Matrix& lhs, const Matrix& rhs);
	
	/** Check if size(A) == size(B) **/
	void checkMatrixDimensions(const Matrix& B) const;
	double cov(const std::vector<double>& vec1, const std::vector<double>& vec2);
	/**
	  * the mean of the values in the double array
	  * @param vec	double values
	  * @return	the mean of the values in vec
	  */
	double mean(const std::vector<double>& vec);
	
	/** Array for internal storage of elements.
	  * @serial internal array storage.
	  */
	std::vector<std::vector<double>> A;
	int m; /** Number of rows. @serial number of rows. */
	int n; /** Number of columns. @serial number of columns. */
};

#endif /* MATRIX_H */