#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>

#include "SingularValueDecomposition.h"

Matrix::Matrix(int m /* = 0 */, int n /* = 0 */)
{
	this->m = m;
	this->n = n;
	A.resize(m, std::vector<double>(n, 0.0));
}

Matrix::Matrix(int m, int n, double s)
{
	this->m = m;
	this->n = n;
	A.resize(m, std::vector<double>(n, s));
}

Matrix::Matrix(const std::vector<std::vector<double>>& A)
{
	m = static_cast<int>(A.size());
	n = static_cast<int>(A[0].size());
	for (int i = 0; i < m; ++i)
	{
		if (static_cast<int>(A[i].size()) != n)
			throw IllegalArgumentException("All rows must have the same length.");
	}
	this->A = A;
}

Matrix::Matrix(const std::vector<std::vector<double>>& A, int m, int n)
{
	this->A = A;
	this->m = m;
	this->n = n;
}

Matrix::Matrix(const std::vector<double>& vals, int m)
{
	this->m = m;
	n = (m != 0 ? static_cast<int>(vals.size()) / m : 0);
	if (m * n != static_cast<int>(vals.size()))
		throw IllegalArgumentException("Array length must be a multiple of m.");
	A.resize(m, std::vector<double>(n, 0.0));
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = vals[i + j * m];
	}
}

double Matrix::operator()(int i, int j) const
{
	return A[i][j];
}

Matrix Matrix::constructWithCopy(const std::vector<std::vector<double>>& A)
{
	int m = static_cast<int>(A.size());
	int n = static_cast<int>(A[0].size());
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		if (static_cast<int>(A[i].size()) != n)
			throw IllegalArgumentException("All rows must have the same length.");
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j];
	}
	return X;
}

Matrix Matrix::copy()
{
	Matrix X {m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j];
	}
	return X;
}

Matrix Matrix::clone()
{
	return copy();
}

std::vector<std::vector<double>>& Matrix::getArray()
{
	return A;
}

std::vector<std::vector<double>> Matrix::getArrayCopy() const
{
	std::vector<std::vector<double>> C;
	C.resize(m, std::vector<double>(n, 0.0));
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j];
	}
	return C;
}

std::vector<double> Matrix::getColumnPackedCopy() const
{
	std::vector<double> vals(m * n, 0.0);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			vals[i + j * m] = A[i][j];
	}
	return vals;
}

std::vector<double> Matrix::getRowPackedCopy() const
{
	std::vector<double> vals(m * n, 0.0);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			vals[i * n + j] = A[i][j];
	}
	return vals;
}

int Matrix::getRowDimension() const
{
	return m;
}

int Matrix::getColumnDimension() const
{
	return n;
}

double Matrix::get(int i, int j) const
{
	return A[i][j];
}

Matrix Matrix::getMatrix(int i0, int i1, int j0, int j1)
{
	Matrix X {i1 - i0 + 1, j1 - j0 + 1};
	std::vector<std::vector<double>>& B = X.getArray();
	try
	{
		for (int i = i0; i <= i1; ++i)
		{
			for (int j = j0; j <= j1; ++j)
				B[i - i0][j - j0] = A[i][j];
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
	return X;
}

Matrix Matrix::getMatrix(const std::vector<int>& r, const std::vector<int>& c)
{
	Matrix X {static_cast<int>(r.size()), static_cast<int>(c.size())};
	std::vector<std::vector<double>>& B = X.getArray();
	try
	{
		for (size_t i = 0; i < r.size(); ++i)
		{
			for (size_t j = 0; j < c.size(); ++j)
				B[i][j] = A[r[i]][c[j]];
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
	return X;
}

Matrix Matrix::getMatrix(int i0, int i1, const std::vector<int>& c)
{
	Matrix X {i1 - i0 + 1, static_cast<int>(c.size())};
	std::vector<std::vector<double>>& B = X.getArray();
	try
	{
		for (int i = i0; i <= i1; ++i)
		{
			for (size_t j = 0; j < c.size(); ++j)
				B[i - i0][j] = A[i][c[j]];
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
	return X;
}

Matrix Matrix::getMatrix(const std::vector<int>& r, int j0, int j1)
{
	Matrix X {static_cast<int>(r.size()), j1 - j0 + 1};
	std::vector<std::vector<double>>& B = X.getArray();
	try
	{
		for (size_t i = 0; i < r.size(); ++i)
		{
			for (int j = j0; j <= j1; ++j)
				B[i][j - j0] = A[r[i]][j];
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
	return X;
}

void Matrix::set(int i, int j, double s)
{
	A[i][j] = s;
}

void Matrix::setMatrix(int i0, int i1, int j0, int j1, const Matrix& X)
{
	try
	{
		for (int i = i0; i <= i1; ++i)
		{
			for (int j = j0; j <= j1; ++j)
				A[i][j] = X.get(i - i0, j - j0);
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
}

void Matrix::setMatrix(const std::vector<int>& r, const std::vector<int>& c, const Matrix& X)
{
	try
	{
		for (int i = 0; i < static_cast<int>(r.size()); ++i)
		{
			for (int j = 0; j < static_cast<int>(c.size()); ++j)
				A[r[i]][c[j]] = X.get(i, j);
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
}

void Matrix::setMatrix(const std::vector<int>& r, int j0, int j1, const Matrix& X)
{
	try
	{
		for (int i = 0; i < static_cast<int>(r.size()); ++i)
		{
			for (int j = j0; j <= j1; ++j)
				A[r[i]][j] = X.get(i, j - j0);
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
}

void Matrix::setMatrix(int i0, int i1, const std::vector<int>& c, const Matrix& X)
{
	try
	{
		for (int i = i0; i <= i1; ++i)
		{
			for (int j = 0; j < static_cast<int>(c.size()); ++j)
				A[i][c[j]] = X.get(i - i0, j);
		}
	}
	catch(...)
	{
		throw std::out_of_range("Submatrix indices");
	}
}

Matrix Matrix::transpose()
{
	Matrix X{n, m};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[j][i] = A[i][j];
	}
	return X;
}

double Matrix::norm1() const
{
	double f = 0;
	for (int j = 0; j < n; ++j)
	{
		double s = 0;
		for (int i = 0; i < m; ++i)
			s += std::abs(A[i][j]);
		f = std::max(f, s);
	}
	return f;
}

double Matrix::norm2() const
{
	SingularValueDecomposition svd(*this);
	return svd.norm2();
}

double Matrix::normInf() const
{
	double f = 0;
	for (int i = 0; i < m; ++i)
	{
		double s = 0;
		for (int j = 0; j < n; ++j)
			s += std::abs(A[i][j]);
		f = std::max(f, s);
	}
	return f;
}

double Matrix::normF() const
{
	double f = 0;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			f = std::hypot(f, A[i][j]);
	}
	return f;
}

Matrix Matrix::uminus()
{
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = -A[i][j];
	}
	return X;
}

Matrix Matrix::plus(const Matrix& B) const
{
	checkMatrixDimensions(B);
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j] + B.A[i][j];
	}
	return X;
}

Matrix& Matrix::plusEquals(const Matrix& B)
{
	checkMatrixDimensions(B);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = A[i][j] + B.A[i][j];
	}
	return *this;
}

Matrix Matrix::minus(const Matrix& B) const
{
	checkMatrixDimensions(B);
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j] - B.A[i][j];
	}
	return X;
}

Matrix& Matrix::minusEquals(const Matrix& B)
{
	checkMatrixDimensions(B);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = A[i][j] - B.A[i][j];
	}
	return *this;
}

Matrix Matrix::arrayTimes(const Matrix& B) const
{
	checkMatrixDimensions(B);
	Matrix X {m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j] * B.A[i][j];
	}
	return X;
}

Matrix& Matrix::arrayTimesEquals(const Matrix& B)
{
	checkMatrixDimensions(B);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = A[i][j] * B.A[i][j];
	}
	return *this;
}

Matrix Matrix::arrayRightDivide(const Matrix& B) const
{
	checkMatrixDimensions(B);
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = A[i][j] / B.A[i][j];
	}
	return X;
}

Matrix& Matrix::arrayRightDivideEquals(const Matrix& B)
{
	checkMatrixDimensions(B);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = A[i][j] / B.A[i][j];
	}
	return *this;
}

Matrix Matrix::arrayLeftDivide(const Matrix& B) const
{
	checkMatrixDimensions(B);
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = B.A[i][j] / A[i][j];
	}
	return X;
}

Matrix& Matrix::arrayLeftDivideEquals(const Matrix& B)
{
	checkMatrixDimensions(B);
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = B.A[i][j] / A[i][j];
	}
	return *this;
}

Matrix Matrix::times(double s) const
{
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = s * A[i][j];
	}
	return X;
}

Matrix& Matrix::timesEquals(double s)
{
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			A[i][j] = s * A[i][j];
	}
	return *this;
}

Matrix Matrix::times(const Matrix& B) const
{
	if (B.m != n)
		throw IllegalArgumentException("Matrix inner dimensions must agree.");
	Matrix X{m, B.n};
	std::vector<std::vector<double>>& C = X.getArray();
	std::vector<double> Bcolj(n, 0.0);
	for (int j = 0; j < B.n; ++j)
	{
		for (int k = 0; k < n; ++k)
			Bcolj[k] = B.A[k][j];
		for (int i = 0; i < m; ++i)
		{
			std::vector<double> Arowi = A[i];
			double s = 0;
			for (int k = 0; k < n; ++k)
				s += Arowi[k] * Bcolj[k];
			C[i][j] = s;
		}
	}
	return X;
}

Matrix Matrix::timesTriangular(Matrix& B) const
{
	if (B.m != n)
		throw IllegalArgumentException("Matrix inner dimensions must agree.");
	
	Matrix X{m, B.n};
	std::vector<std::vector<double>>& c = X.getArray();
	std::vector<std::vector<double>>& b = B.getArray();
	double s = 0;
	std::vector<double> Arowi;
	std::vector<double> Browj;
	
	// Multiply with each row of A
	for (int i = 0; i < m; ++i)
	{
		Arowi = A[i];
		
		// For all columns of B
		for (int j = 0; j < B.n; ++j)
		{
			s = 0;
			Browj = b[j];
			// Since B being triangular, this loop uses k <= j
			for (int k = 0; k <= j; ++k)
				s += Arowi[k] * Browj[k];
			c[i][j] = s;
		}
	}
	return X;
}

void Matrix::diffEquals()
{
	std::vector<double> col;
	for(size_t i = 0; i < A.size(); ++i)
	{
		col = std::vector<double>(A[i].size() - 1);
		
		for(int j = 1; j < static_cast<int>(A[i].size()); ++j)
			col[j - 1] = std::abs(A[i][j] - A[i][j - 1]);
		
		A[i] = col;
	}
	n--;
}

void Matrix::logEquals()
{
	for(size_t i = 0; i < A.size(); ++i)
	{
		for(size_t j = 0; j < A[i].size(); ++j)
			A[i][j] = std::log(A[i][j]);
	}
}

void Matrix::powEquals(double exp)
{
	for(size_t i = 0; i < A.size(); ++i)
	{
		for(size_t j = 0; j < A[i].size(); ++j)
			A[i][j] = std::pow(A[i][j], exp);
	}
}

Matrix Matrix::pow(double exp) const
{
	Matrix X{m, n};
	std::vector<std::vector<double>>& C = X.getArray();
	
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			C[i][j] = std::pow(A[i][j], exp);
	}
	
	return X;
}

void Matrix::thrunkAtLowerBoundary(double value)
{
	for(size_t i = 0; i < A.size(); ++i)
	{
		for(size_t j = 0; j < A[i].size(); ++j)
		{
			if(A[i][j] < value)
				A[i][j] = value;
		}
	}
}

double Matrix::trace() const
{
	double t = 0;
	for (int i = 0; i < std::min(m, n); ++i)
		t += A[i][i];
	return t;
}

Matrix Matrix::random(int m, int n)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);
	
	Matrix A{m, n};
	std::vector<std::vector<double>>& X = A.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			X[i][j] = dis(gen);
	}
	return A;
}

Matrix Matrix::identity(int m, int n)
{
	Matrix A{m, n};
	std::vector<std::vector<double>>& X = A.getArray();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			X[i][j] = (i == j ? 1.0 : 0.0);
	}
	return A;
}

void Matrix::print(int w, int d) const
{
	std::ios::fmtflags f(std::cout.flags());
	std::cout << std::setw(w) << std::setprecision(d) << '\n';
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			std::cout << std::to_string(A[i][j]) << ' ';
		std::cout << '\n';
	}
	std::cout << std::endl;
	std::cout.flags(f);
}

Matrix Matrix::readCSV(const char* file)
{
	std::ifstream ifs(file, std::ifstream::in);
	if (!ifs.is_open())
	{
		std::string ex = file;
		ex += " doesn't exist";
		throw std::runtime_error(ex);
	}
	std::istringstream iss;
	std::string token;
	
	std::vector<std::vector<double>> data;
	
	size_t row = 0;
	size_t col = 0;
	size_t cols = 0;
	for (std::string line; std::getline(ifs, line); )
	{
		col = 0;
		data.push_back(std::vector<double>());
		
		iss.clear();
		iss.str(line);
		
		while(std::getline(iss, token, ','))
		{
			data[row].push_back(std::stod(token));
			++col;
		}
		if (row == 0)
			cols = col;
		else if (cols != col)
			throw std::logic_error("Uneven amount of columns in a matrix");
		
		++row;
	}
	
	return Matrix(data);
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix)
{
	int d = 2;
	int w = 2;
	std::ios::fmtflags f(os.flags());
	os << std::setw(w) << std::setprecision(d) << '\n';
	for (int i = 0; i < matrix.m; ++i)
	{
		for (int j = 0; j < matrix.n; ++j)
			os << std::to_string(matrix.A[i][j]) << ' ';
		os << '\n';
	}
	os << '\n';
	os.flags(f);
	return os;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	return lhs.plus(rhs);
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	return lhs.minus(rhs);
}

Matrix& operator+=(Matrix& lhs, const Matrix& rhs)
{
	return lhs.plusEquals(rhs);
}
 
Matrix& operator-=(Matrix& lhs, const Matrix& rhs)
{
	return lhs.minusEquals(rhs);
}

void Matrix::checkMatrixDimensions(const Matrix& B) const
{
	if (B.m != m || B.n != n)
		throw IllegalArgumentException("Matrix dimensions must agree.");
}

Matrix Matrix::mean(int dim) const
{
	Matrix result;
	switch (dim) {
		case 1:
			result = Matrix{1, n};
			for (int currN = 0; currN < n; ++currN)
			{
				for (int currM = 0; currM < m; ++currM)
					result.A[0][currN] += A[currM][currN];
				result.A[0][currN] /= m;
			}
			return result;
		case 2:
			result = Matrix{m, 1};
			for (int currM = 0; currM < m; ++currM)
			{
				for (int currN = 0; currN < n; ++currN)
					result.A[currM][0] += A[currM][currN];
				result.A[currM][0] /= n;
			}
			return result;
		default:
			std::string ex = "dim must be either 1 or 2, and not: ";
			ex += dim;
			throw IllegalArgumentException(ex);
	}
}

Matrix Matrix::cov()
{
	Matrix transe = transpose();
	Matrix result{transe.m, transe.m};
	for(int currM = 0; currM < transe.m; ++currM)
	{
		for(int currN = currM; currN < transe.m; ++currN)
		{
			double covMN = cov(transe.A[currM], transe.A[currN]);
			result.A[currM][currN] = covMN;
			result.A[currN][currM] = covMN;
		}
	}
	return result;
}

double Matrix::cov(const std::vector<double>& vec1, const std::vector<double>& vec2)
{
	double result = 0;
	int dim = static_cast<int>(vec1.size());
	if(static_cast<int>(vec2.size()) != dim)
		throw IllegalArgumentException("vectors are not of same length");
	double meanVec1 = mean(vec1);
	double meanVec2 = mean(vec2);
	for (int i = 0; i < dim; ++i)
		result += (vec1[i] - meanVec1) * (vec2[i] - meanVec2);
	return result / std::max(1, dim - 1);
}

double Matrix::mean(const std::vector<double>& vec)
{
	double result = 0;
	for (size_t i = 0; i < vec.size(); ++i)
		result += vec[i];
	return result / vec.size();
}

double Matrix::sum() const
{
	double result = 0;
	for(auto& dArr : A)
		for(double d : dArr)
			result += d;
	return result;
}

Matrix Matrix::abs() const
{
	Matrix result{m, n}; // Don't use clone(), as the values are assigned in the loop.
	for(size_t i = 0; i < result.A.size(); ++i)
	{
		for(size_t j = 0; j < result.A[i].size(); ++j)
			result.A[i][j] = std::abs(A[i][j]);
	}
	return result;
}

void Matrix::writeAscii(const char* filename) const
{
	std::ofstream ofs(filename, std::ios::out);
	std::ios::fmtflags f(ofs.flags());
	ofs << std::scientific << std::setprecision(7);
	if (ofs.is_open())
	{
		for(int i = 0; i < m; ++i)
		{
			for(int j = 0; j < n; ++j)
				ofs << A[i][j];
			ofs << "\n";
		}
		ofs.flags(f);
	}
	else
	{
		std::string ex = "there was a problem with the file ";
		ex += filename;
		throw IllegalArgumentException(ex);
	}
}