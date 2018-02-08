#include "Complex.h"

#include <cmath>

Complex::Complex(double real, double imag) : real(real), imag(imag)
{
	// Empty
}

std::string Complex::toString() const
{
	std::string result;
	if (imag == 0)
		return result + std::to_string(real);
	if (real == 0)
		return result + std::to_string(imag) + 'i';
	if (imag < 0)
		return result + std::to_string(real) + " - " + std::to_string(-imag) + 'i';
	return result + std::to_string(real) + " + " + std::to_string(imag) + 'i';
}

double Complex::abs() const
{
	return std::hypot(real, imag); // std::sqrt(real*real + imag*imag)
}

double Complex::phase() const
{
	return std::atan2(imag, real); // Between -pi and pi
}

Complex Complex::plus(const Complex& b) const
{
	return Complex(real + b.real, imag + b.imag);
}

Complex Complex::minus(const Complex& b) const
{
	return Complex(real - b.real, imag - b.imag);
}

Complex Complex::times(const Complex& b) const
{
	return Complex(real * b.real - imag * b.imag, real * b.imag + imag * b.real);
}

Complex Complex::times(double alpha) const
{
	return Complex(alpha * real, alpha * imag);
}

Complex Complex::conjugate() const
{
	return Complex(real, -imag);
}

Complex Complex::reciprocal() const
{
	double scale = real * real + imag * imag;
	return Complex(real / scale, -imag / scale);
}

double Complex::re() const
{
	return real;
}

double Complex::im() const
{
	return imag;
}

Complex Complex::divides(const Complex& b) const
{
	return *this * (b.reciprocal());
}

Complex Complex::exp() const
{
	return Complex(std::exp(real) * std::cos(imag), std::exp(real) * std::sin(imag));
}

Complex Complex::sin() const
{
	return Complex(std::sin(real) * std::cosh(imag), std::cos(real) * std::sinh(imag));
}

Complex Complex::cos() const
{
	return Complex(std::cos(real) * std::cosh(imag), -std::sin(real) * std::sinh(imag));
}

Complex Complex::tan() const
{
	return sin() / cos();
}

std::ostream& operator<<(std::ostream& os, const Complex& complex)
{
	if (complex.imag == 0)
		return os << complex.real;
	if (complex.real == 0)
		return os << complex.imag << 'i';
	if (complex.imag < 0)
		return os << complex.real << " - " << -complex.imag << 'i';
	return os << complex.real << " + " << complex.imag << 'i';
}

Complex Complex::plus(const Complex& a, const Complex& b)
{
	return Complex(a.real + b.real, a.imag + b.imag);
}

Complex operator+(const Complex& lhs, const Complex& rhs)
{
	return Complex(lhs.real + rhs.real, lhs.imag + rhs.imag);
}

Complex operator-(const Complex& lhs, const Complex& rhs)
{
	return Complex(lhs.real - rhs.real, lhs.imag - rhs.imag);
}

Complex operator*(const Complex& lhs, const Complex& rhs)
{
	return Complex(lhs.real * rhs.real - lhs.imag * rhs.imag, lhs.real * rhs.imag + lhs.imag * rhs.real);
}

Complex operator/(const Complex& lhs, const Complex& rhs)
{
	return lhs * (rhs.reciprocal());
}