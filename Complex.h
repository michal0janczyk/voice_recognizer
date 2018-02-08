#ifndef COMPLEX_H
#define COMPLEX_H

#include <ostream>
#include <string>

class Complex // or just use <complex>...
{
public:
	/** Create a new object with the given real and imaginary parts */
	Complex(double real, double imag);
	
	/** Return a string representation of the invoking Complex object */
	std::string toString() const;
	/** Return abs/modulus/magnitude */
	double abs() const;
	/** Return angle/phase/argument */
	double phase() const;
	/** Return a new Complex object whose value is (this + b) */
	Complex plus(const Complex& b) const;
	/** Return a new Complex object whose value is (this - b) */
	Complex minus(const Complex& b) const;
	/** Return a new Complex object whose value is (this * b) */
	Complex times(const Complex& b) const;
	/** 
	  * Scalar multiplication
      * Return a new object whose value is (this * alpha)
	  */
	Complex times(double alpha) const;
	/** Return a new Complex object whose value is the conjugate of this */
	Complex conjugate() const;
	/** Return a new Complex object whose value is the reciprocal of this */
	Complex reciprocal() const;
	/** Return the real part */
	double re() const;
	/** Return the imaginary part */
	double im() const;
	/** Return a / b */
	Complex divides(const Complex& b) const;
	/** Return a new Complex object whose value is the complex exponential of this */
	Complex exp() const;
	/** Return a new Complex object whose value is the complex sine of this */
	Complex sin() const;
	/** Return a new Complex object whose value is the complex cosine of this */
	Complex cos() const;
	/** Return a new Complex object whose value is the complex tangent of this */
	Complex tan() const;
	/** A static version of plus */
	static Complex plus(const Complex& a, const Complex& b);
	
private:
	friend std::ostream& operator<<(std::ostream& os, const Complex& complex);
	friend Complex operator+(const Complex& lhs, const Complex& rhs);
	friend Complex operator-(const Complex& lhs, const Complex& rhs);
	friend Complex operator*(const Complex& lhs, const Complex& rhs);
	friend Complex operator*(const Complex& lhs, double value);
	friend Complex operator/(const Complex& lhs, const Complex& rhs);
	
	const double real; /** the imaginary part */
	const double imag; /** the imaginary part */
};

#endif /* COMPLEX_H */