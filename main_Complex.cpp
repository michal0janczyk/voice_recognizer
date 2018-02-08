#include <iostream>

#include "Complex.h"

int main()
{
	Complex a{5.0, 6.0};
	Complex b{-3.0, 4.0};

	std::cout << "a            = " << a             << std::endl;
	std::cout << "b            = " << b             << std::endl;
	std::cout << "Re(a)        = " << a.re()        << std::endl;
	std::cout << "Im(a)        = " << a.im()        << std::endl;
	std::cout << "b + a        = " << (b + a)       << std::endl;
	std::cout << "a - b        = " << (a - b)       << std::endl;
	std::cout << "a * b        = " << (a * b)       << std::endl;
	std::cout << "b * a        = " << (b * a)       << std::endl;
	std::cout << "a / b        = " << (a / b)       << std::endl;
	std::cout << "(a / b) * b  = " << ((a / b) * b) << std::endl;
	std::cout << "conj(a)      = " << a.conjugate() << std::endl;
	std::cout << "|a|          = " << a.abs()       << std::endl;
	std::cout << "tan(a)       = " << a.tan()       << std::endl;
	
	return 0;
}