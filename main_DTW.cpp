#include <iostream>

#include "DTW.h"

int main()
{
	std::vector<double> n2 = { 1.5, 3.9, 4.1, 3.3 };
	std::vector<double> n1 = { 2.1, 2.45, 3.673, 4.32, 2.05, 1.93, 5.67, 6.01 };
	DTW dtw{n1, n2};
	std::cout << dtw;
	
	return 0;
}