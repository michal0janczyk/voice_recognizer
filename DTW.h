#ifndef DTW_H
#define DTW_H

#include <ostream>
#include <string>
#include <vector>

class DTW
{
public:
	DTW(const std::vector<double>& sample, const std::vector<double>& templete);
	
	void compute();
	/**
	  * Changes the order of the warping path (increasing order)
	  * 
	  * @param path
	  *            the warping path in reverse order
	  */
	void reversePath(const std::vector<std::vector<int>>& path);
	/**
	 * Returns the warping distance
	 * 
	 * @return
	 */
	double getDistance() const;
	/**
	 * Computes a distance between two points
	 * 
	 * @param p1
	 *            the point 1
	 * @param p2
	 *            the point 2
	 * @return the distance between two points
	 */
	double distanceBetween(double p1, double p2) const;
	/**
	 * Finds the index of the minimum element from the given array
	 * 
	 * @param array
	 *            the array containing numeric values
	 * @return the min value among elements
	 */
	int getIndexOfMinimum(const std::vector<double>& array) const;
	/**
	 * Returns a string that displays the warping distance and path
	 */
	std::string toString() const;

private:
	friend std::ostream& operator<<(std::ostream& os, const DTW& dtw);
	
	std::vector<double> seq1;
	std::vector<double> seq2;
	std::vector<std::vector<int>> warpingPath;
	
	int n;
	int m;
	int K;
	
	double warpingDistance;
};

#endif /* DTW_H */