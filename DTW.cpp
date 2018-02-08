#include "DTW.h"

#include <algorithm>

DTW::DTW(const std::vector<double>& sample, const std::vector<double>& templete)
{
	seq1 = sample;
	seq2 = templete;
	
	n = static_cast<int>(seq1.size());
	m = static_cast<int>(seq2.size());
	K = 1;
	
	warpingPath.resize(n + m, std::vector<int>(2, 0)); // max(n, m) <= K < n + m
	warpingDistance = 0.0;
	
	compute();
}

void DTW::compute()
{
	double accumulatedDistance = 0.0;
	
	std::vector<std::vector<double>> d; // local distances
	d.resize(n, std::vector<double>(m, 0.0));
	std::vector<std::vector<double>> D; // global distances
	D.resize(n, std::vector<double>(m, 0.0));
	
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
			d[i][j] = distanceBetween(seq1[i], seq2[j]);
	}
	
	D[0][0] = d[0][0];
	
	for (int i = 1; i < n; ++i)
		D[i][0] = d[i][0] + D[i - 1][0];
	
	for (int j = 1; j < m; ++j)
		D[0][j] = d[0][j] + D[0][j - 1];
	
	for (int i = 1; i < n; ++i)
	{
		for (int j = 1; j < m; ++j)
		{
			accumulatedDistance = std::min(
				std::min(D[i - 1][j], D[i - 1][j - 1]),
				D[i][j - 1]
			);
			accumulatedDistance += d[i][j];
			D[i][j] = accumulatedDistance;
		}
	}
	accumulatedDistance = D[n - 1][m - 1];
	
	int i = n - 1;
	int j = m - 1;
	int minIndex = 1;
	
	warpingPath[K - 1][0] = i;
	warpingPath[K - 1][1] = j;
	
	while ((i + j) != 0)
	{
		if (i == 0)
			j -= 1;
		else if (j == 0)
			i -= 1;
		else // i != 0 && j != 0
		{
			std::vector<double> array = {
				D[i - 1][j],
				D[i][j - 1],
				D[i - 1][j - 1]
			};
			minIndex = getIndexOfMinimum(array);
			
			if (minIndex == 0)
				i -= 1;
			else if (minIndex == 1)
				j -= 1;
			else if (minIndex == 2)
			{
				i -= 1;
				j -= 1;
			}
		}
		K++;
		warpingPath[K - 1][0] = i;
		warpingPath[K - 1][1] = j;
	}
	warpingDistance = accumulatedDistance / K;
	
	reversePath(warpingPath);
}

void DTW::reversePath(const std::vector<std::vector<int>>& path)
{
	std::vector<std::vector<int>> newPath;
	newPath.resize(K, std::vector<int>(2, 0));
	for (int i = 0; i < K; ++i)
	{
		for (int j = 0; j < 2; ++j)
			newPath[i][j] = path[K - i - 1][j];
	}
	warpingPath = newPath;
}

double DTW::getDistance() const
{
	return warpingDistance;
}

double DTW::distanceBetween(double p1, double p2) const
{
	return (p1 - p2) * (p1 - p2);
}

int DTW::getIndexOfMinimum(const std::vector<double>& array) const
{
	int index = 0;
	double val = array[0];
	
	for (size_t i = 1; i < array.size(); ++i)
	{
		if (array[i] < val)
		{
			val = array[i];
			index = static_cast<int>(i);
		}
	}
	return index;
}

std::string DTW::toString() const
{
	std::string retVal = "Warping Distance: " + std::to_string(warpingDistance) + "\n";
	retVal += "Warping Path: {";
	for (int i = 0; i < K; ++i)
	{
		retVal += "(" + std::to_string(warpingPath[i][0]) + ", " + std::to_string(warpingPath[i][1]) + ")";
		retVal += (i == K - 1) ? "}" : ", ";
	}
	return retVal;
}

std::ostream& operator<<(std::ostream& os, const DTW& dtw)
{
	std::string retVal = "Warping Distance: " + std::to_string(dtw.warpingDistance) + "\n";
	retVal += "Warping Path: {";
	for (int i = 0; i < dtw.K; ++i)
	{
		retVal += "(" + std::to_string(dtw.warpingPath[i][0]) + ", " + std::to_string(dtw.warpingPath[i][1]) + ")";
		retVal += (i == dtw.K - 1) ? "}" : ", ";
	}
	return os << retVal;
}