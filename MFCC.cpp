#include "MFCC.h"

#define _USE_MATH_DEFINES

#include <cmath>
#include <limits>

#include "FFT.h"

#if !defined(M_PI)
	#define M_PI 3.14159265358979323846
#endif

MFCC::MFCC(float sampleRate, int windowSize, int numberCoefficients, bool useFirstCoefficient) : MFCC(sampleRate, windowSize, numberCoefficients, useFirstCoefficient, 20.0, 16000.0, 40)
{
	// Empty
}

MFCC::MFCC(float sampleRate, int windowSize, int numberCoefficients, bool useFirstCoefficient, double minFreq, double maxFreq, int numberFilters)
{
	if (windowSize < 32)
		throw IllegalArgumentException("window size must be at least 32");
	else
	{
		int i = 32;
		while (i < windowSize && i < std::numeric_limits<int>::max())
			i = i << 1;
		
		if (i != windowSize)
			throw IllegalArgumentException("window size must be 2^n");
	}
	
	this->sampleRate = sampleRate;
	this->windowSize = windowSize;
	this->hopSize = windowSize / 2; // 50% Overlap
	this->baseFreq = sampleRate / windowSize;
	
	this->numberCoefficients = numberCoefficients;
	this->useFirstCoefficient = useFirstCoefficient;
	
	this->minFreq = minFreq;
	this->maxFreq = maxFreq;
	this->numberFilters = numberFilters;
	
	buffer.resize(windowSize, 0.0);
	
	melFilterBanks = getMelFilterBanks();
	dctMatrix = getDCTMatrix();
	
	normalizedPowerFFT = new FFT(FFT_NORMALIZED_POWER, windowSize, WND_HANNING);
	
	// Compute rescale factor to rescale and normalize at once (default is 96dB = 2^16)
	scale = (std::pow(10, 96 / 20));   
}

MFCC::~MFCC()
{
	delete normalizedPowerFFT;
}

std::vector<double> MFCC::getMelFilterBankBoundaries(double minFreq, double maxFreq, int numberFilters)
{
	// Create return array
	std::vector<double> centers;
	centers.resize(numberFilters + 2, 0.0);
	double maxFreqMel, minFreqMel, deltaFreqMel, nextCenterMel;
	
	// Compute mel min./max. frequency
	maxFreqMel = linToMelFreq(maxFreq);
	minFreqMel = linToMelFreq(minFreq);
	deltaFreqMel = (maxFreqMel - minFreqMel) / (numberFilters + 1);
	
	// Create (numberFilters + 2) equidistant points for the triangles
	nextCenterMel = minFreqMel;
	for (size_t i = 0; i < centers.size(); ++i)
	{
		// Transform the points back to linear scale
		centers[i] = melToLinFreq(nextCenterMel);
		nextCenterMel += deltaFreqMel;
	}
	
	// Ajust boundaries to exactly fit the given min./max. frequency
	centers[0] = minFreq;
	centers[numberFilters + 1] = maxFreq;
	
	return centers;
}

Matrix MFCC::getMelFilterBanks()
{
	// Get boundaries of the different filters
	std::vector<double> boundaries = getMelFilterBankBoundaries(minFreq, maxFreq, numberFilters);
	
	// Ignore filters outside of spectrum
	for (int i = 1; i < static_cast<int>(boundaries.size() - 1); ++i)
	{
		if (boundaries[i] > sampleRate / 2)
		{
			numberFilters = i - 1;
			break;
		}
	}
	
	// Create the filter bank matrix
	std::vector<std::vector<double>> matrix;
	matrix.resize(numberFilters);
	
	// Fill each row of the filter bank matrix with one triangular mel filter
	for (int i = 1; i <= numberFilters; ++i)
	{
		std::vector<double> filter;
		filter.resize((windowSize / 2) + 1);
		
		// For each frequency of the fft
		for (size_t j = 0; j < filter.size(); ++j)
		{
			// Compute the filter weight of the current triangular mel filter
			double freq = baseFreq * j;
			filter[j] = getMelFilterWeight(i, freq, boundaries);
		}
		
		// Add the computed mel filter to the filter bank
		matrix[i - 1] = filter;
	}
	
	// Return the filter bank
	return Matrix(matrix, numberFilters, (windowSize / 2) + 1);
}

double MFCC::getMelFilterWeight(int filterBank, double freq, const std::vector<double>& boundaries) const
{
	// For most frequencies the filter weight is 0
	double result = 0;
	
	// Compute start- , center- and endpoint as well as the height of the filter
	double start = boundaries[filterBank - 1];
	double center = boundaries[filterBank];
	double end = boundaries[filterBank + 1];
	double height = 2.0 / (end - start);
	
	// Is the frequency within the triangular part of the filter
	if (freq >= start && freq <= end)
	{
		// Depending on frequency position within the triangle
		if (freq < center)
		{
			// ...use a ascending linear function
			result = (freq - start) * (height / (center - start));
		}
		else
		{
			// ..use a descending linear function
			result = height + ((freq - center) * (-height / (end - center)));
		}
	}
	
	return result;
}

double MFCC::linToMelFreq(double inputFreq) const
{
	return (2595.0 * (std::log(1.0 + inputFreq / 700.0) / std::log(10.0)));
}

double MFCC::melToLinFreq(double inputFreq) const
{
	return (700.0 * (std::pow(10.0, (inputFreq / 2595.0)) - 1.0));
}

Matrix MFCC::getDCTMatrix()
{
	// Compute constants
	double k = M_PI / numberFilters;
	double w1 = 1.0 / (std::sqrt(numberFilters)); // 1.0 / (std::sqrt(numberFilters / 2));
	double w2 = std::sqrt(2.0 / numberFilters);   // std::sqrt(2.0 / numberFilters) * (std::sqrt(2.0) / 2.0);
	
	// Create new matrix
	Matrix matrix{numberCoefficients, numberFilters};
	
	// Generate dct matrix
	for (int i = 0; i < numberCoefficients; ++i)
	{
		for (int j = 0; j < numberFilters; ++j)
		{
			if (i == 0)
				matrix.set(i, j, w1 * std::cos(k * i * (j + 0.5)));
			else
				matrix.set(i, j, w2 * std::cos(k * i * (j + 0.5)));
		}
	}
	
	// Adjust index if we are using first coefficient
	if (!useFirstCoefficient)
		matrix = matrix.getMatrix(1, numberCoefficients - 1, 0, numberFilters - 1);
	
	return matrix;
}

std::vector<std::vector<double>> MFCC::process(const std::vector<double>& input)
{
	// Check for null
	if (!input.size())
		throw IllegalArgumentException("input data must not be a null value");
	
	// Check for correct array length
	if ((input.size() % hopSize) != 0)
		throw IllegalArgumentException("Input data must be multiple of hop size (windowSize/2).");
	
	// Create return array with appropriate size
	std::vector<std::vector<double>> mfcc;
	mfcc.resize((input.size() / hopSize) - 1, std::vector<double>(numberCoefficients));
	
	// Process each window of this audio segment
	for (int i = 0, pos = 0; pos < static_cast<int>(input.size() - hopSize); ++i, pos += hopSize)
		mfcc[i] = processWindow(input, pos);
	
	return mfcc;
}

int MFCC::getWindowSize() const
{
	return windowSize;
}

std::vector<double> MFCC::processWindow(const std::vector<double>& window, int start)
{
	// Number of unique coefficients, and the rest are symmetrically redundant
	int fftSize = (windowSize / 2) + 1;
	
	// Check start
	if (start < 0)
		throw IllegalArgumentException("start must be a positve value");
	
	// Check window size
	if (!window.size() || static_cast<int>(window.size()) - start < windowSize)
		throw IllegalArgumentException("the given data array must not be a null value and must contain data for one window");
	
	// Just copy to buffer and rescaled the input samples according to the original matlab implementation to 96dB
	for (int j = 0; j < windowSize; ++j)
		buffer[j] = window[j + start] * scale;
	
	// Perform power fft
	std::vector<double> dummyIm; // Won't be used
	normalizedPowerFFT->transform(buffer, dummyIm);
	
	// Use all coefficient up to the nyquist frequency (ceil((fftSize+1)/2))
	Matrix x{buffer, windowSize};
	x = x.getMatrix(0, fftSize - 1, 0, 0); // fftSize - 1 is the index of the nyquist frequency
	
	// Apply mel filter banks
	x = melFilterBanks.times(x);
	
	// To db
	double log10 = 10 * (1 / std::log(10)); // log for base 10 and scale by factor 10
	x.thrunkAtLowerBoundary(1);
	x.logEquals();
	x.timesEquals(log10);
	
	// Compute DCT
	x = dctMatrix.times(x);
	
	return x.getColumnPackedCopy();
}