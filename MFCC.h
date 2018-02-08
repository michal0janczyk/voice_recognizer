#ifndef MFCC_H
#define MFCC_H

#include <vector>

#include "IllegalArgumentException.h"
#include "Matrix.h"

class FFT;

class MFCC
{
public:
	MFCC(float sampleRate, int windowSize, int numberCoefficients, bool useFirstCoefficient);
	MFCC(float sampleRate, int windowSize, int numberCoefficients, bool useFirstCoefficient, double minFreq, double maxFreq, int numberFilters);
	~MFCC();
	
	/**
	  * Performs the transformation of the input data to MFCCs.
	  * This is done by splitting the given data into windows and processing
	  * each of these windows with processWindow().
	  *
	  * @param input double[] input data is an array of samples, must be a multiple
	  *                       of the hop size, must not be a null value
	  * @return double[][] an array of arrays contains a double array of Sone value
	  *                    for each window
	  * @throws IOException if there are any problems regarding the inputstream
	  * @throws IllegalArgumentException raised if method contract is violated
	  */
	std::vector<std::vector<double>> process(const std::vector<double>& input);
	int getWindowSize() const;
	/**
	  * Transforms one window of MFCCs. The following steps are
	  * performed: <br>
	  * <br>
	  * (1) normalized power fft with hanning window function<br>
	  * (2) convert to Mel scale by applying a mel filter bank<br>
	  * (3) Conversion to db<br>
	  * (4) finally a DCT is performed to get the mfcc<br>
	  *<br>
	  * This process is mathematical identical with the process described in [1].
	  *
	  * @param window double[] data to be converted, must contain enough data for
	  *                        one window
	  * @param start int start index of the window data
	  * @return double[] the window representation in Sone
	  * @throws IllegalArgumentException raised if method contract is violated
	  */
	std::vector<double> processWindow(const std::vector<double>& window, int start);
	
private:
	MFCC(const MFCC& mfcc) = delete;
	MFCC(MFCC&& mfcc) = delete;
	
	void operator=(const MFCC& mfcc) = delete;
	void operator=(MFCC&& mfcc) = delete;
	
	/**
	  * Returns the boundaries (start, center, end) of a given number of triangular
	  * mel filters at linear scale. Mel-filters are triangular filters on the
	  * linear scale with an integral (area) of 1. However they are placed
	  * equidistantly on the mel scale, which is non-linear rather logarithmic.
	  * The minimum linear frequency and the maximum linear frequency define the
	  * mel-scaled interval to equidistantly place the filters.
	  * Since mel-filters overlap, an array is used to efficiently store the
	  * boundaries of a filter. For example you can get the boundaries of the k-th
	  * filter by accessing the returned array as follows:
	  *
	  * leftBoundary = boundaries[k-1];
	  * center = boundaries[k];
	  * rightBoundary = boundaries[k+1];
	  *
	  * @param minFreq double frequency used for the left boundary of the first
	  *                       filter
	  * @param maxFreq double frequency used for the right boundary of the last
	  *                       filter
	  * @param numberFilters int number of filters to place within the interval
	  *                          [minFreq, maxFreq]
	  * @return double[] array holding the boundaries
	  */
	std::vector<double> getMelFilterBankBoundaries(double minFreq, double maxFreq, int numberFilters);
	/**
	  * This method creates a matrix containing <code>numberFilters</code>
	  * mel-filters. Each filter is represented by one row of this matrix. Thus all
	  * the filters can be applied at once by a simple matrix multiplication.
	  *
	  * @return Matrix a matrix containing the filter banks
	  */
	Matrix getMelFilterBanks();
	/**
	  * Returns the filter weight of a given mel filter at a given frequency.
	  * Mel-filters are triangular filters on the linear scale with an integral
	  * (area) of 1. However they are placed equidistantly on the mel scale, which
	  * is non-linear rather logarithmic.
	  * Consequently there are lots of high, thin filters at start of the linear
	  * scale and rather few and flat filters at the end of the linear scale.
	  * Since the start-, center- and end-points of the triangular mel-filters on
	  * the linear scale are known, the weights are computed using linear
	  * interpolation.
	  *
	  * @param filterBank int the number of the mel-filter, used to extract the
	  *                       boundaries of the filter from the array
	  * @param freq double    the frequency, at which the filter weight should be
	  *                       returned
	  * @param boundaries double[] an array containing all the boundaries
	  * @return double the filter weight
	  */
	double getMelFilterWeight(int filterBank, double freq, const std::vector<double>& boundaries) const;
	double linToMelFreq(double inputFreq) const;
	double melToLinFreq(double inputFreq) const;
	/**
	  * Generates the DCT matrix for the known number of filters (input vector) and
	  * for the known number of used coefficients (output vector). Therefore the
	  * DCT matrix has the dimensions (numberCoefficients x numberFilters).
	  * If useFirstCoefficient is set to false the matrix dimensions are
	  * (numberCoefficients-1 x numberFilters). This matrix is a submatrix of the
	  * full matrix. Only the first row is missing.
	  *
	  * @return Matrix the appropriate DCT matrix
	  */
	Matrix getDCTMatrix();
	
	int windowSize;
	int hopSize;
	float sampleRate;
	double baseFreq;
	
	double minFreq;
	double maxFreq;
	int numberFilters;
	
	int numberCoefficients;
	bool useFirstCoefficient;
	
	std::vector<double> buffer;
	Matrix dctMatrix;
	Matrix melFilterBanks;
	FFT* normalizedPowerFFT;
	double scale;
};

#endif /* MFCC_H */