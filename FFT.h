#ifndef FFT_H
#define FFT_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <vector>

#include "IllegalArgumentException.h"

#if !defined(M_PI)
	#define M_PI 3.14159265358979323846
#endif

static constexpr double twoPI = 2 * M_PI;

/** used to specify a forward Fourier transform */
static constexpr int FFT_FORWARD            = -1;
/** used to specify an inverse Fourier transform */
static constexpr int FFT_REVERSE            =  1;
/** used to specify a magnitude Fourier transform */
static constexpr int FFT_MAGNITUDE          =  2;
/** used to specify a magnitude phase Fourier transform */
static constexpr int FFT_MAGNITUDE_PHASE    =  3;
/** used to specify a normalized power Fourier transform */
static constexpr int FFT_NORMALIZED_POWER   =  4;
/** used to specify a power Fourier transform */
static constexpr int FFT_POWER              =  5;
/** used to specify a power phase Fourier transform */
static constexpr int FFT_POWER_PHASE        =  6;
/** used to specify a inline power phase Fourier transform */
static constexpr int FFT_INLINE_POWER_PHASE =  7;

/** used to specify a rectangular window function */
static constexpr int WND_NONE         = -1;
/** used to specify a rectangular window function */
static constexpr int WND_RECT         =  0;
/** used to specify a Hamming window function */
static constexpr int WND_HAMMING      =  1;
/** used to specify a 61-dB 3-sample Blackman-Harris window function */
static constexpr int WND_BH3          =  2;
/** used to specify a 74-dB 4-sample Blackman-Harris window function */
static constexpr int WND_BH4          =  3;
/** used to specify a minimum 3-sample Blackman-Harris window function */
static constexpr int WND_BH3MIN       =  4;
/** used to specify a minimum 4-sample Blackman-Harris window function */
static constexpr int WND_BH4MIN       =  5;
/** used to specify a Gaussian window function */
static constexpr int WND_GAUSS        =  6;
/** used to specify a Hanning window function */
static constexpr int WND_HANNING      =  7;
/** used to specify a Hanning window function */
static constexpr int WND_USER_DEFINED =  8;
/** used to specify a Hanning Z window function */
static constexpr int WND_HANNINGZ     =  9;

class FFT
{
public:
	FFT(int transformationType, int windowSize);
	FFT(int transformationType, int windowSize, int windowFunctionType);
	FFT(int transformationType, int windowSize, int windowFunctionType, int support);
	FFT(int transformationType, int windowSize, const std::vector<double>& windowFunction);
	
	void transform(std::vector<double>& re, std::vector<double>& im);
	/**
	  * This method allows to change the window function to one of the predefined
	  * window function types.
	  *
	  * @param windowFunctionType int the type of the window function
	  * @param support int
	  */
	void setWindowFunction(int windowFunctionType, int support);
	int getTransformationType() const;
	int getWindowFunctionType() const;
	
private:
	/** The FFT method. Calculation is inline, for complex data stored
	  * in 2 separate arrays. Length of input data must be a power of two.
	  * @param re        the real part of the complex input and output data
	  * @param im        the imaginary part of the complex input and output data
	  * @param direction the direction of the Fourier transform (FORWARD or
	  * REVERSE)
	  * @throws IllegalArgumentException if the length of the input data is
	  * not a power of 2
	  */
	void fft(std::vector<double>& re, std::vector<double>& im, int direction);
	/** Computes the power spectrum of a real sequence (in place).
	  * @param re the real input and output data; length must be a power of 2
	  */
	void powerFFT(std::vector<double>& re);
	/** Computes the magnitude spectrum of a real sequence (in place).
	  * @param re the real input and output data; length must be a power of 2
	  */
	void magnitudeFFT(std::vector<double>& re);
	/** Computes the power spectrum of a real sequence (in place).
	  *  @param re the real input and output data; length must be a power of 2
	  */
	void normalizedPowerFFT(std::vector<double>& re);
	/** Converts a real power sequence from to magnitude representation,
	  * by computing the square root of each value.
	  * @param re the real input (power) and output (magnitude) data; length
	  * must be a power of 2
	  */
	void toMagnitude(std::vector<double>& re);
	/** Computes a complex (or real if im[] == {0,...}) FFT and converts
	  * the results to polar coordinates (power and phase). Both arrays
	  * must be the same length, which is a power of 2.
	  * @param re the real part of the input data and the power of the output
	  * data
	  * @param im the imaginary part of the input data and the phase of the
	  * output data
	  */
	void powerPhaseFFT(std::vector<double>& re, std::vector<double>& im);
	/** Inline computation of the inverse FFT given spectral input data
	  * in polar coordinates (power and phase).
	  * Both arrays must be the same length, which is a power of 2.
	  * @param pow the power of the spectral input data (and real part of the
	  * output data)
	  * @param ph the phase of the spectral input data (and the imaginary part
	  * of the output data)
	  */
	void powerPhaseIFFT(std::vector<double>& pow, std::vector<double>& ph);
	/** Computes a complex (or real if im[] == {0,...}) FFT and converts
	  * the results to polar coordinates (magnitude and phase). Both arrays
	  * must be the same length, which is a power of 2.
	  * @param re the real part of the input data and the magnitude of the
	  * output data
	  * @param im the imaginary part of the input data and the phase of the
	  * output data
	  */
	void magnitudePhaseFFT(std::vector<double>& re, std::vector<double>& im);
	/** Fill an array with the values of a standard Hamming window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void hamming(int size);
	/** Fill an array with the values of a standard Hanning window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void hanning(int size);
	/** In MATLABTM, picking up the standard hanning window gives an incorrect periodicity,
	  * because the boundary samples are non-zero; in OCTAVETM, both boundary samples 
	  * are zero, which still gives an incorrect periodicity. This is why we use hanningz,
	  * a modified version of the hanning window available with the MATLABTM toolboxes: 
	  * function w = hanningz(n) w = .5*(1 - cos(2*pi*(0:n-1)'/(n)));
	  * 
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void hanningz(int size);
	/** Fill an array with the values of a minimum 4-sample Blackman-Harris
	  * window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void blackmanHarris4sMin(int size);
	/** Fill an array with the values of a 74-dB 4-sample Blackman-Harris
	  * window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void blackmanHarris4s(int size);
	/** Fill an array with the values of a minimum 3-sample Blackman-Harris
	  * window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void blackmanHarris3sMin(int size);
	/** Fill an array with the values of a 61-dB 3-sample Blackman-Harris
	  * window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void blackmanHarris3s(int size);
	/** Fill an array with the values of a Gaussian window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void gauss(int size);
	/** Fill an array with the values of a rectangular window function
	  * @param data the array to be filled
	  * @param size the number of non zero values; if the array is larger than
	  * this, it is zero-padded symmetrically at both ends
	  */
	void rectangle(int size);
	/** Applies a window function to an array of data, storing the result in
	  * the data array.
	  * Performs a dot product of the data and window arrays.
	  * @param data   the array of input data, also used for output
	  * @param window the values of the window function to be applied to data
	  */
	void applyWindowFunction(std::vector<double>& data);
	void calculateWindowFunctionSum();
	
	std::vector<double> windowFunction;
	double windowFunctionSum;
	int windowFunctionType;
	int transformationType;
	int windowSize;
};

#endif /* FFT_H */