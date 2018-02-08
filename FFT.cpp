#include "FFT.h"

FFT::FFT(int transformationType, int windowSize) : FFT(transformationType, windowSize, WND_NONE, windowSize)
{
	// Empty
}

FFT::FFT(int transformationType, int windowSize, int windowFunctionType) : FFT(transformationType, windowSize, windowFunctionType, windowSize)
{
	// Empty
}

FFT::FFT(int transformationType, int windowSize, int windowFunctionType, int support)
{
	// Check and set fft type
	this->transformationType = transformationType;
	if (transformationType < -1 || transformationType > 7)
	{
		transformationType = FFT_FORWARD;
		throw IllegalArgumentException("unknown fft type");
	}
	
	// Check and set windowSize
	this->windowSize = windowSize;
	if (windowSize != (1 << (static_cast<int>(std::rint(std::log(windowSize) / std::log(2))))))
		throw IllegalArgumentException("fft data must be power of 2");
	
	// Create window function buffer and set window function
	this->windowFunction.resize(windowSize, 0.0);
    setWindowFunction(windowFunctionType, support);
}

FFT::FFT(int transformationType, int windowSize, const std::vector<double>& windowFunction) : FFT(transformationType, windowSize, WND_NONE, windowSize)
{
	if (static_cast<int>(windowFunction.size()) != windowSize)
		throw IllegalArgumentException("size of window function match window size");
	else
	{
		this->windowFunction = windowFunction;
		this->windowFunctionType = WND_USER_DEFINED;
		calculateWindowFunctionSum();
	}
}

void FFT::transform(std::vector<double>& re, std::vector<double>& im)
{
	// Check for correct size of the real part data array
	if (static_cast<int>(re.size()) < windowSize)
		throw IllegalArgumentException("data array smaller than fft window size");
	
	// Apply the window function to the real part
	applyWindowFunction(re);
	
	// Perform the transformation
	switch(transformationType)
	{
		case FFT_FORWARD:
			// Check for correct size of the imaginary part data array
			if (static_cast<int>(im.size()) < windowSize)
				throw IllegalArgumentException("data array smaller than fft window size");
			else
				fft(re, im, FFT_FORWARD);
			break;
		case FFT_INLINE_POWER_PHASE:
			if(static_cast<int>(im.size()) < windowSize)
				throw IllegalArgumentException("data array smaller than fft window size");
			else
				powerPhaseIFFT(re, im);
			break;
		case FFT_MAGNITUDE:
			magnitudeFFT(re);
			break;
		case FFT_MAGNITUDE_PHASE:
			if(static_cast<int>(im.size()) < windowSize)
				throw IllegalArgumentException("data array smaller than fft window size");
			else
				magnitudePhaseFFT(re, im);
			break;
		case FFT_NORMALIZED_POWER:
			normalizedPowerFFT(re);
			break;
		case FFT_POWER:
			powerFFT(re);
			break;
		case FFT_POWER_PHASE:
			if(static_cast<int>(im.size()) < windowSize)
				throw IllegalArgumentException("data array smaller than fft window size");
			else
				powerPhaseFFT(re, im);
			break;
		case FFT_REVERSE:
			if(static_cast<int>(im.size()) < windowSize)
				throw IllegalArgumentException("data array smaller than fft window size");
			else
				fft(re, im, FFT_REVERSE);
			break;
	}
}

void FFT::fft(std::vector<double>& re, std::vector<double>& im, int direction)
{
	int n = static_cast<int>(re.size());
	int bits = static_cast<int>(std::rint(std::log(n) / std::log(2)));
	
	if (n != (1 << bits))
		throw IllegalArgumentException("fft data must be power of 2");
	
	int localN;
	int j = 0;
	for (int i = 0; i < n - 1; ++i)
	{
		if (i < j)
		{
			double temp = re[j];
			re[j] = re[i];
			re[i] = temp;
			temp = im[j];
			im[j] = im[i];
			im[i] = temp;
		}
		
		int k = n / 2;
		
		while ((k >= 1) &&  (k - 1 < j))
		{
			j = j - k;
			k = k / 2;
		}
		
		j = j + k;
	}
	
	for (int m = 1; m <= bits; ++m)
	{
		localN = 1 << m;
		double Wjk_r = 1;
		double Wjk_i = 0;
		double theta = twoPI / localN;
		double Wj_r = std::cos(theta);
		double Wj_i = direction * std::sin(theta);
		int nby2 = localN / 2;
		for (j = 0; j < nby2; ++j)
		{
			for (int k = j; k < n; k += localN)
			{
				int id = k + nby2;
				double tempr = Wjk_r * re[id] - Wjk_i * im[id];
				double tempi = Wjk_r * im[id] + Wjk_i * re[id];
				re[id] = re[k] - tempr;
				im[id] = im[k] - tempi;
				re[k] += tempr;
				im[k] += tempi;
			}
			double wtemp = Wjk_r;
			Wjk_r = Wj_r * Wjk_r  - Wj_i * Wjk_i;
			Wjk_i = Wj_r * Wjk_i  + Wj_i * wtemp;
		}
	}
}

void FFT::powerFFT(std::vector<double>& re)
{
	std::vector<double> im;
	im.resize(re.size());
	
	fft(re, im, FFT_FORWARD);
	
	for (size_t i = 0; i < re.size(); ++i)
		re[i] = re[i] * re[i] + im[i] * im[i];
}

void FFT::magnitudeFFT(std::vector<double>& re)
{
	std::vector<double> im;
	im.resize(re.size());
	
	fft(re, im, FFT_FORWARD);
	
	for (size_t i = 0; i < re.size(); ++i)
		re[i] = std::sqrt(re[i] * re[i] + im[i] * im[i]);
}

void FFT::normalizedPowerFFT(std::vector<double>& re)
{
	std::vector<double> im;
	im.resize(re.size());
	double r, i;
	
	fft(re, im, FFT_FORWARD);
	
	for (size_t j = 0; j < re.size(); ++j)
	{
		r = re[j] / windowFunctionSum * 2;
		i = im[j] / windowFunctionSum * 2;
		re[j] = r * r + i * i;
	}
}

void FFT::toMagnitude(std::vector<double>& re)
{
	for (size_t i = 0; i < re.size(); ++i)
		re[i] = std::sqrt(re[i]);
}

void FFT::powerPhaseFFT(std::vector<double>& re, std::vector<double>& im)
{
	fft(re, im, FFT_FORWARD);
	
	for (size_t i = 0; i < re.size(); ++i)
	{
		double pow = re[i] * re[i] + im[i] * im[i];
		im[i] = std::atan2(im[i], re[i]);
		re[i] = pow;
	}
}

void FFT::powerPhaseIFFT(std::vector<double>& pow, std::vector<double>& ph)
{
	toMagnitude(pow);
	
	for (size_t i = 0; i < pow.size(); ++i)
	{
		double re = pow[i] * std::cos(ph[i]);
		ph[i] = pow[i] * std::sin(ph[i]);
		pow[i] = re;
	}
	
	fft(pow, ph, FFT_REVERSE);
}

void FFT::magnitudePhaseFFT(std::vector<double>& re, std::vector<double>& im)
{
	powerPhaseFFT(re, im);
	toMagnitude(re);
}

void FFT::hamming(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double scale = 1.0 / static_cast<double>(size / 0.54);
	double factor = twoPI / static_cast<double>(size);
	
	for (int i = 0; start < stop; ++start, ++i)
		windowFunction[i] = scale * (25.0 / 46.0 - 21.0 / 46.0 * std::cos(factor * i));
}

void FFT::hanning(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double factor = twoPI / (size - 1.0);
	
	for (int i = 0; start < stop; ++start, ++i)
		windowFunction[i] = 0.5 * (1 - std::cos(factor * i));
}

void FFT::hanningz(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	
	for (int i = 0; start < stop; ++start, ++i)
		windowFunction[i] = 0.5 * (1 - std::cos((twoPI * i) / size));
}

void FFT::blackmanHarris4sMin(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double scale = 1.0 / static_cast<double>(size) / 0.36;
	
	for (int i = 0; start < stop; ++start, ++i)
	windowFunction[i] = scale * (0.35875 -
								 0.48829 * std::cos(twoPI * i / size) +
								 0.14128 * std::cos(2 * twoPI * i / size) -
								 0.01168 * std::cos(3 * twoPI * i / size));
}

void FFT::blackmanHarris4s(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double scale = 1.0 / static_cast<double>(size) / 0.4;
	
	for (int i = 0; start < stop; ++start, ++i)
	windowFunction[i] = scale * (0.40217 -
								 0.49703 * std::cos(twoPI * i / size) +
								 0.09392 * std::cos(2 * twoPI * i / size) -
								 0.00183 * std::cos(3 * twoPI * i / size));
}

void FFT::blackmanHarris3sMin(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double scale = 1.0 / static_cast<double>(size) / 0.42;
	
	for (int i = 0; start < stop; ++start, ++i)
	windowFunction[i] = scale * (0.42323 -
								 0.49755 * std::cos(twoPI * i / size) +
								 0.07922 * std::cos(2 * twoPI * i / size));
}

void FFT::blackmanHarris3s(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double scale = 1.0 / static_cast<double>(size) / 0.45;
	
	for (int i = 0; start < stop; ++start, ++i)
	windowFunction[i] = scale * (0.44959 -
								 0.49364 * std::cos(twoPI * i / size) +
								 0.05677 * std::cos(2 * twoPI * i / size));
}

void FFT::gauss(int size)
{
	// ?? between 61/3 and 74/4 BHW
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	double delta = 5.0 / size;
	double x = (1 - size) / 2.0 * delta;
	double c = -M_PI * std::exp(1.0) / 10.0;
	double sum = 0;
	
	for (int i = start; i < stop; ++i)
	{
		windowFunction[i] = std::exp(c * x * x);
		x += delta;
		sum += windowFunction[i];
	}
	
	for (int i = start; i < stop; ++i)
		windowFunction[i] /= sum;
}

void FFT::rectangle(int size)
{
	int start = (static_cast<int>(windowFunction.size()) - size) / 2;
	int stop = (static_cast<int>(windowFunction.size()) + size) / 2;
	
	for (int i = start; i < stop; ++i)
		windowFunction[i] = 1.0 / static_cast<double>(size);
}

void FFT::setWindowFunction(int windowFunctionType, int support)
{
	if (support > windowSize)
		support = windowSize;
	
	switch (windowFunctionType)
	{
		case WND_NONE: 
			break;
		case WND_RECT:
			rectangle(support); break;
		case WND_HAMMING:
			hamming(support); break;
		case WND_HANNING:
			hanning(support); break;
		case WND_BH3:
			blackmanHarris3s(support); break;
		case WND_BH4:
			blackmanHarris4s(support); break;
		case WND_BH3MIN:
			blackmanHarris3sMin(support); break;
		case WND_BH4MIN:
			blackmanHarris4sMin(support); break;
		case WND_GAUSS:
			gauss(support); break;
		case WND_HANNINGZ:
			hanningz(support); break;
		default:
			windowFunctionType = WND_NONE;
			throw IllegalArgumentException("unknown window function specified");
	}
	
	this->windowFunctionType = windowFunctionType;
	calculateWindowFunctionSum();
}

int FFT::getTransformationType() const
{
	return transformationType;
}

int FFT::getWindowFunctionType() const
{
	return windowFunctionType;
}

void FFT::applyWindowFunction(std::vector<double>& data)
{
	if (windowFunctionType != WND_NONE)
	{
		for (size_t i = 0; i < data.size(); ++i)
			data[i] *= windowFunction[i];
	}
}

void FFT::calculateWindowFunctionSum()
{
	windowFunctionSum = 0;
	for(size_t i = 0; i < windowFunction.size(); ++i)
		windowFunctionSum += windowFunction[i];
}