#include "SoundProcessingTask3.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "DTW.h"
#include "MFCC.h"
#include "WavFile.h"

SoundProcessingTask3::SoundProcessingTask3()
{
	files = {"jeden", "dwa", "trzy", "test"};
	audioFileLength = 0;
	sampleRate = 44100.0f;
	windowSize = 512;
	numberCoefficients = 12;
	useFirstCoefficient = true;
	numberFilters = 30;
	nbSamples = 0;
}

void SoundProcessingTask3::load()
{
	for (size_t i = 0; i < files.size(); i++)
	{
		std::ifstream ifs(files[i] + ".txt", std::ifstream::in);
		if (!ifs.is_open())
		{
			std::string ex = files[i];
			ex += ".txt";
			ex += " doesn't exist";
			throw std::runtime_error(ex);
		}
		std::vector<double> temps;
		
		std::istringstream iss;
		std::string token;
		
		for (std::string line; std::getline(ifs, line); )
		{
			iss.clear();
			iss.str(line);
			
			while(std::getline(iss, token, ' '))
			{
				double token1 = std::stod(token);
				temps.push_back(token1);
			}
		}
		ifs.close();
		std::vector<double> array(temps.size());
		
		for (size_t j = 0; j < temps.size(); j++)
			array[j] = temps.at(j);
		database.push_back(array);
		
		//std::cout << database.size() << std::endl;
		//for(auto& currentArray : database)
		//	std::copy(currentArray.begin(), currentArray.end(), std::ostream_iterator<double>(std::cout, " "));
	}
}

void SoundProcessingTask3::dtw()
{
	double min = 999.0;
	std::string result = "?";
	for (size_t i = 0; i < database.size(); i++)
	{
		DTW dtw{mfcc, database.at(i)};
		if (dtw.getDistance() < min)
		{
			result = files[i].substr(0, 1);
			min = dtw.getDistance();
		}
		std::cout << files[i].substr(0, 1) << ":  " << std::to_string(dtw.getDistance()) << '\n';
	}
	std::cout << "RESULT: " << result << std::endl;
}

void SoundProcessingTask3::compare()
{
	std::cout << "Name of file: ";
	std::string option;
	std::cin >> option;
	std::string lines;
	
	std::ifstream ifs(option + ".txt", std::ifstream::in);
	if (!ifs.is_open())
	{
		std::string ex = option;
		ex += ".txt";
		ex += " doesn't exist";
		throw std::runtime_error(ex);
	}
	std::vector<double> temps;
	
	std::istringstream iss;
	std::string token;
	
	for (std::string line; std::getline(ifs, line); )
	{
		iss.clear();
		iss.str(line);
		
		while(std::getline(iss, token, ' '))
		{
			double token1 = std::stod(token);
			temps.push_back(token1);
		}
	}
	ifs.close();
	
	std::vector<double> array(temps.size());
	for (size_t i = 0; i < temps.size(); i++)
		array[i] = temps.at(i);
	
	for (size_t i = 0; i < database.size(); i++)
	{
		DTW dtw{array, database.at(i)};
		std::cout << dtw << '\n';
	}
	std::cout << std::endl;
}

void SoundProcessingTask3::create()
{
	MFCC kons{
		sampleRate,
		windowSize,
		numberCoefficients,
		useFirstCoefficient,
		20.0, 
		16000.0,
		numberFilters
	};
	
	std::cout << "Name of file: ";
	std::string fileName;
	std::cin >> fileName;
	openFile(fileName.c_str());
	toneBuffer = cutBuffer();
	coof = kons.process(toneBuffer);
	std::cout << toneBuffer.size() << std::endl;
	mfcc = normMFCC(coof);
	saveFile(fileName.c_str(), mfcc);
}

void SoundProcessingTask3::saveFile(const char* fileName, const std::vector<double>& array)
{
	try
	{
		std::string file = fileName;
		std::ofstream ofs(file + ".txt", std::ofstream::out);
		
		auto precision = std::cout.precision();
		for (size_t i = 0; i < array.size(); i++)
			ofs << std::setprecision(16) << array[i] << '\n';
		ofs.close();
		std::cout << std::setprecision(precision);
	}
	catch (...)
	{
		std::cerr << "No such file exists.";
	}
}

void SoundProcessingTask3::openFile(const char* fileName)
{
	std::string file = fileName;
	file += ".wav";
	AudioSampleReader sampleReader{file.c_str()};
	nbSamples = sampleReader.getSampleCount();
	toneBuffer.resize(static_cast<int_least32_t>(nbSamples));
	
	samplesInShort.resize(static_cast<int_least32_t>(nbSamples));
	sampleReader.getInterleavedSamples(
		0,
		nbSamples,
		toneBuffer,
		samplesInShort
	);
}

std::vector<double> SoundProcessingTask3::cutBuffer()
{
	std::vector<double> outBuffer;
	if (toneBuffer.size() % 512 != 0)
	{
		int cut = static_cast<int>(toneBuffer.size()) - static_cast<int>(toneBuffer.size()) % 512;
		std::copy(
			toneBuffer.begin(),
			toneBuffer.begin() + cut,
			std::back_inserter(outBuffer)
		);
		return outBuffer;
	}
	return toneBuffer;
}

void SoundProcessingTask3::printMFCC(const std::vector<std::vector<double>>& array)
{
	for (size_t i = 0; i < array.size(); i++)
	{
		std::cout << std::to_string(i) << ": ";
		for (size_t j = 0; j < array[i].size(); j++)
		{
			double roundOff = std::round(array[i][j] * 100.0) / 100.0;
			std::cout << roundOff << ' ';
		}
	}
	std::cout << std::endl;
}

void SoundProcessingTask3::printMFCC(const std::vector<double>& array)
{
	for (size_t i = 0; i < array.size(); i++)
	{
		double roundOff = std::round(array[i] * 100.0) / 100.0;
		std::cout << std::to_string(i) << ": " << roundOff << '\n';
	}
}

std::vector<double> SoundProcessingTask3::normMFCC(const std::vector<std::vector<double>>& array)
{
	std::vector<double> output(array.size());
	for (size_t i = 0; i < array.size(); i++)
		output[i] = norm(array[i]);
	return output;
}

double SoundProcessingTask3::norm(const std::vector<double>& data)
{
	double sum = 0.0;
	for (size_t i = 0; i < data.size(); i++)
		sum = sum + (data[i] * data[i]);
	return std::sqrt(sum);
}

std::vector<std::vector<double>> SoundProcessingTask3::chunkArray(const std::vector<double>& array)
{
	int numOfChunks = static_cast<int>(std::ceil(static_cast<double>(static_cast<int>(array.size()) / SoundProcessingTask3::windowSize)));
	std::vector<std::vector<double>> output;
	output.resize(numOfChunks, std::vector<double>());
	
	for (int i = 0; i < numOfChunks; ++i)
	{
		int start = i * windowSize;
		int length = std::min(static_cast<int>(array.size()) - start, windowSize);
		
		std::vector<double> temp(length);
		std::copy(
			array.begin() + start,
			array.begin() + start + length,
			std::back_inserter(output[i])
		);	
	}
	
	return output;
}