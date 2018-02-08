#ifndef SOUNDPROCESSINGTASK3_H
#define SOUNDPROCESSINGTASK3_H

#include <cstdint>
#include <string>
#include <vector>

#include "AudioSampleReader.h"

class SoundProcessingTask3
{
public:
	SoundProcessingTask3();
	
	void load();
	void dtw();
	void compare();
	void create();
	void saveFile(const char* fileName, const std::vector<double>& array);
	void openFile(const char* fileName);
	std::vector<double> cutBuffer();
	void printMFCC(const std::vector<std::vector<double>>& array);
	void printMFCC(const std::vector<double>& array);
	std::vector<double> normMFCC(const std::vector<std::vector<double>>& array);
	double norm(const std::vector<double>& data);
	std::vector<std::vector<double>> chunkArray(const std::vector<double>& array);
	
private:
	std::vector<double> toneBuffer;
	std::vector<std::string> files;
	
	std::string fileName;
	int_least64_t audioFileLength;
	float sampleRate;
	int windowSize;
	int numberCoefficients;
	bool useFirstCoefficient;
	int numberFilters;
	
	int_least64_t nbSamples;
	std::vector<double> samples;
	std::vector<std::vector<double>> coof;
	std::vector<std::vector<double>> database;
	std::vector<double> mfcc;
	std::vector<int_least16_t> samplesInShort;
};

#endif /* SOUNDPROCESSINGTASK3_H */