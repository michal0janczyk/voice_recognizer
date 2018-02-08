#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "WavFile.h"

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		std::cerr << "Must supply filename";
		return 1;
	}
	
	try
	{
		constexpr int BUF_SIZE = 5001;
		
		for (int i = 1; i < argc; ++i)
		{
			const char* filename = argv[i];
			
			WavFile* readWavFile = WavFile::openWavFile(filename);
			readWavFile->display();
			
			int_least64_t numFrames = readWavFile->getNumFrames();
			uint_least16_t numChannels = readWavFile->getNumChannels();
			uint_least16_t validBits = readWavFile->getValidBits();
			uint_least32_t sampleRate = readWavFile->getSampleRate();
			
			WavFile* writeWavFile = WavFile::newWavFile("out.wav", numChannels, numFrames, validBits, sampleRate);
			
			//std::vector<int_least32_t> buffer(BUF_SIZE * numChannels);
			//std::vector<int_least64_t> buffer(BUF_SIZE * numChannels);
			std::vector<double> buffer(BUF_SIZE * numChannels);
			
			int framesRead = 0;
			int framesWritten = 0;
			
			do
			{
				framesRead = readWavFile->readFrames(buffer, BUF_SIZE);
				framesWritten = writeWavFile->writeFrames(buffer, BUF_SIZE);
				std::printf("%d %d\n", framesRead, framesWritten);
			}
			while (framesRead != 0);
			
			readWavFile->close();
			writeWavFile->close();
			
			delete writeWavFile;
			delete readWavFile;
		}
		
		WavFile* writeWavFile = WavFile::newWavFile("out2.wav", 1, 10, 23, 44100);
		//std::vector<int_least32_t> buffer(10);
		//std::vector<int_least64_t> buffer(10);
		std::vector<double> buffer(10);
		writeWavFile->writeFrames(buffer, 10);
		writeWavFile->close();
		
		delete writeWavFile;
	}
	catch (std::exception& ex)
	{
		std::cerr << ex.what();
		return 1;
	}
	
	return 0;
}