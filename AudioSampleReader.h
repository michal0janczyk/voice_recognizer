#ifndef AUDIOSAMPLEREADER_H
#define AUDIOSAMPLEREADER_H

#include <cstdint>
#include <string>
#include <vector>

#include "IllegalArgumentException.h"

class WavFile;

class AudioSampleReader
{
public:
	explicit AudioSampleReader(const char* file);
	AudioSampleReader(const AudioSampleReader& asr);
	~AudioSampleReader();
	
	void operator=(const AudioSampleReader& asr);
	
	int_least64_t getSampleCount() const;
	void getInterleavedSamples(int_least64_t begin, int_least64_t end, std::vector<double>& samples, std::vector<int_least16_t>& samplesInShort);
	void getChannelSamples(int channel, const std::vector<double>& interleavedSamples, std::vector<double>& channelSamples);
	void getStereoSamples(std::vector<double>& leftSamples, std::vector<double>& rightSamples, std::vector<int_least16_t>& samplesInShort);
	void decodeBytes(const std::vector<int_least8_t>& audioBytes, std::vector<double>& audioSamples, std::vector<int_least16_t>& samplesInShort);
	
private:
	WavFile* wav;
	std::string filename;
};

#endif /* AUDIOSAMPLEREADER_H */