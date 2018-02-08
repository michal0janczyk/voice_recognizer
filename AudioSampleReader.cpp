#include "AudioSampleReader.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "WavFile.h"

AudioSampleReader::AudioSampleReader(const char* file)
{
	wav = WavFile::openWavFile(file);
	filename = file;
}

AudioSampleReader::AudioSampleReader(const AudioSampleReader& asr)
{
	wav = WavFile::openWavFile(asr.filename.c_str());
	filename = asr.filename;
}

AudioSampleReader::~AudioSampleReader()
{
	delete wav;
}

void AudioSampleReader::operator=(const AudioSampleReader& asr)
{
	wav = WavFile::openWavFile(asr.filename.c_str());
	filename = asr.filename;
}

int_least64_t AudioSampleReader::getSampleCount() const
{
	//int_least64_t total = (audioInputStream.getFrameLength() * format.getFrameSize() * 8) / format.getSampleSizeInBits();
	//return total / format.getChannels();
	return (wav->getNumFrames() * wav->getBlockAlign() * 8) / wav->getValidBits() / wav->getNumChannels();
}

void AudioSampleReader::getInterleavedSamples(int_least64_t begin, int_least64_t end, std::vector<double>& samples, std::vector<int_least16_t>& samplesInShort)
{
	int_least64_t nbSamples = end - begin;
	int_least64_t nbBytes = nbSamples * (wav->getValidBits() / 8) * wav->getNumChannels();
	if (nbBytes > std::numeric_limits<int_least32_t>::max())
		throw IllegalArgumentException("too many samples");
	// Allocate a byte buffer
	std::vector<int_least8_t> inBuffer(static_cast<int_least32_t>(nbBytes));
	// Read bytes from audio file
	std::vector<int_least8_t> bytes = wav->getDataBytes();
	std::copy(
		bytes.begin() + static_cast<size_t>((begin * static_cast<int_least64_t>(wav->getValidBits() / 8))),
		bytes.begin() + static_cast<size_t>((begin + (end - begin) * static_cast<int_least64_t>(wav->getValidBits() / 8))),
		inBuffer.begin()
	);
	
	// Decode bytes into samples. Supported encodings are:
	// PCM-SIGNED, PCM-UNSIGNED, A-LAW, U-LAW
	decodeBytes(inBuffer, samples, samplesInShort);
}

void AudioSampleReader::getChannelSamples(int channel, const std::vector<double>& interleavedSamples, std::vector<double>& channelSamples)
{
	uint_least16_t nbChannels = wav->getNumChannels();
	for (int i = 0; i < static_cast<int>(channelSamples.size()); i++)
		channelSamples[i] = interleavedSamples[nbChannels * i + channel];
}

void AudioSampleReader::getStereoSamples(std::vector<double>& leftSamples, std::vector<double>& rightSamples, std::vector<int_least16_t>& samplesInShort)
{
	int_least64_t sampleCount = getSampleCount();
	std::vector<double> interleavedSamples(static_cast<int_least32_t>(sampleCount * 2));
	getInterleavedSamples(0, sampleCount, interleavedSamples, samplesInShort);
	for (size_t i = 0; i < leftSamples.size(); i++)
	{
		leftSamples[i] = interleavedSamples[2 * i];
		rightSamples[i] = interleavedSamples[2 * i + 1];
	}
}

void AudioSampleReader::decodeBytes(const std::vector<int_least8_t>& audioBytes, std::vector<double>& audioSamples, std::vector<int_least16_t>& samplesInShort)
{
	int sampleSizeInBytes = wav->getValidBits() / 8;
	std::vector<int> sampleBytes(sampleSizeInBytes);
	int k = 0; // Index in audioBytes
	for (size_t i = 0; i < audioSamples.size(); i++)
	{
		// Collect sample byte in big-endian order
		bool bigEndian = false;
		if (bigEndian)
		{
			// Bytes start with MSB - most significant bit
			for (int j = 0; j < sampleSizeInBytes; j++)
				sampleBytes[j] = audioBytes[k++];
		}
		else
		{
			// Bytes start with LSB
			for (int j = sampleSizeInBytes - 1; j >= 0; j--)
			{
				sampleBytes[j] = audioBytes[k++];
				if (sampleBytes[j] != 0)
					j = j + 0;
			}
		}
		samplesInShort.resize(sampleBytes.size());
		for (size_t n = 0; n < sampleBytes.size(); n += 2)
			samplesInShort[n] = static_cast<int_least16_t>(sampleBytes[n] | sampleBytes[n + 1] << 8);
		
		// Get integer value from bytes
		int ival = 0;
		for (int j = 0; j < sampleSizeInBytes; j++)
		{
			ival += sampleBytes[j];
			if (j < sampleSizeInBytes - 1)
				ival <<= 8;
		}
		// Decode value
		double ratio = std::pow(2.0, static_cast<double>(wav->getValidBits() - 1));
		double val = (static_cast<double>(ival)) / ratio;
		audioSamples[i] = val;
	}
}