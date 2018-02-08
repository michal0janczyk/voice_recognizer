#ifndef WAVFILE_H
#define WAVFILE_H

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "WavFileException.h"

class SoundProcessingTask3;

// C++ port of:
// Wav file IO class
// A.Greensted
// http://www.labbookpages.co.uk

// File format is based on the information from
// http://www.sonicspot.com/guide/wavefiles.html
// http://www.blitter.com/~russtopia/MIDI/~jglatt/tech/wave.htm

// Version 1.0

class WavFile
{
public:
	~WavFile();
	
	std::istream& ifstream();
	std::vector<int_least8_t> getDataBytes();
	std::vector<int_least8_t> getHeaderBytes();
	uint_least16_t getBlockAlign() const;
	uint_least16_t getNumChannels() const;
	int_least64_t getNumFrames() const;
	int_least64_t getFramesRemaining() const;
	uint_least32_t getSampleRate() const;
	uint_least32_t getByteRate() const;
	uint_least16_t getValidBits() const;
	static WavFile* newWavFile(const char* file, uint_least16_t numChannels, int_least64_t numFrames, uint_least16_t validBits, uint_least32_t sampleRate);
	static WavFile* openWavFile(const char* file);
	// Integer 32-bit
	int readFrames(std::vector<int_least32_t>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<int_least32_t>& sampleBuffer, int offset, int numFramesToRead);
	int readFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int offset, int numFramesToRead);
	int writeFrames(std::vector<int_least32_t>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<int_least32_t>& sampleBuffer, int offset, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int offset, int numFramesToWrite);
	// Long 64-bit
	int readFrames(std::vector<int_least64_t>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<int_least64_t>& sampleBuffer, int offset, int numFramesToRead);
	int readFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int offset, int numFramesToRead);
	int writeFrames(std::vector<int_least64_t>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<int_least64_t>& sampleBuffer, int offset, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int offset, int numFramesToWrite);
	// Double 64-bit
	int readFrames(std::vector<double>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<double>& sampleBuffer, int offset, int numFramesToRead);
	int readFrames(std::vector<std::vector<double>>& sampleBuffer, int numFramesToRead);
	int readFrames(std::vector<std::vector<double>>& sampleBuffer, int offset, int numFramesToRead);
	int writeFrames(std::vector<double>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<double>& sampleBuffer, int offset, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<double>>& sampleBuffer, int numFramesToWrite);
	int writeFrames(std::vector<std::vector<double>>& sampleBuffer, int offset, int numFramesToWrite);
	void close();
	void display();
	
private:
	friend class SoundProcessingTask3;
	
	enum IOState
	{
		READING,
		WRITING,
		CLOSED
	};
	
	/** Cannot instantiate WavFile directly, must either use newWavFile() or openWavFile() */
	WavFile();
	WavFile(const WavFile& wavFile) = delete;
	WavFile(WavFile&& wavFile) = delete;
	
	void operator=(const WavFile& wavFile) = delete;
	void operator=(WavFile&& wavFile) = delete;
	
	/** Get little endian data from local buffer */
	static int_least64_t getLE(std::vector<int_least8_t>& buffer, int pos, int numBytes);
	/** Put little endian data to local buffer */
	static void putLE(int_least64_t val, std::vector<int_least8_t>& buffer, int pos, int numBytes);
	/** Sample Writing */
	void writeSample(int_least64_t val);
	/** Sample Reading */
	int_least64_t readSample();
	
	constexpr static int BUFFER_SIZE = 4096;
	
	constexpr static int FMT_CHUNK_ID = 0x20746D66;
	constexpr static int DATA_CHUNK_ID = 0x61746164;
	constexpr static int RIFF_CHUNK_ID = 0x46464952;
	constexpr static int RIFF_TYPE_ID = 0x45564157;
	
	std::string filename;    /** Filename */
	IOState ioState;	     /** Specifies the IO State of the Wav File (used for snaity checking) */
	int bytesPerSample;	     /** Number of bytes required to store a single sample */
	int_least64_t numFrames; /** Number of frames within the data section */
	std::ofstream oStream;   /** Output stream used for writting data */
	std::ifstream iStream;   /** Input stream used for reading data */
	double floatScale;	     /** Scaling factor used for int <-> float conversion */
	double floatOffset;	     /** Offset factor used for int <-> float conversion */
	bool wordAlignAdjust;    /** Specify if an extra byte at the end of the data chunk is required for word alignment */
	
	// Wav Header
	uint_least16_t numChannels; /** 2 bytes unsigned, 0x0001 (1) to 0xFFFF (65,535) */
	uint_least32_t sampleRate;  /** 4 bytes unsigned, 0x00000001 (1) to 0xFFFFFFFF (4,294,967,295) */
	uint_least32_t byteRate;
	uint_least16_t blockAlign;  /** 2 bytes unsigned, 0x0001 (1) to 0xFFFF (65,535) */
	uint_least16_t validBits;   /** 2 bytes unsigned, 0x0002 (2) to 0xFFFF (65,535) */
	
	// Buffering
	std::vector<int_least8_t> headerBytes;
	std::vector<int_least8_t> dataBytes;
	std::vector<int_least8_t> buffer; /** Local buffer used for IO */
	int bufferPointer;		          /** Points to the current position in local buffer */
	int bytesRead;			          /** Bytes read after last read into local buffer */
	int_least64_t frameCounter;	      /** Current number of frames read or written */
};

#endif /* WAVFILE_H */