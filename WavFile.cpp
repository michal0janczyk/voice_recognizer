#include "WavFile.h"

#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>

WavFile::~WavFile()
{
	// Empty
}

std::istream& WavFile::ifstream()
{
	return iStream;
}

std::vector<int_least8_t> WavFile::getHeaderBytes()
{
	return headerBytes;
}

std::vector<int_least8_t> WavFile::getDataBytes()
{
	return dataBytes;
}

uint_least16_t WavFile::getBlockAlign() const
{
	return blockAlign;
}

uint_least16_t WavFile::getNumChannels() const
{
	return numChannels;
}

int_least64_t WavFile::getNumFrames() const
{
	return numFrames;
}

int_least64_t WavFile::getFramesRemaining() const
{
	return numFrames - frameCounter;
}

uint_least32_t WavFile::getSampleRate() const
{
	return sampleRate;
}

uint_least32_t WavFile::getByteRate() const
{
	return byteRate;
}

uint_least16_t WavFile::getValidBits() const
{
	return validBits;
}

WavFile* WavFile::newWavFile(const char* file, uint_least16_t numChannels, int_least64_t numFrames, uint_least16_t validBits, uint_least32_t sampleRate)
{
	// Instantiate new Wavfile and initialise
	WavFile* wavFile = new WavFile();
	wavFile->filename = file;
	wavFile->numChannels = numChannels;
	wavFile->numFrames = numFrames;
	wavFile->sampleRate = sampleRate;
	wavFile->bytesPerSample = (validBits + 7) / 8;
	wavFile->blockAlign = wavFile->bytesPerSample * numChannels;
	wavFile->validBits = validBits;
	
	// Sanity check arguments
	if (numChannels < 1 || numChannels > 65535)
	{
		delete wavFile;
		throw WavFileException("Illegal number of channels, valid range 1 to 65536");
	}
	if (numFrames < 0)
	{
		delete wavFile;
		throw WavFileException("Number of frames must be positive");
	}
	if (validBits < 2 || validBits > 65535)
	{
		delete wavFile;
		throw WavFileException("Illegal number of valid bits, valid range 2 to 65536");
	}
	if (sampleRate < 0)
	{
		delete wavFile;
		throw WavFileException("Sample rate must be positive");
	}
	
	// Create output stream for writing data
	wavFile->oStream.open(file, std::ofstream::out | std::ofstream::binary);
	
	// Calculate the chunk sizes
	int_least64_t dataChunkSize = wavFile->blockAlign * numFrames;
	int_least64_t mainChunkSize = 4  +	// Riff Type
								  8  +	// Format ID and size
								  16 +	// Format data
								  8  + 	// Data ID and size
								  dataChunkSize;
	
	// Chunks must be word aligned, so if odd number of audio data bytes
	// adjust the main chunk size
	if (dataChunkSize % 2 == 1)
	{
		mainChunkSize += 1;
		wavFile->wordAlignAdjust = true;
	}
	else
		wavFile->wordAlignAdjust = false;
	
	// Set the main chunk size
	putLE(RIFF_CHUNK_ID, wavFile->buffer, 0, 4);
	putLE(mainChunkSize, wavFile->buffer, 4, 4);
	putLE(RIFF_TYPE_ID,	wavFile->buffer, 8, 4);
	
	// Write out the header
	wavFile->oStream.write(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 12);
	
	// Put format data in buffer
	long averageBytesPerSecond = sampleRate * wavFile->blockAlign;
	
	putLE(FMT_CHUNK_ID,			 wavFile->buffer, 0, 4);  // Chunk ID
	putLE(16,					 wavFile->buffer, 4, 4);  // Chunk Data Size
	putLE(1,					 wavFile->buffer, 8, 2);  // Compression Code (Uncompressed)
	putLE(numChannels,			 wavFile->buffer, 10, 2); // Number of channels
	putLE(sampleRate,			 wavFile->buffer, 12, 4); // Sample Rate
	putLE(averageBytesPerSecond, wavFile->buffer, 16, 4); // Average Bytes Per Second
	putLE(wavFile->blockAlign,	 wavFile->buffer, 20, 2); // Block Align
	putLE(validBits,			 wavFile->buffer, 22, 2); // Valid Bits
	
	// Write Format Chunk
	wavFile->oStream.write(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 24);
	
	// Start Data Chunk
	putLE(DATA_CHUNK_ID, wavFile->buffer, 0, 4); // Chunk ID
	putLE(dataChunkSize, wavFile->buffer, 4, 4); // Chunk Data Size
	
	// Write Format Chunk
	wavFile->oStream.write(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 8);
	
	// Calculate the scaling factor for converting to a normalised double
	if (wavFile->validBits > 8)
	{
		// If more than 8 validBits, data is signed
		// Conversion required multiplying by magnitude of max positive value
		wavFile->floatOffset = 0;
		wavFile->floatScale = static_cast<double>(std::numeric_limits<int_least64_t>::max() >> (64 - wavFile->validBits));
	}
	else
	{
		// Else if 8 or less validBits, data is unsigned
		// Conversion required dividing by max positive value
		wavFile->floatOffset = 1;
		wavFile->floatScale = 0.5 * ((1 << wavFile->validBits) - 1);
	}
	
	// Finally, set the IO State
	wavFile->bufferPointer = 0;
	wavFile->bytesRead = 0;
	wavFile->frameCounter = 0;
	wavFile->ioState = IOState::WRITING;
	
	return wavFile;
}

WavFile* WavFile::openWavFile(const char* file)
{
	// Instantiate new Wavfile and store the file reference
	WavFile* wavFile = new WavFile();
	wavFile->filename = file;
	
	// Create a new file input stream for reading file data
	wavFile->iStream.open(file, std::ifstream::in | std::ifstream::binary);
	if (!wavFile->iStream.is_open())
	{
		std::string ex = wavFile->filename;
		ex += " doesn't exist";
		delete wavFile;
		throw WavFileException(ex);
	}
	wavFile->iStream.seekg(0, wavFile->iStream.end);
	size_t fileLength = static_cast<size_t>(wavFile->iStream.tellg());
	wavFile->iStream.seekg(0, wavFile->iStream.beg);
	
	wavFile->headerBytes.resize(44u);
	wavFile->iStream.read(static_cast<char*>(static_cast<void*>(wavFile->headerBytes.data())), 44u);
	
	wavFile->dataBytes.resize(fileLength - 44u);
	wavFile->iStream.read(static_cast<char*>(static_cast<void*>(wavFile->dataBytes.data())), fileLength - 44u);
	
	wavFile->iStream.seekg(0, wavFile->iStream.beg);
	
	// Read the first 12 bytes of the file
	wavFile->iStream.read(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 12);
	if (wavFile->iStream.gcount() != 12)
	{
		delete wavFile;
		throw WavFileException("Not enough wav file bytes for header");
	}
	
	// Extract parts from the header
	int_least64_t riffChunkID = getLE(wavFile->buffer, 0, 4);
	int_least64_t chunkSize = getLE(wavFile->buffer, 4, 4);
	int_least64_t riffTypeID = getLE(wavFile->buffer, 8, 4);
	
	// Check the header bytes contains the correct signature
	if (riffChunkID != RIFF_CHUNK_ID)
	{
		delete wavFile;
		throw WavFileException("Invalid Wav Header data, incorrect riff chunk ID");
	}
	if (riffTypeID != RIFF_TYPE_ID)
	{
		delete wavFile;
		throw WavFileException("Invalid Wav Header data, incorrect riff type ID");
	}
	
	// Check that the file size matches the number of bytes listed in header
	if (static_cast<int>(fileLength) != chunkSize + 8)
	{
		std::string ex = "Header chunk size (";
		ex += std::to_string(chunkSize);
		ex += ") does not match file size (";
		ex += std::to_string(fileLength);
		ex += ")";
		delete wavFile;
		throw WavFileException(ex);
	}
	
	bool foundFormat = false;
	bool foundData = false;
	
	// Search for the Format and Data Chunks
	while (true)
	{
		// Read the first 8 bytes of the chunk (ID and chunk size)
		wavFile->iStream.read(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 8);
		auto bytesRead = wavFile->iStream.gcount();
		if (bytesRead == -1)
		{
			delete wavFile;
			throw WavFileException("Reached end of file without finding format chunk");
		}
		if (bytesRead != 8)
		{
			delete wavFile;
			throw WavFileException("Could not read chunk header");
		}
		
		// Extract the chunk ID and Size
		int_least64_t chunkID = getLE(wavFile->buffer, 0, 4);
		chunkSize = getLE(wavFile->buffer, 4, 4);
		
		// Word align the chunk size
		// chunkSize specifies the number of bytes holding data. However,
		// the data should be word aligned (2 bytes) so we need to calculate
		// the actual number of bytes in the chunk
		int_least64_t numChunkBytes = (chunkSize % 2 == 1) ? chunkSize + 1 : chunkSize;
		
		if (chunkID == FMT_CHUNK_ID)
		{
			// Flag that the format chunk has been found
			foundFormat = true;
			
			// Read in the header info
			wavFile->iStream.read(static_cast<char*>(static_cast<void*>(wavFile->buffer.data())), 16);
			bytesRead = wavFile->iStream.gcount();
			
			// Check this is uncompressed data
			int compressionCode = static_cast<int>(getLE(wavFile->buffer, 0, 2));
			if (compressionCode != 1)
			{
				delete wavFile;
				std::string ex = "Compression Code ";
				ex += compressionCode;
				ex += " not supported";
				throw WavFileException(ex);
			}
			
			// Extract the format information
			wavFile->numChannels = static_cast<uint_least16_t>(getLE(wavFile->buffer, 2, 2));
			wavFile->sampleRate = static_cast<uint_least32_t>(getLE(wavFile->buffer, 4, 4));
			wavFile->byteRate = static_cast<uint_least32_t>(getLE(wavFile->buffer, 8, 4));
			wavFile->blockAlign = static_cast<uint_least16_t>(getLE(wavFile->buffer, 12, 2));
			wavFile->validBits = static_cast<uint_least16_t>(getLE(wavFile->buffer, 14, 2));
			
			if (wavFile->numChannels == 0)
			{
				delete wavFile;
				throw WavFileException("Number of channels specified in header is equal to zero");
			}
			if (wavFile->blockAlign == 0)
			{
				delete wavFile;
				throw WavFileException("Block Align specified in header is equal to zero");
			}
			if (wavFile->validBits < 2)
			{
				delete wavFile;
				throw WavFileException("Valid Bits specified in header is less than 2");
			}
			if (wavFile->validBits > 64)
			{
				delete wavFile;
				throw WavFileException("Valid Bits specified in header is greater than 64, this is greater than a long can hold");
			}
			
			// Calculate the number of bytes required to hold 1 sample
			wavFile->bytesPerSample = (wavFile->validBits + 7) / 8;
			if (wavFile->bytesPerSample * wavFile->numChannels != wavFile->blockAlign)
			{
				delete wavFile;
				throw WavFileException("Block Align does not agree with bytes required for validBits and number of channels");
			}
			
			// Account for number of format bytes and then skip over
			// any extra format bytes
			numChunkBytes -= 16;
			if (numChunkBytes > 0)
				wavFile->iStream.ignore(numChunkBytes);
		}
		else if (chunkID == DATA_CHUNK_ID)
		{
			// Check if we've found the format chunk,
			// If not, throw an exception as we need the format information
			// before we can read the data chunk
			if (foundFormat == false)
			{
				delete wavFile;
				throw WavFileException("Data chunk found before Format chunk");
			}
			
			// Check that the chunkSize (wav data length) is a multiple of the
			// block align (bytes per frame)
			if (chunkSize % wavFile->blockAlign != 0)
			{
				delete wavFile;
				throw WavFileException("Data Chunk size is not multiple of Block Align");
			}
			
			// Calculate the number of frames
			wavFile->numFrames = chunkSize / wavFile->blockAlign;
			
			// Flag that we've found the wave data chunk
			foundData = true;
			
			break;
		}
		else
		{
			// If an unknown chunk ID is found, just skip over the chunk data
			wavFile->iStream.ignore(numChunkBytes);
		}
	}
	
	// Throw an exception if no data chunk has been found
	if (!foundData)
	{
		delete wavFile;
		throw WavFileException("Did not find a data chunk");
	}
	
	// Calculate the scaling factor for converting to a normalised double
	if (wavFile->validBits > 8)
	{
		// If more than 8 validBits, data is signed
		// Conversion required dividing by magnitude of max negative value
		wavFile->floatOffset = 0;
		wavFile->floatScale = 1 << (wavFile->validBits - 1);
	}
	else
	{
		// Else if 8 or less validBits, data is unsigned
		// Conversion required dividing by max positive value
		wavFile->floatOffset = -1;
		wavFile->floatScale = 0.5 * ((1 << wavFile->validBits) - 1);
	}
	
	wavFile->bufferPointer = 0;
	wavFile->bytesRead = 0;
	wavFile->frameCounter = 0;
	wavFile->ioState = IOState::READING;
	
	return wavFile;
}

int WavFile::readFrames(std::vector<int_least32_t>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<int_least32_t>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			sampleBuffer[offset] = static_cast<int_least32_t>(readSample());
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::readFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			sampleBuffer[c][offset] = static_cast<int_least32_t>(readSample());
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::writeFrames(std::vector<int_least32_t>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<int_least32_t>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			writeSample(sampleBuffer[offset]);
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToWrite;
}

int WavFile::writeFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<std::vector<int_least32_t>>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			writeSample(sampleBuffer[c][offset]);
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToWrite;
}

int WavFile::readFrames(std::vector<int_least64_t>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<int_least64_t>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			sampleBuffer[offset] = readSample();
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::readFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			sampleBuffer[c][offset] = readSample();
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::writeFrames(std::vector<int_least64_t>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<int_least64_t>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			writeSample(sampleBuffer[offset]);
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToWrite;
}

int WavFile::writeFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<std::vector<int_least64_t>>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			writeSample(sampleBuffer[c][offset]);
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToWrite;
}

int WavFile::readFrames(std::vector<double>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<double>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			sampleBuffer[offset] = floatOffset + static_cast<double>(readSample()) / floatScale;
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::readFrames(std::vector<std::vector<double>>& sampleBuffer, int numFramesToRead)
{
	return readFrames(sampleBuffer, 0, numFramesToRead);
}

int WavFile::readFrames(std::vector<std::vector<double>>& sampleBuffer, int offset, int numFramesToRead)
{
	if (ioState != IOState::READING)
		throw std::runtime_error("Cannot read from WavFile instance");
	
	for (int f = 0; f < numFramesToRead; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			sampleBuffer[c][offset] = floatOffset + static_cast<double>(readSample()) / floatScale;
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToRead;
}

int WavFile::writeFrames(std::vector<double>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<double>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
		{
			writeSample(static_cast<int_least64_t>(floatScale * (floatOffset + sampleBuffer[offset])));
			offset++;
		}
		
		frameCounter++;
	}
	
	return numFramesToWrite;
}

int WavFile::writeFrames(std::vector<std::vector<double>>& sampleBuffer, int numFramesToWrite)
{
	return writeFrames(sampleBuffer, 0, numFramesToWrite);
}

int WavFile::writeFrames(std::vector<std::vector<double>>& sampleBuffer, int offset, int numFramesToWrite)
{
	if (ioState != IOState::WRITING)
		throw std::runtime_error("Cannot write to WavFile instance");
	
	for (int f = 0; f < numFramesToWrite; f++)
	{
		if (frameCounter == numFrames)
			return f;
		
		for (int c = 0; c < numChannels; c++)
			writeSample(static_cast<int_least64_t>(floatScale * (floatOffset + sampleBuffer[c][offset])));
		
		offset++;
		frameCounter++;
	}
	
	return numFramesToWrite;
}

void WavFile::close()
{
	// Close the input stream
	if (iStream.is_open())
		iStream.close();
	
	if (oStream.is_open()) 
	{
		// Write out anything still in the local buffer
		if (bufferPointer > 0)
			oStream.write(static_cast<char*>(static_cast<void*>(buffer.data())), bufferPointer);
		
		// If an extra byte is required for word alignment, add it to the end
		if (wordAlignAdjust)
		{
			char zero = static_cast<char>(0);
			const char* zeroPtr = &zero;
			oStream.write(zeroPtr, 1);
		}
		
		// Close the stream
		oStream.close();
	}

	// Flag that the stream is closed
	ioState = IOState::CLOSED;
}

void WavFile::display()
{
	std::cout << "File: " << filename << '\n';
	std::cout << "Channels: " << std::to_string(numChannels) << '\n';
	std::cout << "Frames: " << std::to_string(numFrames) << std::endl;
}

WavFile::WavFile()
{
	buffer.resize(BUFFER_SIZE);
	bytesPerSample = 0;
	numFrames = 0;
	floatScale = 0.0;
	floatOffset = 0.0;
	wordAlignAdjust = false;
	numChannels = 0;
	sampleRate = 0;
	blockAlign = 0;
	validBits = 0;
	bufferPointer = 0;
	bytesRead = 0;
	frameCounter = 0;
}

int_least64_t WavFile::getLE(std::vector<int_least8_t>& buffer, int pos, int numBytes)
{
	numBytes--;
	pos += numBytes;
	
	int_least64_t val = buffer[pos] & 0xFF;
	for (int b = 0; b < numBytes; b++)
		val = (val << 8) + (buffer[--pos] & 0xFF);
	
	return val;
}

void WavFile::putLE(int_least64_t val, std::vector<int_least8_t>& buffer, int pos, int numBytes)
{
	for (int b = 0; b < numBytes; b++)
	{
		buffer[pos] = static_cast<int_least8_t>(val & 0xFF);
		val >>= 8;
		pos++;
	}
}

void WavFile::writeSample(int_least64_t val)
{
	for (int b = 0; b < bytesPerSample; b++)
	{
		if (bufferPointer == BUFFER_SIZE)
		{
			auto test = static_cast<char*>(static_cast<void*>(buffer.data()));
			oStream.write(test, BUFFER_SIZE);
			bufferPointer = 0;
		}
		buffer[bufferPointer] = static_cast<int_least8_t>(val & 0xFF);
		val >>= 8;
		bufferPointer++;
	}
}

int_least64_t WavFile::readSample()
{
	int_least64_t val = 0;
	
	for (int b = 0; b < bytesPerSample; b++)
	{
		if (bufferPointer == bytesRead)
		{
			iStream.read(static_cast<char*>(static_cast<void*>(buffer.data())), BUFFER_SIZE);
			if (iStream.fail() && !iStream.eof())
				throw WavFileException("Not enough data available");
			bytesRead = static_cast<int>(iStream.gcount());
			bufferPointer = 0;
		}
		
		int v = buffer[bufferPointer];
		if (b < bytesPerSample - 1 || bytesPerSample == 1)
			v &= 0xFF;
		val += v << (b * 8);
		
		bufferPointer++;
	}
	
	return val;
}