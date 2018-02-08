#ifndef WAVFILEEXCEPTION_H
#define WAVFILEEXCEPTION_H

#include <exception>
#include <ostream>
#include <string>

class WavFileException : public std::exception
{
public:
	explicit WavFileException(const char* exception);
	explicit WavFileException(const std::string& exception);
	
	const char* what() const noexcept override;

private:
	friend std::ostream& operator<<(std::ostream& os, const WavFileException& wfe);
	
	std::string exception;
};

#endif /* WAVFILEEXCEPTION_H */