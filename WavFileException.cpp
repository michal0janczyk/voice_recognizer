#include "WavFileException.h"

WavFileException::WavFileException(const char* exception) : exception(exception)
{
	// Empty
}

WavFileException::WavFileException(const std::string& exception) : exception(exception)
{
	// Empty
}

const char* WavFileException::what() const noexcept
{
	return exception.c_str();
}

std::ostream& operator<<(std::ostream& os, const WavFileException& wfe)
{
	return os << wfe.what();
}