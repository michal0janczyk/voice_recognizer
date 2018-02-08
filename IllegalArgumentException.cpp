#include "IllegalArgumentException.h"

IllegalArgumentException::IllegalArgumentException(const char* exception) : exception(exception)
{
	// Empty
}

IllegalArgumentException::IllegalArgumentException(const std::string& exception) : exception(exception)
{
	// Empty
}

const char* IllegalArgumentException::what() const noexcept
{
	return exception.c_str();
}

std::ostream& operator<<(std::ostream& os, const IllegalArgumentException& iae)
{
	return os << iae.what();
}