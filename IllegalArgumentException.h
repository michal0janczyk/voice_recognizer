#ifndef ILLEGALARGUMENTEXCEPTION_H
#define ILLEGALARGUMENTEXCEPTION_H

#include <exception>
#include <ostream>

class IllegalArgumentException : public std::exception
{
public:
	explicit IllegalArgumentException(const char* exception);
	explicit IllegalArgumentException(const std::string& exception);
	
	const char* what() const noexcept override;

private:
	friend std::ostream& operator<<(std::ostream& os, const IllegalArgumentException& iae);
	
	std::string exception;
};

#endif /* ILLEGALARGUMENTEXCEPTION_H */