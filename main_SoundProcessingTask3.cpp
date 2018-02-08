#include <exception>
#include <iostream>
#include <stdexcept>

#include "SoundProcessingTask3.h"
#include "IllegalArgumentException.h"
#include "WavFileException.h"

int main()
{
	try
	{
		SoundProcessingTask3 task;
		task.create();
		task.load();
		//task.compare();
		task.dtw();
	}
	catch(std::exception& ex)
	{
		std::cerr << ex.what();
		return 1;
	}
	
	return 0;
}