cmake_minimum_required(VERSION 3.1.3 FATAL_ERROR)

project(Complex)

add_executable(Complex
	Complex.h Complex.cpp
	main_Complex.cpp
)

project(DTW)

add_executable(DTW DTW.h DTW.cpp main_DTW.cpp)

project(WavFile)

add_executable(WavFile
	WavFileException.h WavFileException.cpp
	WavFile.h WavFile.cpp 
	main_WavFile.cpp
)

project(SoundProcessingTask3)

add_executable(SoundProcessingTask3
	DTW.h DTW.cpp
	IllegalArgumentException.h IllegalArgumentException.cpp
	SingularValueDecomposition.h SingularValueDecomposition.cpp
	Matrix.h Matrix.cpp
	FFT.h FFT.cpp
	MFCC.h MFCC.cpp
	WavFileException.h WavFileException.cpp
	WavFile.h WavFile.cpp
	AudioSampleReader.h AudioSampleReader.cpp
	SoundProcessingTask3.h SoundProcessingTask3.cpp
	main_SoundProcessingTask3.cpp
)