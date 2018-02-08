GXX = g++ -std=c++17 -Wall -Werror -pedantic
#GXX = g++ -m64 -std=c++17 -Wall -Werror -pedantic

##################################################

all: main_Complex.exe main_DTW.exe main_WavFile.exe main_SoundProcessingTask3.exe

##################################################

main_Complex.exe: Complex.o main_Complex.o
	$(GXX) Complex.o main_Complex.o -o main_Complex.exe

main_DTW.exe: DTW.o main_DTW.o
	$(GXX) DTW.o main_DTW.o -o main_DTW.exe

main_WavFile.exe: WavFileException.o WavFile.o main_WavFile.o
	$(GXX) WavFileException.o WavFile.o main_WavFile.o -o main_WavFile.exe

main_SoundProcessingTask3.exe: DTW.o IllegalArgumentException.o SingularValueDecomposition.o Matrix.o FFT.o MFCC.o WavFileException.o WavFile.o AudioSampleReader.o SoundProcessingTask3.o main_SoundProcessingTask3.o
	$(GXX) DTW.o IllegalArgumentException.o SingularValueDecomposition.o Matrix.o FFT.o MFCC.o WavFileException.o WavFile.o AudioSampleReader.o SoundProcessingTask3.o main_SoundProcessingTask3.o -o main_SoundProcessingTask3.exe

##################################################

AudioSampleReader.o:
	$(GXX) -c AudioSampleReader.cpp

Complex.o:
	$(GXX) -c Complex.cpp

DTW.o:
	$(GXX) -c DTW.cpp

FFT.o:
	$(GXX) -c FFT.cpp

IllegalArgumentException.o:
	$(GXX) -c IllegalArgumentException.cpp
	
Matrix.o:
	$(GXX) -c Matrix.cpp

MFCC.o:
	$(GXX) -c MFCC.cpp

main_Complex.o:
	$(GXX) -c main_Complex.cpp

main_DTW.o:
	$(GXX) -c main_DTW.cpp

main_SoundProcessingTask3.o:
	$(GXX) -c main_SoundProcessingTask3.cpp

main_WavFile.o:
	$(GXX) -c main_WavFile.cpp

SingularValueDecomposition.o:
	$(GXX) -c SingularValueDecomposition.cpp
	
SoundProcessingTask3.o:
	$(GXX) -c SoundProcessingTask3.cpp
	
WavFile.o:
	$(GXX) -c WavFile.cpp

WavFileException.o:
	$(GXX) -c WavFileException.cpp