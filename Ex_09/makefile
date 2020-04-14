CFLAGS = -Wall -O3 --std=c++11
main.exe: main.o funzioni.o
	g++ main.o funzioni.o -o main.exe
main.o: main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)
funzioni.o: funzioni.cpp funzioni.h
	g++ -c funzioni.cpp -o funzioni.o $(CFLAGS)
