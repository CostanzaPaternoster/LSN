CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn.exe : MolDyn.o random.o
	$(CC) random.o MolDyn.o -o MolDyn.exe
MolDyn.o : MolDyn.cpp
	$(CC) -c MolDyn.cpp -o MolDyn.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o MolDyn.exe seed.out
