CC = g++
CFLAGS = -Wall -O3 --std=c++11

Variational_Monte_Carlo.exe : Variational_Monte_Carlo.o random.o
	$(CC) random.o Variational_Monte_Carlo.o -o Variational_Monte_Carlo.exe
Variational_Monte_Carlo.o : Variational_Monte_Carlo.cpp
	$(CC) -c Variational_Monte_Carlo.cpp -o Variational_Monte_Carlo.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Variational_Monte_Carlo.exe output.*
