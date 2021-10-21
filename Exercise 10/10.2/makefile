MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
CC = g++
CFLAGS = -Wall -O3 --std=c++11
OBJECTS = main.o random.o TSP.o


main.exe : $(OBJECTS) 
	$(CC) $(MPI_COMPILE_FLAGS) $(OBJECTS) $(MPI_LINK_FLAGS) -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TSP.o : TSP.cpp TSP.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
