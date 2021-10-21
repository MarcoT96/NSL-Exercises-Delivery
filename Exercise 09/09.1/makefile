CC = g++
CFLAGS = -Wall -O3 --std=c++11
OBJECTS = main.o random.o GA_TSP.o

main.exe : $(OBJECTS) 
	$(CC) $(OBJECTS) -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
GA_TSP.o : GA_TSP.cpp GA_TSP.h
	$(CC) -c GA_TSP.cpp -o GA_TSP.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
