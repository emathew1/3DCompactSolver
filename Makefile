CC=g++

# CFLAGS will be the options passed to the compiler. 
CFLAGS= -O3 -std=c++11 
#CFLAGS= -g -std=c++11 -fopenmp 
OBJECTS  = 3DCompact.o Utils.o CSolver.o IdealGas.o SpongeBC.o TimeStepping.o BC.o Macros.o

all: 3DCompact

2D_HOCFD:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ 

3DCompact.o: 3DCompact.cpp Utils.hpp BC.hpp TimeStepping.hpp CSolver.hpp
	$(CC) $(CFLAGS) -c $< 

Utils.o: Utils.cpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

CSolver.o: CSolver.hpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp IdealGas.hpp SpongeBC.hpp
	$(CC) $(CFLAGS) -c $<

IdealGas.o: Macros.hpp Domain.hpp IdealGas.cpp 
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf   *.o 
