#General...
CC=g++


# CFLAGS will be the options passed to the compiler. 
#CFLAGS= -O3 -std=c++11  
CFLAGS= -O3 -std=c++11 -Xpreprocessor -fopenmp -lomp
OBJECTS  = 3DCompact.o Utils.o CSolver.o Derivatives.o Filter.o

all: 3D_HOCFD

3D_HOCFD:  $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $@ 

3DCompact.o: 3DCompact.cpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp CSolver.hpp 
	$(CC) $(CFLAGS) -c $< 

Utils.o: Utils.cpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

CSolver.o: CSolver.cpp CSolver.hpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp IdealGas.hpp SpongeBC.hpp Derivatives.hpp Filter.hpp
	$(CC) $(CFLAGS) -c $<

Derivatives.o: Derivatives.cpp Derivatives.hpp Macros.hpp Utils.hpp Domain.hpp BC.hpp
	$(CC) $(CFLAGS) -c $<

Filter.o: Filter.cpp Filter.hpp Derivatives.hpp BC.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf   *.o 
