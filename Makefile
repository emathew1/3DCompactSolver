
include Makefile.in

OBJECTS  = 3DCompact.o Utils.o CSolver.o CSolver_AWS.o Derivatives.o Filter.o AeroOptics.o
POSTPROOBJ = Utils.o Derivatives.o PostProcess.o VisitWriter.o

all: 3D_HOCFD POST_HOCFD

3DCompact.o: 3DCompact.cpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp CSolver.hpp CSolver_AWS.hpp AbstractCSolver.hpp AbstractRK.hpp RK4.hpp TVDRK3.hpp
	$(CC) $(CFLAGS) -c $< 

Utils.o: Utils.cpp Utils.hpp Domain.hpp
	$(CC) $(CFLAGS) -c $<

CSolver.o: CSolver.cpp CSolver.hpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp IdealGas.hpp SpongeBC.hpp Derivatives.hpp Filter.hpp AbstractCSolver.hpp
	$(CC) $(CFLAGS) -c $<

CSolver_AWS.o: CSolver_AWS.cpp CSolver_AWS.hpp Macros.hpp Utils.hpp BC.hpp TimeStepping.hpp IdealGas.hpp SpongeBC.hpp Derivatives.hpp Filter.hpp AbstractCSolver.hpp PngWriter.hpp
	$(CC) $(CFLAGS) -c $<

Derivatives.o: Derivatives.cpp Derivatives.hpp Macros.hpp Utils.hpp Domain.hpp BC.hpp
	$(CC) $(CFLAGS) -c $<

Filter.o: Filter.cpp Filter.hpp Derivatives.hpp BC.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

AeroOptics.o: AeroOptics.cpp AeroOptics.hpp Macros.hpp Utils.hpp Domain.hpp IdealGas.hpp
	$(CC) $(CFLAGS) -c $<

VisitWriter.o: VisitWriter.cpp VisitWriter.hpp
	$(CC) $(CFLAGS) -c $<

PostProcess.o:	PostProcess.cpp Macros.hpp Utils.hpp Domain.hpp Derivatives.hpp BC.hpp IdealGas.hpp VisitWriter.hpp
	$(CC) $(CFLAGS) -c $<

3D_HOCFD:  $(OBJECTS)
	$(CC) $(CFLAGS) -I$(INC) $(OBJECTS) -o $@ -L$(LIB) $(LIBF)

POST_HOCFD: $(POSTPROOBJ)
	$(CC) $(CFLAGS) -I$(INC) $(POSTPROOBJ) -o $@ -L$(LIB) $(LIBF) 

clean: 
	rm -rf   *.o 


