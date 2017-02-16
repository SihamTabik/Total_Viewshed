CC=g++
CFLAGS= -openmp -ipo
CFLAGS= -openmp -O3 -no-prec-div -ip -xAVX
CFLAGS= -fopenmp -D_OPENMP -Wno-unused-result -Wno-write-strings -g
CFLAGS= -fopenmp -Ofast -D_OPENMP -Wno-unused-result -Wno-write-strings 
CFLAGS= -fopenmp -O3 -D_OPENMP  -Wno-write-strings 
SECFLAGS=
C++11= -std=c++11
C++11= 
SECFLAGS=-mrecip -funsafe-math-optimizations -ffinite-math-only  -fno-trapping-math
SECFLAGS=-ffast-math 
LFLAGS=-fopenmp -L.  -lboost_system -lboost_filesystem
LFLAGS=-fopenmp
OBJDIR=obj
OBJECTS= $(OBJDIR)/LinkedList.o $(OBJDIR)/Sector.o $(OBJDIR)/AuxFunc.o $(OBJDIR)/main.o $(OBJDIR)/earth.o $(OBJDIR)/lodepng.o 
OBJECTS_P=$(OBJDIR)/LinkedList.o $(OBJDIR)/Sector.o $(OBJDIR)/precomp.o 
OBJECTS_T=$(OBJDIR)/AuxFunc.o $(OBJDIR)/test.o  $(OBJDIR)/earth.o $(OBJDIR)/lodepng.o

all: compute test

compute: $(OBJECTS) Makefile
	$(CC) -o comp.out $(OBJECTS)  $(LFLAGS)

precompute: $(OBJECTS_P) Makefile
	$(CC) -o precomp.out $(OBJECTS_P) $(LFLAGS)
	
test: $(OBJECTS_T) Makefile
	$(CC) -o test.out $(OBJECTS_T) $(LFLAGS)
	
$(OBJDIR)/Sector.o: Sector.cpp Makefile Defs.h Closest.h Sector.h
	$(CC) -c $(CFLAGS)  $(SECFLAGS) -o $@ $<
	
$(OBJDIR)/precomp.o: precomp.cpp Makefile Defs.h AuxFunc.h
	$(CC) -c $(CFLAGS) -o $@ $< 
	
$(OBJDIR)/main.o: main.cpp Makefile  Defs.h AuxFunc.h
	$(CC) -c $(CFLAGS) $(C++11) -o $@ $<

$(OBJDIR)/AuxFunc.o: AuxFunc.cpp Makefile Defs.h AuxFunc.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/LinkedList.o: LinkedList.cpp Makefile Defs.h
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o: %.c Makefile Defs.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: %.cpp Makefile Defs.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:	cleanx cleanp

cleanx:
	rm -rf $(OBJDIR)/*.o *.out

cleanp:
	rm -rf *.out


