EXEC   = simple_particle_projection

OPTIMIZE =  -O2 -DATHENA4 -Wno-write-strings 

OBJS   = main.o read_athena_header.o load_tracers.o

#CXX    = mpiicpc
#CC     = mpiicpc
CXX    = mpicxx
CC     = mpicxx

INCL   = read_athena_header.hpp load_tracers.hpp

#LIBS   = -lm -lgsl -lgslcblas -lmpich
LIBS   = -lm -lgsl -lgslcblas -lmpi  -stdlib=libstdc++
#LIBS   = -lm -lgsl -lgslcblas -lmpi++

CFLAGS = $(OPTIMIZE) -stdlib=libstdc++
CXXFLAGS = $(OPTIMIZE) -stdlib=libstdc++

$(EXEC): $(OBJS) 
	 $(CXX) $(OBJS) $(LIBS) -o $(EXEC)   
         

$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

