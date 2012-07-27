CXX = g++
# CXX = $(GCC47)/bin/g++

CC  = gcc
F90 = gfortran

LD  = $(CXX)


ifeq ($(OMP), 1)
	OMPFLAGS  = -fopenmp
endif
ifeq ($(OMPGLIBC), 1)
	OMPFLAGS += -D_GLIBCXX_PARALLEL
endif

OFLAGS = -O3 -g 
ifeq ($(GDB), 1)
	OFLAGS = -O0 -g
endif
	
OFLAGS  += -funroll-all-loops
CXXFLAGS = $(OFLAGS) 
# -Wstrict-aliasing=2 $(OMPFLAGS)

ifeq ($(SSE), 1)
	CXXFLAGS += -D__mySSE__ -msse4.1
endif
ifeq ($(AVX), 1)
	CXXFLAGS += -D__myAVX__ -mavx -mtune=corei7-avx -march=corei7-avx
endif
ifeq ($(SSE1), 1)
	CXXFLAGS += -D__mySSE1__ -msse4.1
endif
ifeq ($(RELEASE), 1)
	CXXFLAGS += -DNDEBUG
endif

# OFLAGS += -pg


#CXXFLAGS = -I$(CCP4)/include -I$(CCP4)/include/leiden-software/bp3 -I$(CCP4)/include/ccp4 -I$(CCP4)/lib/cctbx/cctbx_sources/cctbx_project/scitbx/sparse
CXXFLAGS += -I$(CCP4)/include -I$(CCP4)/include/ccp4 

SRCPATH = ./
SRC = crystal.C  gauss.C  input.C  likelihood.C  minimizer.C  misc.C  model.C  pavol.C  sfdata.C
OBJ = $(SRC:%.C=%.o) 

TRGT = libbp3

RM = /bin/rm

all:	  $(TRGT)


$(TRGT): $(OBJ) 
	/bin/rm -f $(TRGT).a
	ar qv $(TRGT).a $^
	ranlib $(TRGT).a

%.o: $(SRCPATH)/%.C
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(TRGT).a

clean_all:
	/bin/rm -rf *.o $(TRGT).a *~

$(OBJ): bp3_version.h  covmat.h  crystal.h  gauss.h  gaussian.h  input.h  likelihood.h  llhood.h  matrix.h  minimizer.h  misc.h  model.h  sfdata.h  tabfunc.h

