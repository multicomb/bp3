CXX = g++
CC  = gcc
LD  = g++
F90 = gfortran

CXX = $(GCC47)/bin/g++ 
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
#-Wstrict-aliasing=2 $(OMPFLAGS)

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
CXXFLAGS += -I./bp3lib -I$(CCP4)/include -I$(CCP4)/include/ccp4 
LDFLAGS  =  -L./bp3lib -lbp3  -g -pg
LDFLAGS += -L$(CCP4)/lib -lccp4f -lccp4c  -lclipper-ccp4 -lclipper-contrib -lclipper-minimol -lclipper-mmdb -lclipper-core -lccif -lmmdb  -lrfftw -lfftw
LDFLAGS += -llapack -lblas 

SRCPATH = ./
SRC = bp3main.C functions.C bp3likelihood.C bp3input.C
OBJ = $(SRC:%.C=%.o) 

PROG = bp3

RM = /bin/rm

all:  $(PROG) 


$(PROG): $(OBJ) 
	$(LD) $^ -o $@ $(OMPFLAGS) 	$(LDFLAGS) 

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@



clean:
	/bin/rm -rf *.o $(PROG) 

clean_all:
	/bin/rm -rf *.o $(PROG) *~

$(OBJ): bp3input.h bp3likelihood.h 
$(PROG): bp3lib/libbp3.a

