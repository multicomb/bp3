bp3_SOURCES = bp3input.C	bp3likelihood.C	bp3main.C	functions.C \
bp3input.h	bp3likelihood.h

bin_PROGRAMS = bp3

INCLUDES = -I../../bp3lib/src $(PTHREAD_CFLAGS) $(CLIPPER_CXXFLAGS) $(CCP4_CXXFLAGS) $(MMDB_CXXFLAGS) $(FFTW_CXXFLAGS)

bp3_LDADD = ../../bp3lib/src/libbp3.la $(PTHREAD_CFLAGS) $(CLIPPER_LIBS) $(CCP4_LIBS) $(FFTW_LIBS) $(PTHREAD_LIBS) $(LAPACK_LIBS)
