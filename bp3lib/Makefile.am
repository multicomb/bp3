INCLUDES = $(CCP4_CXXFLAGS) $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS)

#libtool
lib_LTLIBRARIES = libbp3.la
libbp3_la_SOURCES = \
crystal.C gauss.C input.C likelihood.C minimizer.C misc.C model.C pavol.C sfdata.C

library_include_HEADERS = \
crystal.h gauss.h input.h likelihood.h minimizer.h misc.h model.h sfdata.h tabfunc.h

library_includedir = $(pkgincludedir)/bp3

libbp3_a_AR = $(AR) $(AR_FLAGS)

libbp3_la_LDFLAGS = -version-info 1:1
