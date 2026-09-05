SUBDIRS = src utilities

AM_LDFLAGS += -L$(top_builddir)/src
AM_CXXFLAGS += -I$(top_srcdir)/src

bin_PROGRAMS = HadronsMILC

HadronsMILC_SOURCES = main.cpp
# Name the library as a file, not as -lHadronsMILC: a -l flag creates no make
# dependency, so the binary was never relinked when a module changed and
# silently kept the previously compiled module code.
HadronsMILC_LDADD = $(top_builddir)/src/libHadronsMILC.a
