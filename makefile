LIBNAME = PlotProcessedData

#Necessary to use shell built-in commands
SHELL=bash

USERINCLUDES += -Iinclude -I$(RATROOT)/include -I$(ROOTSYS)/include -I$(RATROOT)/src/stlplusINCLUDE
USERLIBS += $(shell root-config --cflags) -L$(RATROOT)/lib -lRATEvent_$(RATSYSTEM) -L$(ROOTSYS)/lib $(ROOTLIBS) 

CXXFLAGS = -Wall -W -O2 -std=c++0x -Wno-deprecated-declarations -Wno-unused-parameter
LDFLAGS = -shared -Wall -W 

CXX=g++
LD=g++
RM = rm -f

CXXFLAGS += $(USERINCLUDES)
LIBS += $(USERLIBS)

# A list of directories
BASEDIR = $(shell pwd)
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

# Build a list of srcs and bins to build
SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))
BINS=$(subst $(TESTDIR), $(EXEDIR),$(subst .$(TEST_EXT),,$(EXES)))

all: 
	@$(MAKE) --no-print-directory lib $(BINS)

docs: all
	doxygen Doxyfile
    
$(EXEDIR)/%:  $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/src/*.cc) 
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) 


$(OBJDIR)/%.$(OBJ_EXT):  $(SRCDIR)/%.cc $(BASEDIR)/interface/%.hh
	$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so:  $(OBJS) $(LIBDEPENDS) 
	$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so $(OBJS) $(LIBS)
    
lib: $(LIBDIR)/lib$(LIBNAME).so

clean:
	rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)
