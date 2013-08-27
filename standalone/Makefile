DEPDIR = .deps
SRCDIR = src
INCDIR = include
BUILDDIR = .build
LIBDIR = lib
DICTDIR = dict

df = $(DEPDIR)/$(*F)

LD  = $(CXX)
CYTHON = cython
ROOTCONFIG = root-config
PYTHONCONFIG = python-config
COMMONFLAGS = -Wall -O2 -pipe -ggdb -fPIC -DPIC
CXXFLAGS := $(COMMONFLAGS) $(shell $(ROOTCONFIG) --cflags) -I$(INCDIR) -I. 
LDFLAGS  := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit -lRooFitCore -lMathCore -lMathMore
## Reflex dictionary testing flags:
## LDFLAGS  := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit -lRooFitCore -lMathCore -lMathMore -lReflex -lCintex

HAVE_WORKING_PYTHONCONFIG := $(shell python-config  --prefix 2>&1 1> /dev/null ; if [[ "$?" -eq  "0" ]]  ; then echo 'yes' ; else echo 'no' ; fi )
ifeq ($(HAVE_WORKING_PYTHONCONFIG), yes)
PYTHON_CFLAGS := $(shell $(PYTHONCONFIG) --cflags | sed -e "s=/afs/cern.ch/sw=/cvmfs/lhcb.cern.ch/lib=")
PYTHON_LDFLAGS := $(shell $(PYTHONCONFIG) --ldflags | sed -e "s=/afs/cern.ch/sw=/cvmfs/lhcb.cern.ch/lib=")
else
PYTHON = $(shell which python)
PYTHON_VERSION := $(shell $(PYTHON) --version 2>&1 | grep -o "2\.[567]")
PYTHON_DIR := $(shell echo `dirname $(PYTHON)`/..)
PYTHON_CFLAGS := -I$(PYTHON_DIR)/include/python$(PYTHON_VERSION)
PYTHON_LDFLAGS := -L$(PYTHON_DIR)/lib -lpython$(PYTHON_VERSION)
endif



CFLAGS := $(COMMONFLAGS)
PROF_LDFLAGS := $(PYTHON_LDFLAGS) -lprofiler

SOURCES = $(wildcard $(SRCDIR)/*.cxx)

BINDINGSDIR = src/bindings
DILUTION_SOURCES = $(wildcard $(BINDINGSDIR)/dilution*.c)

OBJECTS = $(SOURCES:$(SRCDIR)/%.cxx=$(BUILDDIR)/%.o) $(BUILDDIR)/P2VV_dict.o
DILUTION_OBJECTS = $(DILUTION_SOURCES:$(BINDINGSDIR)/%.c=$(BUILDDIR)/%.o)

vpath %.cxx $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.c   $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.h   $(INCDIR)/P2VV:$(DICTDIR):$(BUILDDIR)
vpath %.o   $(BUILDDIR)

.PHONY: all clean profiler

all: $(DEPDIR) $(LIBDIR) $(BUILDDIR) .deps $(LIBDIR)/libP2VV.so $(LIBDIR)/dilution.so

profiler: $(LIBDIR)/profiler.so

$(BUILDDIR)/P2VV_dict.o: $(BUILDDIR)/P2VV_dict.cxx

$(BUILDDIR)/%.o : %.cxx %.h
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/%.o : bindings/%.c
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/P2VV_dict.cxx: $(wildcard $(INCDIR)/P2VV/*.h) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h
	rootcint -f $@ -c -I$(INCDIR) -I$(DICTDIR) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h
##	genreflex $(DICTDIR)/P2VV_dict.h -o $@ -s $(DICTDIR)/P2VV_dict.xml --gccxmlopt="-I${ROOTSYS}/include -I$(INCDIR)"

$(LIBDIR)/libP2VV.so: $(OBJECTS) $(BUILDDIR)/P2VV_dict.o
	$(LD) $(LDFLAGS) -shared -o $@ $^

$(BINDINGSDIR)/profiler.c: $(SRCDIR)/profiler.pyx
	$(CYTHON) -o $@ $^

$(LIBDIR)/profiler.so: $(BUILDDIR)/profiler.o
	$(LD) $(PROF_LDFLAGS) -shared -o $@ $^

$(LIBDIR)/dilution.so: $(DILUTION_OBJECTS)
	$(LD) $(PYTHON_LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf $(LIBDIR)/*.so $(DILUTION_OBJECTS) $(OBJECTS) P2VV_dict.* *.pyc *.bak *.aux $(BUILDDIR)/* texput.log

$(DEPDIR) $(LIBDIR) $(BUILDDIR):
	@if [ ! -d $@ ]; then mkdir $@; fi
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
