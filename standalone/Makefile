DEPDIR = .deps
SRCDIR = src
INCDIR = include
BUILDDIR = .build
LIBDIR = lib
DICTDIR = dict
df = $(DEPDIR)/$(*F)

#CXX = g++
LD  = $(CXX)
CYTHON = cython
PYTHON = python
ROOTCONFIG = root-config
PKGCONFIG = pkg-config
COMMONFLAGS = -Wall -O2 -pipe -ggdb -fPIC -DPIC
PYTHONVERSION := $(shell $(PYTHON) --version 2>&1 | grep -o "2\.[567]")
CXXFLAGS := $(COMMONFLAGS) $(shell $(ROOTCONFIG) --cflags) -I$(INCDIR) -I.
LDFLAGS  := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit -lRooFitCore -lMathCore -lMathMore

HAVE_PYTHON := $(shell if `$(PKGCONFIG) --exists python-$(PYTHONVERSION)`; then echo "yes"; else echo "no"; fi)
ifeq ($(HAVE_PYTHON), yes)
PYTHON_CFLAGS := $(shell $(PKGCONFIG) --cflags python-$(PYTHONVERSION))
else
PYTHON_CFLAGS = ""
endif

CFLAGS   := $(COMMONFLAGS) $(PYTHON_CFLAGS)
PROF_LDFLAGS := -lpython$(PYTHONVERSION) -lprofiler

SOURCES = $(wildcard $(SRCDIR)/*.cxx)

OBJECTS = $(SOURCES:$(SRCDIR)/%.cxx=$(BUILDDIR)/%.o) $(BUILDDIR)/P2VV_dict.o

vpath %.cxx $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.c   $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.h   $(INCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.o   $(BUILDDIR)

.PHONY: all clean profiler

all: $(DEPDIR) $(LIBDIR) $(BUILDDIR) .deps $(LIBDIR)/libP2VV.so 

profiler: $(LIBDIR)/profiler.so

$(BUILDDIR)/P2VV_dict.o: $(BUILDDIR)/P2VV_dict.cxx

$(BUILDDIR)/%.o : %.cxx %.h
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/%.o : %.c
	$(CC) $(CFLAGS) -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/P2VV_dict.cxx: $(wildcard $(INCDIR)/*.h) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h
	rootcint -f $@ -c -I$(INCDIR) -I$(DICTDIR) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h

$(LIBDIR)/libP2VV.so: $(OBJECTS) $(BUILDDIR)/P2VV_dict.o
	$(LD) $(LDFLAGS) -shared -o $@ $^

$(BUILDDIR)/profiler.c: $(SRCDIR)/profiler.pyx
	$(CYTHON) -o $@ $^

$(LIBDIR)/profiler.so: $(BUILDDIR)/profiler.o
	$(LD) $(PROF_LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf $(LIBDIR)/libP2VV.so $(OBJECTS) P2VV_dict.* *.pyc *.bak *.aux $(BUILDDIR)/* texput.log

$(DEPDIR) $(LIBDIR) $(BUILDDIR):
	@if [ ! -d $@ ]; then mkdir $@; fi
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
