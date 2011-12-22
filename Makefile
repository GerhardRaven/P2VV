DEPDIR = .deps
SRCDIR = src
INCDIR = include
BUILDDIR = build
LIBDIR = lib
DICTDIR = dict
df = $(DEPDIR)/$(*F)

#CXX = g++
LD  = $(CXX)
ROOTCONFIG = root-config
CXXFLAGS := $(shell $(ROOTCONFIG) --cflags) -Wall -O2 -pipe -ggdb -I$(INCDIR) -I.
LDFLAGS := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit -lRooFitCore -lMathCore -lMathMore

SOURCES = $(wildcard $(SRCDIR)/*.cxx)

OBJECTS = $(SOURCES:$(SRCDIR)/%.cxx=$(BUILDDIR)/%.o) $(BUILDDIR)/P2VV_dict.o

vpath %.cxx $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.h   $(INCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.o   $(BUILDDIR)

.PHONY: all clean

all: $(DEPDIR) $(LIBDIR) $(BUILDDIR) .deps $(LIBDIR)/libP2VV.so 

$(BUILDDIR)/P2VV_dict.o: $(BUILDDIR)/P2VV_dict.cxx

$(BUILDDIR)/%.o : %.cxx %.h
	$(CXX) $(CXXFLAGS) -fPIC -DPIC -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/P2VV_dict.cxx: $(wildcard $(INCDIR)/*.h) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h
	rootcint -f $@ -c -I$(INCDIR) -I$(DICTDIR) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h

$(LIBDIR)/libP2VV.so: $(OBJECTS) $(BUILDDIR)/P2VV_dict.o
	$(LD) $(LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf $(LIBDIR)/libP2VV.so $(OBJECTS) P2VV_dict.* *.pyc *.bak *.aux $(BUILDDIR)/* texput.log

$(DEPDIR) $(LIBDIR) $(BUILDDIR):
	@if [ ! -d $@ ]; then mkdir $@; fi
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
