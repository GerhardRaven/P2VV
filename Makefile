DEPDIR = .deps
SRCDIR = src
INCDIR = include
BUILDDIR = build
LIBDIR = lib
DICTDIR = dict
df = $(DEPDIR)/$(*F)

CPP = g++
LD  = g++
ROOTCONFIG = root-config
CPPFLAGS := $(shell $(ROOTCONFIG) --cflags) -Wall -O2 -pipe -ggdb -I$(INCDIR) -I$(DICTDIR)
LDFLAGS := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit \
	-lRooFitCore -lMathCore -lMathMore

SOURCES =				\
    RooP2VVAngleBasis.cxx		\
    RooMultiCatGenerator.cxx		\
    RooBTagDecay.cxx		\
    RooThresholdPdf.cxx		\
    RooGammaPdf.cxx		\
    RooEffHistProd.cxx          

OBJECTS = $(SOURCES:%.cxx=$(BUILDDIR)/%.o) $(BUILDDIR)/p2vv_dict.o

vpath %.cxx $(SRCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.h   $(INCDIR):$(DICTDIR):$(BUILDDIR)
vpath %.o   $(BUILDDIR)

.PHONY: all clean

all: $(DEPDIR) $(LIBDIR) $(BUILDDIR) .deps $(LIBDIR)/libp2vv.so 

$(BUILDDIR)/p2vv_dict.o: $(BUILDDIR)/p2vv_dict.cxx $(BUILDDIR)/p2vv_dict.h

$(BUILDDIR)/%.o : %.cxx %.h
	$(CPP) $(CPPFLAGS) -fPIC -DPIC -MMD -c $< -o $@
	@cp $(BUILDDIR)/$*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$*.d >> $(df).P; \
	rm -f $(BUILDDIR)/$*.d

$(BUILDDIR)/p2vv_dict.cxx: $(wildcard $(INCDIR)/*.h) $(INCDIR)/p2vv.h $(DICTDIR)/p2vv_Linkdef.h
	rootcint -f $@ -c -I$(INCDIR) -I$(DICTDIR) p2vv.h $(DICTDIR)/p2vv_LinkDef.h

$(LIBDIR)/libp2vv.so: $(OBJECTS) $(BUILDDIR)/p2vv_dict.o
	$(LD) $(LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf libp2vv.so $(OBJECTS) p2vv_dict.* *.pyc *.bak *.aux $(BUILDDIR)/*

$(DEPDIR) $(LIBDIR) $(BUILDDIR):
	mkdir $@
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
