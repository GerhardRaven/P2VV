DEPDIR = .deps
df = $(DEPDIR)/$(*F)

TOPDIR = ..
SRCDIR = $(TOPDIR)/src
DICTDIR = $(TOPDIR)/dict
INCDIR = $(TOPDIR)/P2VV

CPP = g++
LD  = g++
ROOTCONFIG = root-config
CPPFLAGS := $(shell $(ROOTCONFIG) --cflags) -Wall -O2 -pipe -ggdb -I$(INCDIR)
LDFLAGS := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit \
	-lRooFitCore -lMathCore -lMathMore

SOURCES =				\
	P2VV_dict.cxx			\
	Moments.cxx			\
	P2VV.cxx			\
	RooBTagDecay.cxx		\
	RooGammaPdf.cxx			\
	RooMultiCatGenerator.cxx	\
	RooMultiMultinomial.cxx		\
	RooP2VVAngleBasis.cxx		\
	RooThresholdPdf.cxx

OBJECTS = $(SOURCES:%.cxx=%.o)

.PHONY: all clean

all: .deps libp2vv.so

%.o : %.cxx
	$(CPP) $(CPPFLAGS) -fPIC -DPIC -MMD -c $<
	@cp $*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

P2VV_dict.cxx: $(wildcard $(INCDIR)/*.h) $(DICTDIR)/P2VVDict.h $(DICTDIR)/P2VV_LinkDef.h
	rootcint -f $@ -c $(DICTDIR)/P2VVDict.h $(DICTDIR)/P2VV_LinkDef.h

libp2vv.so: $(OBJECTS)
	$(LD) $(LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf libP2VV.so $(OBJECTS) P2VV_dict.* *.bak *.aux

$(DEPDIR):
	mkdir $@
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
