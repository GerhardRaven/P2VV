DEPDIR = .deps
df = $(DEPDIR)/$(*F)

TOPDIR = ..
SRCDIR = $(TOPDIR)/src
DICTDIR = $(TOPDIR)/dict
INCDIR = $(TOPDIR)/include

CPP = g++
LD  = g++
ROOTCONFIG = root-config
CPPFLAGS := $(shell $(ROOTCONFIG) --cflags) -I$(INCDIR) -Wall -O2 -pipe -ggdb
LDFLAGS := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit \
	-lRooFitCore -lMathCore -lMathMore

SOURCES =				\
	Moments.cxx			\
	P2VV.cxx			\
	RooBTagDecay.cxx		\
	RooGammaPdf.cxx			\
	RooMultiCatGenerator.cxx	\
	RooMultiMultinomial.cxx		\
	RooP2VVAngleBasis.cxx		\
	RooEffHistProd.cxx		\
	RooThresholdPdf.cxx

OBJECTS = $(SOURCES:%.cxx=%.o)

.PHONY: all clean

all: .deps libP2VV.so

%.o : $(SRCDIR)/%.cxx
	$(CPP) $(CPPFLAGS) -fPIC -DPIC -MMD -c $<
	@cp $*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

%.o : %.cxx
	$(CPP) $(CPPFLAGS) -fPIC -DPIC -MMD -c $<
	@cp $*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

P2VV_dict.cxx: $(wildcard $(INCDIR)/*.h) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h
	rootcint -f $@ -c -I$(INCDIR) $(DICTDIR)/P2VV_dict.h $(DICTDIR)/P2VV_LinkDef.h

libP2VV.so: $(OBJECTS) P2VV_dict.o
	$(LD) $(LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf libP2VV.so $(OBJECTS) P2VV_dict.* *.bak *.aux

$(DEPDIR):
	@mkdir -p $@
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
