DEPDIR = .deps
df = $(DEPDIR)/$(*F)

CPP = g++
LD  = g++
ROOTCONFIG = root-config
CPPFLAGS := $(shell $(ROOTCONFIG) --cflags) -Wall -O2 -pipe -ggdb
LDFLAGS := $(shell $(ROOTCONFIG) --libs) -lRooFit -lFoam -lMinuit \
	-lRooFitCore -lMathCore -lMathMore

SOURCES =				\
    RooAddition_.cxx		\
    RooLegendre.cxx			\
    RooP2VVAngleBasis.cxx		\
    RooSpHarmonic.cxx		\
    RooMultiCatGenerator.cxx		\
    RooBTagDecay.cxx		\
    RooThresholdPdf.cxx		\
    RooGammaPdf.cxx		\
    p2vv_dict.cxx

OBJECTS = $(SOURCES:%.cxx=%.o)

.PHONY: all clean

all: .deps libp2vv.so

%.o : %.cxx
	$(CPP) $(CPPFLAGS) -fPIC -DPIC -MMD -c $<
	@cp $*.d $(df).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	rm -f $*.d

p2vv_dict.cxx: $(wildcard *.h)
	rootcint -f p2vv_dict.cxx -c p2vv.h p2vv_LinkDef.h

libp2vv.so: $(OBJECTS)
	$(LD) $(LDFLAGS) -shared -o $@ $^

clean:
	-rm -rf libp2vv.so $(OBJECTS) p2vv_dict.* *.pyc *.bak *.aux

$(DEPDIR):
	mkdir $@
-include $(SOURCES:%.cxx=$(DEPDIR)/%.P)
