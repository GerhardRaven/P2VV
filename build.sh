#!/bin/bash

source setP2VVEnv.sh

LIBNAME=libP2VV
DICTNAME=P2VVDict

SRCDIR=src
INCDIR=P2VV
DICTDIR=dict

CPP=g++
LD=g++
ROOTCONFIG=root-config

CPPFLAGS="$($ROOTCONFIG --cflags)"" -I$P2VVROOT/$INCDIR -Wall -O2 -pipe -ggdb"
LDFLAGS="$($ROOTCONFIG --libs)"" -lRooFit -lFoam -lMinuit -lRooFitCore\
    -lMathCore -lMathMore"

mkdir -p build
cd build
rm -f *.so *.o *.d *.h *.cxx

# compile code
for SRCFILE in\
    Moments.cxx\
    RooP2VVAngleBasis.cxx\
    RooMultiCatGenerator.cxx\
    RooBTagDecay.cxx\
    RooThresholdPdf.cxx\
    RooGammaPdf.cxx
do
  $CPP $CPPFLAGS -fPIC -DPIC -MMD -c $P2VVROOT/$SRCDIR/$SRCFILE
done

# generate dictionary source
rootcint -f $DICTNAME.cxx -c -I$P2VVROOT/$INCDIR\
    $P2VVROOT/$DICTDIR/P2VV.h $P2VVROOT/$DICTDIR/P2VVLinkDef.h

# compile dictionary
$CPP $CPPFLAGS -fPIC -DPIC -MMD -c $DICTNAME.cxx

# build shared library
$LD $LDFLAGS -shared -o $LIBNAME.so *.o

mkdir -p ../lib
rm -f ../lib/$LIBNAME.so
cp -f $LIBNAME.so ../lib/

cd ..

