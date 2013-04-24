/*****************************************************************************
 * Project: RooFit                                                           *
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROO_B_DATA_SET_TO_TREE
#define ROO_B_DATA_SET_TO_TREE

#include "Rtypes.h"

class RooDataSet;
class TTree;

TTree* RooDataSetToTree(const RooDataSet& dataSet, const char* name = 0,
      const char* title = 0, const char* branchList = 0,
      Bool_t RooFitFormat = kTRUE);

#endif

