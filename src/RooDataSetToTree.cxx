/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * Authors:                                                                  *
 *   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           *
 *                                                                           *
 * Copyright (c) 2012, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include <memory>
#include <set>
#include <vector>

#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TTree.h"

TTree* RooDataSetToTree(const RooDataSet& dataSet, const char* name,
      const char* title, const char* branchList, Bool_t RooFitFormat)
{
  // get branch names
  std::set<TString> branches;
  TString branchStr(branchList);
  while (branchStr.Length() > 0) {
    Ssiz_t sep = branchStr.First(',');
    if (sep < 0) sep = branchStr.Length();

    TString name(((TString)branchStr(0, sep)).Strip(TString::kBoth));
    if (name.Length() > 0) branches.insert(name);
    branchStr = branchStr(sep + 1, branchStr.Length() - sep - 1);
  }

  // get tree name and title
  TString treeName(name);
  TString treeTitle(title);
  if (treeName.Length() < 1) treeName = dataSet.GetName();
  if (treeTitle.Length() < 1) treeTitle = dataSet.GetTitle();

  // build tree
  TTree* tree = new TTree(treeName, treeTitle);
  std::vector<RooRealVar*>  realVars;
  std::vector<RooCategory*> categories;
  std::vector<Double_t*>    realVarAdds;
  std::vector<char*>        catLabelAdds;
  std::vector<Int_t*>       catIndexAdds;

  std::auto_ptr<TIterator> varIter(dataSet.get(0)->createIterator());
  RooAbsArg* arg = 0;
  while ((arg = static_cast<RooAbsArg*>(varIter->Next()))) {
    TString name(arg->GetName());
    if (branches.size() > 0 && branches.count(name) == 0) continue;

    RooRealVar*  realVar  = dynamic_cast<RooRealVar*>(arg);
    RooCategory* category = dynamic_cast<RooCategory*>(arg);
    if (realVar != 0) {
      // create RooRealVar branch
      Double_t* var = new Double_t(0.);
      tree->Branch(name, var, name + "/D");

      realVars.push_back(realVar);
      realVarAdds.push_back(var);
    } else if (category != 0) {
      // create RooCategory branch
      char* catLabel = new char[64];
      Int_t* catInd  = new Int_t(0);
      if (RooFitFormat) {
        tree->Branch(name + "_lbl", catLabel, name + "_lbl/C");
        tree->Branch(name + "_idx", catInd, name + "_idx/I");
        catIndexAdds.push_back(catInd);
        catLabelAdds.push_back(catLabel);
      } else {
        tree->Branch(name, catInd, name + "/I");
        catIndexAdds.push_back(catInd);
      }

      categories.push_back(category);
    }
  }

  for (Int_t event = 0; event < dataSet.numEntries(); ++event) {
    dataSet.get(event);

    // set RooRealVar values
    for (Int_t realVarIter = 0; realVarIter < (Int_t)realVars.size();
        ++realVarIter) {
      *realVarAdds.at(realVarIter) = realVars.at(realVarIter)->getVal();
    }

    // set RooCategory values
    for (Int_t catIter = 0; catIter < (Int_t)categories.size(); ++catIter) {
      if (RooFitFormat) {
        char* labelAdd = catLabelAdds.at(catIter);
        TString label(categories.at(catIter)->getLabel());
        Ssiz_t labelLength = label.Length();
        if (labelLength > 63) labelLength = 63;
        for (Ssiz_t charIter = 0; charIter < labelLength; ++charIter)
          labelAdd[charIter] = label(charIter);
        labelAdd[labelLength] = 0;
      }

      *catIndexAdds.at(catIter) = categories.at(catIter)->getIndex();
    }

    tree->Fill();
  }

  return tree;
}
