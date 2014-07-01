// low level stuff
#include <sys/types.h>

// STD & STL
#include <memory>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <list>
#include <typeinfo>

// ROOT
#include <TH1.h>
#include <TTree.h>
#include <TMatrixT.h>
#include <TEntryList.h>
#include <TDirectory.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TEfficiency.h>

// RooFit
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooLinkedListIter.h>
#include <RooNumber.h>

#include "P2VV/Functions.h"

namespace {
   using std::cout;
   using std::endl;
   using std::vector;
   using std::list;
}

void addSWeightToTree(const std::vector<double>& weights, TTree& tree, const std::string& branch_name)
{
   assert(Long64_t(weights.size()) == tree.GetEntries());

   Double_t w = 0;
   //#       Create the output Branch
   std::string branch_def = branch_name + "/D";
   TBranch* branch = tree.Branch(branch_name.c_str(), &w, branch_def.c_str());

   for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
      w = weights[i];
      branch->Fill();
   }
   tree.FlushBaskets();

}

void addSWeightToTree(double* weights, size_t n, TTree& tree, const std::string& branch_name)
{
   std::vector<double> w(weights, weights + n);
   addSWeightToTree(w, tree, branch_name);
}

void addSWeightToTree(const RooDataSet& ds, TTree& tree, const std::string& branch_name,
                      const std::string& cut)
{
   tree.Draw(">>elist", cut.c_str(), "entrylist");
   TEntryList *cut_list = static_cast<TEntryList*>(gDirectory->Get( "elist" ));

   Long64_t nds(ds.numEntries());
   if (nds != cut_list->GetN()) {
      std::cout << "addSWeightToTree(): ERROR: number of entries in data set does not match number of selected entries in tree"
            << std::endl;
      assert(nds == cut_list->GetN());
   }

   std::vector<Double_t> weights(tree.GetEntries(), 0.);

   Int_t j = 0;
   for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
      if (cut_list->Contains(i)) {
         ds.get(j++);
         weights[i] = ds.weight();
      }
   }
   addSWeightToTree(weights, tree, branch_name);
}

void addIntegerToTree(TTree& tree, Int_t value, const char* branchName) {
  TString branchNameStr(branchName);
  Int_t* output = new Int_t(value);
  TBranch* branch = tree.Branch(branchNameStr, output, branchNameStr + "/I");
  for (Long64_t it = 0; it < tree.GetEntries(); ++it) branch->Fill();
  tree.FlushBaskets();
}

void addFloatToTree(TTree& tree, Double_t value, const char* branchName) {
  TString branchNameStr(branchName);
  Double_t* output = new Double_t(value);
  TBranch* branch = tree.Branch(branchNameStr, output, branchNameStr + "/D");
  for (Long64_t it = 0; it < tree.GetEntries(); ++it) branch->Fill();
  tree.FlushBaskets();
}

void addProductToTree(TTree& tree, vector<TString> inBranches,
    const char* outBranch) {
  vector<Float_t*> inputF;
  vector<Double_t*> inputD;
  for (vector<TString>::const_iterator brIt = inBranches.begin();
      brIt != inBranches.end(); ++brIt) {
    TObjArray* lfList = tree.GetBranch(*brIt)->GetListOfLeaves();
    TString brType = lfList->GetEntries() == 1 ?
        ((TLeaf*)lfList->At(0))->GetTypeName() : "";
    if (brType == "Double_t") {
      Double_t* input = new Double_t(0.);
      tree.SetBranchAddress(*brIt, input);
      inputD.push_back(input);
    } else if (brType == "Float_t") {
      cout << "P2VV - INFO: copyFloatInTree(): values from Float_t branch \""
        << *brIt << "\" will be converted to double precision" << endl;
      Float_t* input = new Float_t(0.);
      tree.SetBranchAddress(*brIt, input);
      inputF.push_back(input);
    } else {
      cout << "P2VV - ERROR: addProductToTree(): branch \"" << *brIt
        << "\" has unknown type \"" << brType << "\"" << endl;
      assert(0);
    }
  }
  if (inputF.size() + inputD.size() <= 0.) {
      cout << "P2VV - ERROR: addProductToTree(): no input branches found"
          << endl;
      assert(0);
  }

  TString branchNameStr(outBranch);
  Double_t* output = new Double_t(0.);
  TBranch* branch = tree.Branch(branchNameStr, output, branchNameStr + "/D");

  for (Long64_t it = 0; it < tree.GetEntries(); ++it) {
    tree.GetEntry(it);
    *output = 1.;
    for (vector<Float_t*>::const_iterator fIt = inputF.begin();
        fIt != inputF.end(); ++fIt) {
      *output *= *(*fIt);
    }
    for (vector<Double_t*>::const_iterator dIt = inputD.begin();
        dIt != inputD.end(); ++dIt) {
      *output *= *(*dIt);
    }
    branch->Fill();
  }
  tree.FlushBaskets();
}

void copyFloatInTree(TTree& tree, const char* inBranch, const char* outBranch) {
  Float_t*  inputF = 0;
  Double_t* inputD = 0;
  TObjArray* lfList = tree.GetBranch(inBranch)->GetListOfLeaves();
  TString brType = lfList->GetEntries() == 1 ?
      ((TLeaf*)lfList->At(0))->GetTypeName() : "";
  if (brType == "Double_t") {
    inputD = new Double_t(0.);
    tree.SetBranchAddress(inBranch, inputD);
  } else if (brType == "Float_t") {
    cout << "P2VV - INFO: copyFloatInTree(): values from Float_t branch \""
      << inBranch << "\" will be converted to double precision" << endl;
    inputF = new Float_t(0.);
    tree.SetBranchAddress(inBranch, inputF);
  } else {
    cout << "P2VV - ERROR: copyFloatInTree(): branch \"" << inBranch
      << "\" has unknown type \"" << brType << "\"" << endl;
    assert(0);
  }

  TString branchNameStr(outBranch);
  Double_t* output = new Double_t(0.);
  TBranch* branch = tree.Branch(branchNameStr, output, branchNameStr + "/D");

  for (Long64_t it = 0; it < tree.GetEntries(); ++it) {
    tree.GetEntry(it);
    *output = inputD != 0 ? *inputD : *inputF;
    branch->Fill();
  }
  tree.FlushBaskets();
}

void addCategoryToTree(TTree& tree, const char* inBranch,
    const char* catBranch, std::vector<Double_t> boundaries,
    std::vector<Int_t> indices) {
  if (indices.size() != boundaries.size() + 1) {
    std::cerr << "ERROR: addCategoryToTree(): number of boundaries ("
        << boundaries.size() << ") and number of indices (" << indices.size()
        << ") do not match" << std::endl;
    assert(indices.size() == boundaries.size() + 1);
  }

  TString branchNameStr(catBranch);
  Int_t* output = new Int_t(0);
  TBranch* branch = tree.Branch(branchNameStr, output, branchNameStr + "/I");

  TObjArray* lfList = tree.GetBranch(inBranch)->GetListOfLeaves();
  TString brType = lfList->GetEntries() == 1 ?
      ((TLeaf*)lfList->At(0))->GetTypeName() : "";
  Double_t* inputD = 0; Float_t* inputF = 0;
  Int_t* inputI = 0; UInt_t* inputUI = 0;
  Long64_t* inputL = 0; ULong64_t* inputUL = 0;
  if (brType == "Double_t") {
    inputD = new Double_t(0.);
    tree.SetBranchAddress(inBranch, inputD);
  } else if (brType == "Float_t") {
    inputF = new Float_t(0.);
    tree.SetBranchAddress(inBranch, inputF);
  } else if (brType == "Int_t") {
    inputI = new Int_t(0);
    tree.SetBranchAddress(inBranch, inputI);
  } else if (brType == "UInt_t") {
    inputUI = new UInt_t(0);
    tree.SetBranchAddress(inBranch, inputUI);
  } else if (brType == "Long64_t") {
    inputL = new Long64_t(0);
    tree.SetBranchAddress(inBranch, inputL);
  } else if (brType == "ULong64_t") {
    inputUL = new ULong64_t(0);
    tree.SetBranchAddress(inBranch, inputUL);
  } else {
    cout << "P2VV - ERROR: addCategoryToTree(): branch \"" << inBranch
      << "\" has unknown type \"" << brType << "\"" << endl;
    assert(0);
  }

  for (Long64_t it = 0; it < tree.GetEntries(); ++it) {
    tree.GetEntry(it);
    Int_t pos(0);
    for (std::vector<Double_t>::const_iterator boundIt = boundaries.begin();
        boundIt != boundaries.end(); ++boundIt) {
      if (brType == "Double_t" && *inputD < *boundIt) break;
      else if (brType == "Float_t" && (Double_t)*inputF < *boundIt) break;
      else if (brType == "Int_t" && (Double_t)*inputI < *boundIt) break;
      else if (brType == "UInt_t" && (Double_t)*inputUI < *boundIt) break;
      else if (brType == "Long64_t" && (Double_t)*inputL < *boundIt) break;
      else if (brType == "ULong64_t" && (Double_t)*inputUL < *boundIt) break;
      ++pos;
    }
    *output = indices[pos];
    branch->Fill();
  }
  tree.FlushBaskets();
}

void addVertexErrors(TTree* tree, const std::list<RooDataSet*>& dss, const std::string& cut) {

   cout << "Reading tree " << tree->GetName() << " to get vertex errors." << endl;

   tree->Draw(">>vxerr_elist", cut.c_str(), "entrylist");
   TEntryList *cut_list = static_cast<TEntryList*>(gDirectory->Get( "vxerr_elist" ));

   Double_t px = 0., py = 0., pz = 0.;
   Double_t m = 5366.7;
   Double_t p = 0.;
   Float_t cov_sv[3][3];
   Float_t cov_pv[3][3];
   Float_t cov_jpsi[3][3];

   Double_t pvx = 0., pvy = 0., pvz = 0.;

   tree->SetBranchAddress("B_s0_OWNPV_X", &pvx);
   tree->SetBranchAddress("B_s0_OWNPV_Y", &pvy);
   tree->SetBranchAddress("B_s0_OWNPV_Z", &pvz);

   tree->SetBranchAddress("B_s0_PX", &px);
   tree->SetBranchAddress("B_s0_PY", &py);
   tree->SetBranchAddress("B_s0_PZ", &pz);

   tree->SetBranchAddress("B_s0_P", &p);
   
   tree->SetBranchAddress("B_s0_ENDVERTEX_COV_", &cov_sv);
   tree->SetBranchAddress("B_s0_OWNPV_COV_", &cov_pv);
   tree->SetBranchAddress("J_psi_1S_ENDVERTEX_COV_", &cov_jpsi);

   TMatrixT<float> P(3, 1);
   TMatrixT<float> dx(3, 1);
   TMatrixT<float> r(1, 1);
   TMatrixT<float> tmp(3, 1);

   RooRealVar* sv_err = new RooRealVar("sv_err", "sv_err", 0, 0.1);
   RooRealVar* pv_err = new RooRealVar("pv_err", "pv_err", 0, 0.1);
   RooRealVar* psi_err = new RooRealVar("jpsi_vx_err", "jpsi_vx_err", 0, 0.1);

   RooDataSet* ds = new RooDataSet("vertex_errors", "vertex_errors", RooArgSet(*sv_err, *pv_err, *psi_err));
   const RooArgSet* obs = ds->get();
   sv_err = static_cast<RooRealVar*>(obs->find(sv_err->GetName()));
   pv_err = static_cast<RooRealVar*>(obs->find(pv_err->GetName()));
   psi_err = static_cast<RooRealVar*>(obs->find(psi_err->GetName()));

   Long64_t n = tree->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      if (i != 0 && i != n && i % (n / 20) == 0) {
         cout << int(double(i + 20) / double(n) * 100) << "% ";
         cout.flush();
      }
      if (!cut_list->Contains(i)) {
         continue;
      } else {
         tree->GetEntry(i);
      }

      P(0, 0) = px;
      P(1, 0) = py;
      P(2, 0) = pz;
      P *= m / (p * p);
         
      TMatrixT<float> P_T(P);
      P_T.T();

      TMatrixT<float> csv(3, 3, &cov_sv[0][0]);
      TMatrixT<float> cpv(3, 3, &cov_pv[0][0]);
      TMatrixT<float> cjpsi(3, 3, &cov_jpsi[0][0]);

      tmp.Mult(csv, P);
      r.Mult(P_T, tmp);

      // Result is (c * tau)^2, set value in ps
      sv_err->setVal(sqrt(r(0, 0)) / 0.299792458);

      tmp.Mult(cpv, P);
      r.Mult(P_T, tmp);

      // Result is (c * tau)^2, set value in ps
      pv_err->setVal(sqrt(r(0, 0)) / 0.299792458);

      tmp.Mult(cjpsi, P);
      r.Mult(P_T, tmp);

      // Result is (c * tau)^2, set value in ps
      psi_err->setVal(sqrt(r(0, 0)) / 0.299792458);
      ds->fill();
   }
   cout << endl;

   cout << "Adding vertex errors to RooDataSets" << endl;
   for (list<RooDataSet*>::const_iterator it = dss.begin(), end = dss.end(); it != end; ++it) {
      (*it)->merge(ds);
   }
}

void addJpsiDLS(TTree* tree, const std::list<RooDataSet*>& dss, const std::string& cut) {

   cout << "Reading tree " << tree->GetName() << " to get vertex errors." << endl;

   tree->Draw(">>dls_elist", cut.c_str(), "entrylist");
   TEntryList *cut_list = static_cast<TEntryList*>(gDirectory->Get( "dls_elist" ));

   Double_t x = 0., y = 0., z = 0.;
   Double_t pv_x = 0., pv_y = 0., pv_z = 0.;

   Float_t cov_pv[3][3];
   Float_t cov_jpsi[3][3];

   tree->SetBranchAddress("J_psi_1S_OWNPV_X", &pv_x);
   tree->SetBranchAddress("J_psi_1S_OWNPV_Y", &pv_y);
   tree->SetBranchAddress("J_psi_1S_OWNPV_Z", &pv_z);

   tree->SetBranchAddress("J_psi_1S_ENDVERTEX_X", &x);
   tree->SetBranchAddress("J_psi_1S_ENDVERTEX_Y", &y);
   tree->SetBranchAddress("J_psi_1S_ENDVERTEX_Z", &z);

   tree->SetBranchAddress("J_psi_1S_OWNPV_COV_", &cov_pv);
   tree->SetBranchAddress("J_psi_1S_ENDVERTEX_COV_", &cov_jpsi);

   TMatrixT<float> X(3, 1);
   TMatrixT<float> r(1, 1);
   TMatrixT<float> tmp(3, 1);

   RooRealVar* dl = new RooRealVar("jpsi_DL", "jpsi_DL", -1000, 1000);
   RooRealVar* dle = new RooRealVar("jpsi_DLE", "jpsi_DLE", -1000, 1000);
   RooRealVar* dls = new RooRealVar("jpsi_DLS", "jpsi_DLS", -1000, 1000);

   RooDataSet* ds = new RooDataSet("dls", "dls", RooArgSet(*dl, *dle, *dls));
   const RooArgSet* obs = ds->get();
   std::string dln = dl->GetName();
   std::string dlen = dle->GetName();
   std::string dlsn = dls->GetName();
   delete dl;
   delete dls;
   delete dle;
   dl = static_cast<RooRealVar*>(obs->find(dln.c_str()));
   dle = static_cast<RooRealVar*>(obs->find(dlen.c_str()));
   dls = static_cast<RooRealVar*>(obs->find(dlsn.c_str()));

   Long64_t n = tree->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      if (i != 0 && i != n && i % (n / 20) == 0) {
         cout << int(double(i + 20) / double(n) * 100) << "% ";
         cout.flush();
      }
      if (!cut_list->Contains(i)) {
         continue;
      } else {
         tree->GetEntry(i);
      }

      double D = sqrt((x - pv_x) * (x - pv_x) + (y- pv_y) * (y- pv_y) + (z- pv_z) + (z- pv_z));
      dl->setVal(D);

      X(0, 0) = x - pv_x;
      X(1, 0) = y - pv_y;
      X(2, 0) = z - pv_z;
      X *= 1 / D;

      TMatrixT<float> X_T(TMatrixT<float>::kTransposed, X);

      TMatrixT<float> cjpsi(3, 3, &cov_jpsi[0][0]);

      tmp.Mult(cjpsi, X);
      r.Mult(X_T, tmp);

      double dle_jpsi = sqrt(r(0, 0));

      TMatrixT<float> cpv(3, 3, &cov_pv[0][0]);
      tmp.Mult(cpv, X);
      r.Mult(X_T, tmp);

      double dle_pv = sqrt(r(0, 0));

      dle->setVal(sqrt(dle_jpsi * dle_jpsi + dle_pv * dle_pv));
      dls->setVal(dl->getVal() / dle->getVal());
      ds->fill();
   }
   cout << endl;

   cout << "Adding vertex errors to RooDataSets" << endl;
   for (list<RooDataSet*>::const_iterator it = dss.begin(), end = dss.end(); it != end; ++it) {
      (*it)->merge(ds);
   }
}

TTree* RooDataSetToTree(const RooDataSet& dataSet, const char* name,
      const char* title, const char* weightName, const char* branchList,
      Bool_t RooFitFormat)
{
  // get branch names
  std::set<TString> branches;
  TString branchStr(branchList);
  while (branchStr.Length() > 0) {
    Ssiz_t sep = branchStr.First(',');
    if (sep < 0) sep = branchStr.Length();

    TString brName(((TString)branchStr(0, sep)).Strip(TString::kBoth));
    if (brName.Length() > 0) branches.insert(brName);
    branchStr = branchStr(sep + 1, branchStr.Length() - sep - 1);
  }

  // get tree name and title
  TString treeName(name);
  TString treeTitle(title);
  if (treeName.Length() < 1) treeName = dataSet.GetName();
  if (treeTitle.Length() < 1) treeTitle = dataSet.GetTitle();

  // create tree
  TTree* tree = new TTree(treeName, treeTitle);

  // create branch in tree for event weights
  TString weightNameStr(weightName);
  Double_t* weightVarAdd = 0;
  if (weightNameStr.Length() > 0) {
    weightVarAdd = new Double_t(0.);
    tree->Branch(weightNameStr, weightVarAdd, weightNameStr + "/D");
  }

  // create branches in tree for data set columns
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

  // loop over events in data set
  for (Int_t event = 0; event < dataSet.numEntries(); ++event) {
    dataSet.get(event);

    // set weight
    if (weightVarAdd != 0) {
      *weightVarAdd = dataSet.weight();
    }

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

RooDataSet* TreeToRooDataSet(TTree& tree, const RooArgSet& observables,
      const char* name, const char* title, const char* cuts,
      const char* indexName, const char* weightName, RooDataSet* origDataSet)
{
  // get tree name and title
  TString dsName(name);
  TString dsTitle(title);
  if (dsName.Length() < 1) dsName = tree.GetName();
  if (dsTitle.Length() < 1) dsTitle = tree.GetTitle();

  TString selStr(cuts);
  TString indStr(indexName);
  TString wStr(weightName);
  if (selStr.Length() < 1 && indStr.Length() < 1 && origDataSet == 0) {
    // default: import tree with RooDataSet constructor
    if (wStr.Length() > 0) {
      return new RooDataSet(dsName, dsTitle, observables,
          RooFit::Import(tree), RooFit::WeightVar(wStr));
    } else {
      return new RooDataSet(dsName, dsTitle, observables,
          RooFit::Import(tree));
    }
  }

  // check number of entries in tree and original data set
  if (origDataSet != 0
      && tree.GetEntries() < (Long64_t)origDataSet->numEntries()) {
    cout << "P2VV - ERROR: TreeToRooDataSet(): less entries in tree than in original data set"
        << endl;
    return 0;
  }

  // get number of data set entries
  Long64_t numEntr = origDataSet == 0 ?
      tree.GetEntries() : (Long64_t)origDataSet->numEntries();
  if (indStr.Length() > 0 && numEntr > (Long64_t)1.e15) {
    cout << "P2VV - ERROR: TreeToRooDataSet(): number of entries with index variable limited to 10^15"
        << endl;
    return 0;
  }

  // get observables in data sets
  const RooArgSet* origObsSet = origDataSet == 0 ?
      new RooArgSet() : origDataSet->get();
  RooArgSet obsSet(observables);

  // find index variable
  RooRealVar* index(0);
  RooRealVar* origIndex(0);
  if (indStr.Length() > 0) {
    index = dynamic_cast<RooRealVar*>(obsSet.find(indStr));
    if (index == 0 && obsSet.find(indStr) != 0) {
      cout << "P2VV - ERROR: TreeToRooDataSet(): index variable is not a RooRealVar"
          << endl;
      return 0;
    }
    origIndex = dynamic_cast<RooRealVar*>(origObsSet->find(indStr));
    if (origIndex == 0 && origObsSet->find(indStr) != 0) {
      cout << "P2VV - ERROR: TreeToRooDataSet(): original index variable is not a RooRealVar"
          << endl;
      return 0;
    }
  }

  // get/set branch status
  std::vector<Bool_t> brStatusList;
  std::auto_ptr<TIterator> brIt(tree.GetListOfBranches()->MakeIterator());
  TBranch* br(0);
  while ((br = dynamic_cast<TBranch*>(brIt->Next())) != 0) {
    brStatusList.push_back(tree.GetBranchStatus(br->GetName()));
    br->SetStatus(kFALSE);
  }

  // initialize observables
  std::map<TString,Double_t*> doubleMap;
  std::map<TString,Float_t*> floatMap;
  std::map<TString,Int_t*> intMap;
  RooAbsArg* arg = 0;
  RooLinkedListIter obsSetIter = obsSet.iterator();
  while ((arg = static_cast<RooAbsArg*>(obsSetIter.Next())) != 0) {
    // check type of observables
    RooRealVar*  real = dynamic_cast<RooRealVar*>(arg);
    RooCategory* cat  = dynamic_cast<RooCategory*>(arg);
    if (real == 0 && cat == 0) {
      cout << "P2VV - WARNING: TreeToRooDataSet(): variable \""
          << arg->GetName() << "\" is not RooRealVar nor RooCategory" << endl;
      obsSet.remove(*arg);
      continue;
    }

    // check if observable is in original data set
    RooAbsArg* origArg = origObsSet->find(arg->GetName());
    if (origArg != 0) {
      if (real != 0) {
        // check if original observable is RooRealVar with correct range
        RooRealVar* origReal = dynamic_cast<RooRealVar*>(origArg);
        if (origReal == 0) {
          cout << "P2VV - WARNING: TreeToRooDataSet(): variable \""
              << real->GetName() << "\" is not a RooRealVar in original data"
              << endl;
          obsSet.remove(*arg);
        } else if (origReal->getMin() != real->getMin()
            || origReal->getMax() != real->getMax()) {
          cout << "P2VV - WARNING: TreeToRooDataSet(): variable \""
              << real->GetName()
              << "\" has not the same range in original data" << endl;
          obsSet.remove(*arg);
        }
      } else {
        // check if original observable is RooCategory with correct states
        RooCategory* origCat = dynamic_cast<RooCategory*>(origArg);
        if (origCat == 0) {
          cout << "P2VV - WARNING: TreeToRooDataSet(): variable \""
              << cat->GetName() << "\" is not a RooCategory in original data"
              << endl;
          obsSet.remove(*arg);
        } else {
          Bool_t checkStates = kTRUE;
          std::auto_ptr<TIterator> catIt(cat->typeIterator());
          RooCatType* type(0);
          while ((type = dynamic_cast<RooCatType*>(catIt->Next())) != 0) {
            if (!origCat->isValidIndex(type->getVal())) checkStates = kFALSE;
          }
          if (!checkStates) {
            cout << "P2VV - WARNING: TreeToRooDataSet(): variable \""
                << cat->GetName()
                << "\" has not the same range in original data" << endl;
            obsSet.remove(*arg);
          }
        }
      }
      continue;
    }

    // don't try to find branch if this is the index variable
    if (index && arg->GetName() == index->GetName()) continue;

    // set branch status
    UInt_t brFound = 0;
    tree.SetBranchStatus(arg->GetName(), kTRUE, &brFound);
    if (brFound < 1) {
      cout << "P2VV - WARNING: TreeToRooDataSet(): branch \""
          << arg->GetName() << "\" not found in tree" << endl;
      obsSet.remove(*arg);
      continue;
    } else if (brFound > 1) {
      cout << "P2VV - WARNING: TreeToRooDataSet(): " << brFound
          << " branches \"" << arg->GetName()
          << "\" found in tree: not reading this branch" << endl;
      obsSet.remove(*arg);
      continue;
    }

    // set branch address and build selection string
    if (real != 0) {
      // set branch address of floating point branch
      TObjArray* lfList = tree.GetBranch(real->GetName())->GetListOfLeaves();
      TString brType = lfList->GetEntries() == 1 ?
          ((TLeaf*)lfList->At(0))->GetTypeName() : "";
      Int_t setAdd(0);
      if (brType == "Float_t") {
        // Float_t branch
        cout << "P2VV - INFO: TreeToRooDataSet(): values from Float_t branch \""
          << real->GetName() << "\" will be converted to double precision"
          << endl;
        Float_t* brAdd = new Float_t(0.);
        setAdd = tree.SetBranchAddress(real->GetName(), brAdd);
        if (setAdd == 0) floatMap[TString(real->GetName())] = brAdd;
        else delete brAdd;
      } else if (brType == "Double_t") {
        // Double_t branch
        Double_t* brAdd = new Double_t(0.);
        setAdd = tree.SetBranchAddress(real->GetName(), brAdd);
        if (setAdd == 0) doubleMap[TString(real->GetName())] = brAdd;
        else delete brAdd;
      }

      // check if branch address was set
      if (setAdd != 0 || (brType != "Float_t" && brType != "Double_t")) {
        cout << "P2VV - WARNING: TreeToRooDataSet(): branch \""
          << real->GetName() << "\" is not of type Double_t or Float_t"
          << endl;
        obsSet.remove(*arg);
        continue;
      }

      // build RealVar selection string
      selStr += TString("&&") + real->GetName() + ">=";
      selStr += real->getMin();
      selStr += TString("&&") + real->GetName() + "<=";
      selStr += real->getMax();
    } else if (cat != 0) {
      // set branch address to Int_t pointer
      Int_t* brAdd = new Int_t(0);
      Int_t setAdd = tree.SetBranchAddress(cat->GetName(), brAdd);
      if (setAdd == 0) {
        intMap[TString(cat->GetName())] = brAdd;
      } else {
        cout << "P2VV - WARNING: TreeToRooDataSet(): branch \""
          << cat->GetName() << "\" is not of type Int_t" << endl;
        delete brAdd;
        obsSet.remove(*arg);
        continue;
      }

      // build category selection string
      TString catSelStr("");
      std::auto_ptr<TIterator> catIt(cat->typeIterator());
      RooCatType* type(0);
      while ((type = dynamic_cast<RooCatType*>(catIt->Next())) != 0) {
        catSelStr += TString("||") + cat->GetName() + "==";
        catSelStr += type->getVal();
      }
      selStr += TString("&&(") + catSelStr(2, catSelStr.Length() - 2) + ")";
    }
  }

  // enable branches for selection
  std::vector<TBranch*> disabledBrList;
  brIt->Reset();
  while ((br = dynamic_cast<TBranch*>(brIt->Next())) != 0) {
    if (!tree.GetBranchStatus(br->GetName())
        && selStr.Contains(br->GetName())) {
      disabledBrList.push_back(br);
      br->SetStatus(kTRUE);
    }
  }

  // get selected entries
  selStr = selStr.BeginsWith("&&") ? selStr(2, selStr.Length() - 2) : selStr;
  tree.Draw(">>elist", selStr.Data(), "entrylist");
  TEntryList *entryList = static_cast<TEntryList*>(gDirectory->Get("elist"));
  cout << "P2VV - INFO: TreeToRooDataSet(): " << entryList->GetN()
      << "/" << tree.GetEntries() << " entries selected in tree";
  if (origDataSet != 0)
    cout << ", " << numEntr << " entries in data set" << endl;
  else
    cout << endl;
  if (origDataSet != 0 && ((origIndex != 0 && numEntr > entryList->GetN())
      || (origIndex == 0 && numEntr != entryList->GetN()))) {
    cout << "P2VV - ERROR: TreeToRooDataSet(): number of entries in original data set ("
        << numEntr
        << ") is not equal to number of selected entries in the tree ("
        << entryList->GetN() << ")" << endl;
    return 0;
  }

  // disable branches that were only enabled for selection
  for (std::vector<TBranch*>::iterator it = disabledBrList.begin();
      it != disabledBrList.end(); ++it)
    (*it)->SetStatus(kFALSE);
  disabledBrList.clear();

  // get weight variable
  RooRealVar* weight = 0;
  if (wStr.Length() > 0) {
    weight = dynamic_cast<RooRealVar*>(obsSet.find(wStr));
    if (weight == 0)
      cout << "P2VV - WARNING: TreeToRooDataSet(): no RooRealVar named "
          << wStr << " found in set of observables: no event weights applied"
          << endl;
  }

  // create data set
  RooDataSet* dataSet = 0;
  if (weight != 0) {
    dataSet = new RooDataSet(dsName, dsTitle, obsSet, wStr);
  } else {
    dataSet = new RooDataSet(dsName, dsTitle, obsSet);
  }

  for (Long64_t it = 0; it < numEntr; ++it) {
    // get entries in tree and original data set
    if (origDataSet != 0) {
      origDataSet->get((Int_t)it);
      Long64_t entry(0);
      if (origIndex != 0) entry = (Long64_t)origIndex->getVal();
      else entry = it;
      if (!entryList->Contains(entry)) {
        cout << "P2VV - ERROR: TreeToRooDataSet(): entry " << entry
            << " in tree is not selected" << endl;
        delete dataSet;
        return 0;
      }
      tree.GetEntry(entry);
    } else {
      if (!entryList->Contains(it)) continue;
      tree.GetEntry(it);
    }

    // set values of observables in data set
    obsSetIter.Reset();
    while ((arg = static_cast<RooAbsArg*>(obsSetIter.Next())) != 0) {
      RooRealVar* real = dynamic_cast<RooRealVar*>(arg);
      RooAbsArg* origArg = origObsSet->find(arg->GetName());
      if (origArg != 0) {
        // get value from original data set
        if (real != 0) {
          real->setVal(dynamic_cast<RooRealVar*>(origArg)->getVal());
        } else {
          RooCategory* cat = dynamic_cast<RooCategory*>(arg);
          cat->setIndex(dynamic_cast<RooCategory*>(origArg)->getIndex());
        }
      } else if (index == 0 || arg->GetName() != index->GetName()) {
        // get value from tree
        if (real != 0) {
          std::map<TString,Float_t*>::iterator
              fIt(floatMap.find(arg->GetName()));
          if (fIt != floatMap.end()) real->setVal(*fIt->second);
          else real->setVal(*doubleMap[arg->GetName()]);
        } else {
          RooCategory* cat = dynamic_cast<RooCategory*>(arg);
          cat->setIndex(*intMap[arg->GetName()]);
        }
      } else if (index != 0) {
        index->setVal((Double_t)it);
      }
    }

    // add row to data set
    if (weight != 0)
      dataSet->add(obsSet, weight->getVal());
    else
      dataSet->add(obsSet);
  }

  // delete branch addresses
  tree.ResetBranchAddresses();
  for (std::map<TString,Double_t*>::iterator it = doubleMap.begin();
      it != doubleMap.end(); ++it)
    delete it->second;
  for (std::map<TString,Float_t*>::iterator it = floatMap.begin();
      it != floatMap.end(); ++it)
    delete it->second;
  for (std::map<TString,Int_t*>::iterator it = intMap.begin();
      it != intMap.end(); ++it)
    delete it->second;

  // delete dummy observables set
  if (origDataSet == 0) delete origObsSet;

  // reset branch status
  brIt->Reset();
  std::vector<Bool_t>::iterator brStatIt(brStatusList.begin());
  while ((br = dynamic_cast<TBranch*>(brIt->Next())) != 0) {
    br->SetStatus(*brStatIt);
    ++brStatIt;
  }

  // return data set
  return dataSet;
}

void addHelicityAnglesToTree(TTree& tree, 
			     const TString& posHadrName, const TString& negHadrName, 
			     const TString& posLeptName, const TString& negLeptName,
			     const Double_t& posHadrMass, const Double_t& negHadrMass, 
			     const Double_t& posLeptMass, const Double_t& negLeptMass,
			     const TString& helcosthetaK_brName, const TString& helcosthetaL_brName, const TString& helphi_brName,
			     const TString& trackMomType)
{// TODO: Use defaul values for the hel angles branch names

  // new branch addresses
  Double_t* helcosthetaK_address = new Double_t(0);
  Double_t* helcosthetaL_address = new Double_t(0) ;
  Double_t* helphi_address       = new Double_t(0);

  // create new branches
  TBranch* helcosthetaK_branch = tree.Branch(helcosthetaK_brName, helcosthetaK_address, helcosthetaK_brName + "/D");
  TBranch* helcosthetaL_branch = tree.Branch(helcosthetaL_brName, helcosthetaL_address, helcosthetaL_brName + "/D");
  TBranch* helphi_branch       = tree.Branch(helphi_brName,       helphi_address,       helphi_brName + "/D"      );

  // existing branch addresses
  // if(trackMomType=='F'){
  //   std::vector<Float_t> posHadrMomentum(3);
  //   std::vector<Float_t> negHadrMomentum(3);
  //   std::vector<Float_t> posLeptMomentum(3);
  //   std::vector<Float_t> negLeptMomentum(3);
  // };
  // else {
  //   std::vector<Double_t> posHadrMomentum(3);
  //   std::vector<Double_t> negHadrMomentum(3);
  //   std::vector<Double_t> posLeptMomentum(3);
  //   std::vector<Double_t> negLeptMomentum(3);
  // };

  std::vector<Float_t> posHadrMomentum(3);
  std::vector<Float_t> negHadrMomentum(3);
  std::vector<Float_t> posLeptMomentum(3);
  std::vector<Float_t> negLeptMomentum(3);

  // set branch addresses
  //TODO: Make the suffixes more flexible
  tree.SetBranchAddress( posHadrName + "_PX", &posHadrMomentum[0]);
  tree.SetBranchAddress( posHadrName + "_PY", &posHadrMomentum[1]);
  tree.SetBranchAddress( posHadrName + "_PZ", &posHadrMomentum[2]);

  tree.SetBranchAddress( negHadrName + "_PX", &negHadrMomentum[0]);
  tree.SetBranchAddress( negHadrName + "_PY", &negHadrMomentum[1]);
  tree.SetBranchAddress( negHadrName + "_PZ", &negHadrMomentum[2]);

  tree.SetBranchAddress( posLeptName + "_PX", &posLeptMomentum[0]);
  tree.SetBranchAddress( posLeptName + "_PY", &posLeptMomentum[1]);
  tree.SetBranchAddress( posLeptName + "_PZ", &posLeptMomentum[2]);

  tree.SetBranchAddress( negLeptName + "_PX", &negLeptMomentum[0]);
  tree.SetBranchAddress( negLeptName + "_PY", &negLeptMomentum[1]);
  tree.SetBranchAddress( negLeptName + "_PZ", &negLeptMomentum[2]);

  
  // calculate helicity angles and fill branches
  for (Long64_t it = 0; it < tree.GetEntries(); ++it)
    { 
      tree.GetEntry(it);
      // TODO: Do all these a lambda functions.
    
      // positive Hadron four momentum
      TVector3 posHadrMomVector = TVector3(posHadrMomentum[0], posHadrMomentum[1], posHadrMomentum[2]);
      Double_t posHadronEnergy = TMath::Sqrt( posHadrMass*posHadrMass + posHadrMomVector.Mag2() );
      TLorentzVector posHadrFourVect = TLorentzVector(posHadrMomVector, posHadronEnergy);

      // negative Hadron four momentum
      TVector3 negHadrMomVector = TVector3(negHadrMomentum[0], negHadrMomentum[1], negHadrMomentum[2]);
      Double_t negHadronEnergy = TMath::Sqrt( negHadrMass*negHadrMass + negHadrMomVector.Mag2() );
      TLorentzVector negHadrFourVect = TLorentzVector(negHadrMomVector, negHadronEnergy);

      // positive Lepton four momentum
      TVector3 posLeptMomVector = TVector3(posLeptMomentum[0], posLeptMomentum[1], posLeptMomentum[2]);
      Double_t posLeptonEnergy = TMath::Sqrt( posLeptMass*posLeptMass + posLeptMomVector.Mag2() );
      TLorentzVector posLeptFourVect = TLorentzVector(posLeptMomVector, posLeptonEnergy);

      // negative Lepton four momentum
      TVector3 negLeptMomVector = TVector3(negLeptMomentum[0], negLeptMomentum[1], negLeptMomentum[2]);
      Double_t negLeptonEnergy = TMath::Sqrt( negLeptMass*negLeptMass + negLeptMomVector.Mag2() );
      TLorentzVector negLeptFourVect = TLorentzVector(negLeptMomVector, negLeptonEnergy);
    
      // call helicity angles calculation function
      std::vector<Double_t> helAngles = HelicityAngles(posHadrFourVect, negHadrFourVect, posLeptFourVect, negLeptFourVect);    
    
      *helcosthetaK_address = helAngles[0];
      *helcosthetaL_address = helAngles[1];
      *helphi_address       = helAngles[2];
    
      helcosthetaK_branch->Fill();
      helcosthetaL_branch->Fill();
      helphi_branch->Fill();

    }
  tree.FlushBaskets();
}


std::vector<double> HelicityAngles(TLorentzVector Kplus_P, TLorentzVector Kminus_P, TLorentzVector muplus_P,
                                   TLorentzVector muminus_P )
{ // Calculation based on the ANA-2012-067-v3
  // IMPORTANT NOTE: To Generalisze remember that: Kplus/Kminus MUST be the postive/hadron. Same oes for the leptons.

  // Bs, KK, mm momenta 4 vectors 
  TLorentzVector KK_P   = Kplus_P + Kminus_P;
  TLorentzVector mm_P   = muplus_P + muminus_P;
  TLorentzVector KKmm_P = KK_P + mm_P;

  // Unit vector along mumu direction in the KK mass r.f.
  muplus_P.Boost( - KK_P.BoostVector() );
  muminus_P.Boost( - KK_P.BoostVector() );
  TVector3 e_KK = - (muplus_P + muminus_P).Vect().Unit();

  // Boost the muons back to lab frame
  muplus_P.Boost( KK_P.BoostVector() );
  muminus_P.Boost( KK_P.BoostVector() );

  // Unit vector along KK direction in the mm mass r.f.
  Kplus_P.Boost( - mm_P.BoostVector() );
  Kminus_P.Boost( - mm_P.BoostVector() );
  TVector3 e_mm = - (Kplus_P+Kminus_P).Vect().Unit();
  // Boost the Kaons back to lab frame
  Kplus_P.Boost( mm_P.BoostVector() );
  Kminus_P.Boost( mm_P.BoostVector() );

  // Unit vector along KK direction in the mm mass r.f.
  Kplus_P.Boost( - KKmm_P.BoostVector() );
  Kminus_P.Boost( - KKmm_P.BoostVector() );
  muplus_P.Boost( - KKmm_P.BoostVector() );
  muminus_P.Boost( - KKmm_P.BoostVector() );
  TVector3  e_KKmm = (muplus_P + muminus_P).Vect().Unit();

  // Perpenticular vectors to KK and mm planes in the KKmmm r.f.
  TVector3 eta_KK = ( Kplus_P.Vect().Cross( Kminus_P.Vect()) ).Unit();
  TVector3 eta_mm = ( muplus_P.Vect().Cross( muminus_P.Vect()) ).Unit();
  
  Kplus_P.Boost( KKmm_P.BoostVector() );
  Kminus_P.Boost( KKmm_P.BoostVector() );
  muplus_P.Boost( KKmm_P.BoostVector() );
  muminus_P.Boost( KKmm_P.BoostVector() );

  // Helicity angles. 
  Kplus_P.Boost( - KK_P.BoostVector() );
  muplus_P.Boost( - mm_P.BoostVector() );

  std::vector<Double_t> angles(3);
  angles[0] = ( Kplus_P.Vect().Unit()  ).Dot(e_KK);
  angles[1] = ( muplus_P.Vect().Unit() ).Dot(e_mm);
      
  if ( eta_KK.Cross(eta_mm).Dot(e_KKmm) > 0 ) {         // sinphi = eta_KK.Cross(eta_mm).Dot(e_KKmm); 
    angles[2] = + TMath::ACos( eta_KK.Dot(eta_mm) ) ;}  // cosphi = eta_KK.Dot(eta_mm);
  else                                        { 
    angles[2] = - TMath::ACos( eta_KK.Dot(eta_mm) ) ;}

  return angles;
}

bool GetOwnership(_object* obj) {
   return (reinterpret_cast<TFakeObjectProxy*>(obj))->fFlags & 0x0001;
}

Double_t getRooRealMaxVal(RooAbsReal& function, const RooArgList& scanVars,
    const vector<Int_t>& numPoints) {
  Double_t maxVal = -RooNumber::infinity();

  // check number of scan variables
  Int_t numVars(scanVars.getSize());
  if (numVars < 1) {
    std::cerr << "P2VV - ERROR: getRooRealMaxVal(): no scan variables found"
         << std::endl;
    assert(numVars > 0);
    return maxVal;
  }
  if ((Int_t)numPoints.size() != numVars) {
    std::cerr << "P2VV - ERROR: getRooRealMaxVal(): number of variables in \"numPoints\" ("
        << numPoints.size() << ") does not match number of scan variables ("
         << numVars << ")" << std::endl;
    assert((Int_t)numPoints.size() == numVars);
    return maxVal;
  }

  // create vector of variable values
  RooRealVar** varAdds = new RooRealVar*[numVars];
  vector<Double_t>* varVals = new vector<Double_t>[numVars];
  vector<Double_t>::const_iterator* varValsIts
      = new vector<Double_t>::const_iterator[numVars];
  for (Int_t varIt = 0; varIt < numVars; ++varIt) {
    // get variable
    RooRealVar* var = dynamic_cast<RooRealVar*>(scanVars.at(varIt));
    if (var == 0) {
      std::cerr << "P2VV - ERROR: getRooRealMaxVal(): variable \""
          << scanVars.at(varIt)->GetName() << "\" is not a RooRealVar"
          << std::endl;
      assert(var != 0);
      return maxVal;
    }

    // get number of points for variable
    Int_t nPts = numPoints.at(varIt);
    if (nPts < 1) {
      std::cerr << "P2VV - ERROR: getRooRealMaxVal(): number of scan points for variable \""
          << var->GetName() << "\" is smaller than one" << std::endl;
      assert(nPts > 0);
      return maxVal;
    }

    // create vector of variable values
    if (!var->hasMin() || !var->hasMax()) {
      std::cerr << "P2VV - ERROR: getRooRealMaxVal(): range of variable \""
          << var->GetName() << "\" is infinite" << std::endl;
      assert(var->hasMin() && var->hasMax());
      return maxVal;
    }
    vector<Double_t> vals;
    if (nPts == 1) {
      vals.push_back(0.5 * (var->getMin() + var->getMax()));
    } else {
      vals.push_back(var->getMin());
      vals.push_back(var->getMax());
    }
    for (Int_t valIt = 1; valIt < nPts - 1; ++valIt) {
      vals.push_back(var->getMin() + ((Double_t)valIt) / (Double_t)(nPts - 1)
          * (var->getMax() - var->getMin()));
    }

    // fill vectors
    varAdds[varIt]    = var;
    varVals[varIt]    = vals;
    varValsIts[varIt] = varVals[varIt].begin();
    var->setVal(vals.front());
  }

  // find maximum of function
  std::cout << "P2VV - INFO: getRooRealMaxVal(): find maximum of function \""
      << function.GetName() << "\":"<< std::endl << "    ";
  function.Print();

  Int_t varIt(0);
  Long64_t itCount(0);
  Long64_t numIts0(varVals[0].size());
  while (varIt < numVars) {
    // print number of iterations
    if (itCount % 100000 == 0) {
      std::cout << "    iteration " << itCount << ": maximum = " << maxVal
          << std::endl;
    }

    // loop over values of first variable and find function maximum
    while (varValsIts[0] != varVals[0].end()) {
      varAdds[0]->setVal(*varValsIts[0]);
      Double_t funcVal(function.getVal());
      if (funcVal > maxVal) maxVal = funcVal;
      ++varValsIts[0];
    }
    varValsIts[0] -= varVals[0].size();

    // set values of variables other than the first
    varIt = 1;
    while (varIt < numVars) {
      if (++varValsIts[varIt] != varVals[varIt].end()) {
        varAdds[varIt]->setVal(*varValsIts[varIt]);
        break;
      }
      varValsIts[varIt] -= varVals[varIt].size();
      varAdds[varIt]->setVal(*varValsIts[varIt]);
      ++varIt;
    }

    itCount += numIts0;
  }
  std::cout << "    maximum = " << maxVal << std::endl;

  delete[] varAdds;
  delete[] varVals;
  delete[] varValsIts;
  return maxVal;
}
