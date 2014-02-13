// $Id: $
#ifndef FUNCTIONS_H 
#define FUNCTIONS_H 1

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include "Rtypes.h"
#include <TLorentzVector.h>


class TH1l;
class TTree;
class RooDataSet;

void addSWeightToTree(const std::vector<double>& weights, TTree& tree, const std::string& branch_name);
void addSWeightToTree(double* weights, size_t n, TTree& tree, const std::string& branch_name);
void addSWeightToTree(const RooDataSet& ds, TTree& tree, const std::string& branchname,
                      const std::string& cut = std::string("1"));
void addIntegerToTree(TTree& tree, Int_t value, const char* branchName);
void addFloatToTree(TTree& tree, Double_t value, const char* branchName);
void copyFloatInTree(TTree& tree, const char* inBranch, const char* outBranch);
void addCategoryToTree(TTree& tree, const char* inBranch,
    const char* catBranch, std::vector<Double_t> boundaries,
    std::vector<Int_t> indices);

void addVertexErrors(TTree* tree, const std::list<RooDataSet*>& dss, const std::string& cut);

TTree* RooDataSetToTree(const RooDataSet& dataSet, const char* name = 0,
      const char* title = 0, const char* weightName = 0,
      const char* branchList = 0, Bool_t RooFitFormat = kTRUE);

RooDataSet* TreeToRooDataSet(TTree& tree, const RooArgSet& observables,
      const char* name = 0, const char* title = 0, const char* cuts = 0,
      const char* indexName = 0, const char* weightName = 0,
      RooDataSet* origDataSet = 0);

std::vector<double> HelicityAngles(TLorentzVector Kplus_P, TLorentzVector Kminus_P, 
				   TLorentzVector muplus_P, TLorentzVector muminus_P);


struct _object;

struct TFakeObjectProxy {
   ssize_t fRefCnt;          // PyObject_HEAD
   void* fPyType;            // PyObject_HEAD
   void* fRootObj;
   int fFlags;
};

bool GetOwnership(_object* obj);

Double_t getRooRealMaxVal(RooAbsReal& function, const RooArgList& scanVars,
    const std::vector<Int_t>& numPoints);

#endif // FUNCTIONS_H
