// $Id: $
#ifndef FUNCTIONS_H 
#define FUNCTIONS_H 1

#include <iostream>
#include <string>
#include <list>

class TH1l;
class TTree;
class RooDataSet;

double sigmaFromFT( const TH1& h1, const double dMs, const double dMsErr, std::ostream& out = std::cout );

void addSWeightToTree(const RooDataSet& ds, TTree& tree, const std::string& branchname,
                      const std::string& cut = std::string("1"));

void addVertexErrors(TTree* tree, const std::list<RooDataSet*>& dss, const std::string& cut);

#endif // FUNCTIONS_H
