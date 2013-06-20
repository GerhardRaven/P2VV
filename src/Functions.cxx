// STD & STL
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <list>

// ROOT
#include <TH1.h>
#include <TTree.h>
#include <TMatrixT.h>
#include <TEntryList.h>
#include <TDirectory.h>

// RooFit
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>

namespace {
   using std::cout;
   using std::endl;
   using std::list;
}

double sigmaFromFT( const TH1& h1, const double dMs, const double dMsErr, std::ostream& out)
{
   double sum(0), sumdt(0), sumdterr(0), sumcos(0), sumx2(0), error(0) ;
   for(int i=1; i<=h1.GetNbinsX() ; ++i) {
      double c = h1.GetBinContent( i ) ;
      double e = h1.GetBinError( i ) ;
      double x = h1.GetBinCenter( i ) ;
      double dt = h1.GetBinWidth( i ) ;
      sum += c ;
      sumdt += c * dt ;
      sumdterr += e * e * dt * dt;

      double cosv = cos( - dMs * x );
      sumcos += c * cosv * dt ;
      sumx2 += c * x * x ;

      double coserr = fabs(sin(dMsErr)) * dMsErr;
      error += (c * c * coserr * coserr +  cosv * cosv * e * e) * dt * dt;
   }

   error = sqrt(error / (sumdt * sumdt) + sumdterr * sumcos * sumcos / pow(sumdt, 4));

   double rms = sqrt(sumx2/sum) ;
   double D = sumcos / sumdt ;
   double sigma = sqrt( -2*log(D) ) / dMs ;

   out << sum << " " << sumdt << " " << sumcos << " " << sumx2 << std::endl;
   out << "RMS of input histogram: " << rms<< std::endl
       << "If distribution were Gaussian, dilution is: " << exp(-0.5*rms*rms*dMs*dMs) << std::endl 
       << "Dilution from FT: " << D << " +- " << error << std::endl
       << "Corresponding Gaussian resolution: " << sigma << std::endl;

   return D;
}

void addSWeightToTree(const RooDataSet& ds, TTree& tree, const std::string& branch_name,
                      const std::string& cut)
{
   Double_t w = 0;

   tree.Draw(">>elist", cut.c_str(), "entrylist");
   TEntryList *cut_list = static_cast<TEntryList*>(gDirectory->Get( "elist" ));

   //#       Create the output Branch
   std::string branch_def = branch_name + "/D";
   TBranch* branch = tree.Branch(branch_name.c_str(), &w, branch_def.c_str());

   Long64_t nds(ds.numEntries());
   assert(nds == cut_list->GetN());

   std::vector<Double_t> weights(tree.GetEntries(), 0.);

   Int_t j = 0;
   for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
      if (cut_list->Contains(i)) {
         ds.get(j++);
         weights[i] = ds.weight();
      }
   }

   for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
      w = weights[i];
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
