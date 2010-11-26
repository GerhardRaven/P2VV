#include <iostream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TROOT.h>
#include <TKey.h>
#include <TEventList.h>

// ROOFIT CLASSES
#include <RooFit.h>
#include <RooConstVar.h>
#include <RooHist.h>
#include <RooCategory.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooAddModel.h"
#include "RooResolutionModel.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooBDecay.h"
#include "TPaveLabel.h"
#include "RooStats/SPlot.h"
#include "RooWorkspace.h"

#include "TRegexp.h"

using namespace RooFit;

double inline sqr(double x) { return x*x ; }
// -------------------------------------------------------------------------------------------------

TTree* findFirstTree( TDirectory* dir )
{
  TList* keys = dir->GetListOfKeys() ;
  TIterator* it = keys->MakeIterator() ;
  TTree* tree(0) ;
  TObject* obj ;
  TDirectory* subdir(0) ;
  while( (obj = it->Next() ) && tree==0 ) {
    TKey* key = dynamic_cast<TKey*>(obj) ;
    obj = key->ReadObj() ;
    tree = dynamic_cast<TTree*>(obj) ;
    if( TString(tree->GetName()) != TString("dataset") ) tree = 0 ;
    if( !tree && (subdir = dynamic_cast<TDirectory*>(obj)) ) 
      tree = findFirstTree(subdir) ;
  }
  delete it ;
  return tree ;
}

enum DecayMode { UnknownMode, JPsiK, JPsiKs, JPsiPhi, JPsiKstar } ;

void defineObservables(RooWorkspace& ws)
{
  RooArgSet obs ;
  obs.add(*new RooRealVar("m","B mass",5200,5450) ) ;
  obs.add(*new RooRealVar("t","B proper time",-2,10)) ;
  obs.add(*new RooRealVar("sigmat","B proper time error",0.005,0.2)) ;
  obs.add(*new RooRealVar("mdau1","dau1 mass",3097-60,3097+60)) ;
  //obs.add(*new RooRealVar("mdau2","dau2 mass",497-25,497+25)) ;
  obs.add(*new RooRealVar("mdau2","dau2 mass",0,2000)) ;
  obs.add(*new RooRealVar("mKS","KS mass",497-25,497+25)) ;
  obs.add(*new RooRealVar("mKK","KK mass",1019.455-20,1019.455+20)) ;

  obs.add(*new RooRealVar("pt","B transverse momentum",1000,0,500000)) ;
  obs.add(*new RooRealVar("chi2dof","B chi2/dof",0,5)) ;
  obs.add(*new RooRealVar("trcostheta","cos(theta_tr)",-1,1)) ;
  obs.add(*new RooRealVar("trcospsi","cos(theta_psi)",-1,1)) ;
  obs.add(*new RooRealVar("trphi","phi_tr",-M_PI,M_PI)) ;


  RooCategory * cat = new RooCategory("decaytype","B decay type") ;
  cat->defineType("J/psiKsLL",31) ;
  cat->defineType("J/psiKsDD",30) ;
  cat->defineType("J/psiKplus",10) ;
  cat->defineType("J/psiKminus",11) ;
  cat->defineType("J/psiPhi",40) ;
  cat->defineType("J/psiKstar",20) ;
  cat->defineType("J/psiKstarbar",21) ;
  
  obs.add( *cat ) ;
  //obs.add( *(new RooCategory("runnr","Run number") ) ) ;
  
  ws.import( obs ) ;
  //ws.defineSet("observables",obs) ;
  ws.defineSet("observables","m,t,sigmat,mdau1,mdau2,pt,chi2dof,decaytype,mKK,mKS,trcostheta,trcospsi,trphi") ;
  ws.var("sigmat")->setBins(40) ;
  
  cat = new RooCategory("decaymode","B decay mode") ;
  cat->defineType("Unknown",UnknownMode) ;
  cat->defineType("JPsiKs",JPsiKs) ;
  cat->defineType("JPsiK",JPsiK) ;
  cat->defineType("JPsiKstar",JPsiKstar) ;
  cat->defineType("JPsiPhi",JPsiPhi) ;
  ws.import(*cat) ;
}

void definePDF(RooWorkspace& ws)
{  

  //RooMsgService::instance().addStream(INFO,Topic(Integration)) ;

  // define some constants
  ws.factory("const_zero[0]") ;
  ws.factory("const_one[1]") ;

  // signal B mass pdf
  ws.factory("Gaussian::m_sig(m,m_sig_mean[5280,5200,5400],m_sig_sigma[10,3,30])");

  // background B mass pdf
  ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.1,0.1])") ;

  // signal phi mass pdf
  ws.factory("Voigtian::mphi_phisig(mKK,mphi_phi_mean[1919.455],mphi_phi_width[4.26],mphi_phi_sigma[1,0.1,10])") ;

  // background phi mass pdf (would like to use a nice threshold function here ... RooDstD0BG seems a bit complicated
  ws.factory("DstD0BG::mphi_combbkg(mKK,mphi_bkg_m0[987.4],mphi_bkg_C[10,1,30],mphi_bkg_B[1,0.1,10],const_zero)") ;

  // signal J/psi mass pdf
  ws.factory("CBShape::mpsi_sig(mdau1,mpsi_sig_mean[3097,3085,3110],mpsi_sig_sigma[13.1,5,20],"
	     "mpsi_sig_alpha[1.36,0.5,3],mpsi_sig_n[3])");
  //ws.factory("BifurGauss::mpsi_sig(mdau1,mpsi_sig_mean[3097,3085,3110],mpsi_sig_sigmaL[18,5,30],"
  // "m_psi_sigmaR[13,5,30])");
  //ws.factory("Gaussian::mpsi_sig(mdau1,mpsi_sig_mean[3097,3085,3110],m_sigpsi_sigma[13,5,30])");

  // background J/psi mass pdf
  ws.factory("Exponential::mpsi_bkg(mdau1,mpsi_bkg_exp[-0.0002,-0.01,0.01])") ;

  // signal resolution model
  ws.factory("GaussModel::tres_sig_g1(t,tres_sig_m1[0,-0.2,0.2],tres_sig_s1[1.1,0.3,2], const_one, sigmat)");
  ws.factory("GaussModel::tres_sig_g2(t,tres_sig_m1,            tres_sig_s2[2.0,1.5,10],const_one, sigmat)");
  //ws.factory("GaussModel::tres_sig_g2(t,tres_sig_m1,          2*tres_sig_s1,const_one, sigmat)");
  ws.factory("GaussModel::tres_sig_g3(t,tres_sig_m1,            tres_sig_s3[0.54,0.1,3.0])");
  ws.factory("AddModel::tres_sig({tres_sig_g3,tres_sig_g2,tres_sig_g1},{tres_sig_f3[0.001,0.00,0.01],tres_sig_f2[0.2,0.01,1]})");

  // background resolution model
  ws.factory("GaussModel::tres_bkg_g1(t,tres_bkg_m1[0,-0.2,0.2],tres_bkg_s1[1.5,0.3,3], const_one,sigmat)");
  ws.factory("GaussModel::tres_bkg_g2(t,tres_bkg_m1,            tres_bkg_s2[3.0,0.9,10],const_one,sigmat)");
  //ws.factory("GaussModel::tres_bkg_g2(t,tres_bkg_m1,            2*tres_bkg_s2,const_one,sigmat)");
  ws.factory("GaussModel::tres_bkg_g3(t,tres_bkg_m1,            tres_bkg_s3[0.54,0.1,3.0])");
  ws.factory("AddModel::tres_bkg({tres_bkg_g3,tres_bkg_g2,tres_bkg_g1},{tres_bkg_f3[0.002,0.00,0.01],tres_bkg_f2[0.2,0.01,1]})");

  // define sigmat pdf
  ws.import(*(new RooDataHist("data_sigt","hist Err Per Ev",
			      *ws.var("sigmat"),*ws.data("data")))) ;
  ws.factory("HistPdf::pdf_stE(sigmat,data_sigt)");
 
  // signal propertime
  ws.factory("Decay::t_sig(t,t_sig_tau[1.5,1.0,2.0],tres_sig,SingleSided)");
  
  // prompt psi propertime
  ws.factory("Decay::t_prompt_sl(t,const_zero,tres_sig,SingleSided)");
  ws.factory("Decay::t_prompt_ml(t,t_prompt_ml_tau[0.23,0.05,1.5],tres_sig,SingleSided)");
  ws.factory("Decay::t_prompt_ll(t,t_prompt_ll_tau[1.5,1.5,1.5],tres_sig,SingleSided)");
  ws.factory("SUM::t_prompt(t_prompt_fll[0.004,0,1]*t_prompt_ll,t_prompt_fml[0.02,0,1]*t_prompt_ml,t_prompt_sl)");
  ws.factory("PROD::tconv_prompt(t_sig|sigmat)");
  
  // background propertime
  ws.factory("Decay::t_bkg_sl(t,const_zero,tres_bkg,SingleSided)");
  ws.factory("Decay::t_bkg_ml(t,t_bkg_ml_tau[0.23,0.05,1.5],tres_bkg,SingleSided)");
  ws.factory("Decay::t_bkg_ll(t,t_bkg_ll_tau[1.5,1.5,1.5],tres_bkg,SingleSided)");
  ws.factory("SUM::t_bkg(t_bkg_fll[0.004,0,1]*t_bkg_ll,t_bkg_fml[0.02,0,1]*t_bkg_ml,t_bkg_sl)");

  // now combine ...
  //ws.factory("PROD::sig_pdf( m_sig, mpsi_sig, PROD(t_sig|sigmat,pdf_stE))");
  //ws.factory("PROD::prompt_pdf( m_bkg, mpsi_sig, PROD(t_prompt|sigmat,pdf_stE))");
  //ws.factory("PROD::bkg_pdf( m_bkg, mpsi_bkg, PROD(t_bkg|sigmat,pdf_stE))");
  ws.factory("PROD::sig_pdf( m_sig, mpsi_sig, PROD(t_sig|sigmat)");
  ws.factory("PROD::prompt_pdf( m_bkg, mpsi_sig, PROD(t_prompt|sigmat))");
  ws.factory("PROD::bkg_pdf( m_bkg, mpsi_bkg, PROD(t_bkg|sigmat))");
  ws.factory("SUM::pdf(Nsig[0,100000]*sig_pdf,Nprompt[0,1000000]*prompt_pdf,Nbkg[0,1000000]*bkg_pdf)");

  // also build the mass-only pdf
  ws.factory("SUM::masspdf(Nsig*PROD( m_sig, mpsi_sig),"
	     "Nprompt*PROD( m_bkg, mpsi_sig),Nbkg*PROD( m_bkg, mpsi_bkg))") ;

  // also build the mass-only pdf
  ws.factory("SUM::psimasspdf(Nprompt*mpsi_sig,Nbkg*mpsi_bkg)") ;
  
  // since we haev the data, make sure that the yields are reasonably well initialized
  int numentries = (static_cast<RooDataSet*>(ws.data("data")))->numEntries() ;
  double signal  = 1  ;
  double psifrac = 0.54 ; 
  ws.var("Nsig")->setVal(signal) ;
  ws.var("Nprompt")->setVal(psifrac*(numentries-signal)) ;
  ws.var("Nbkg")->setVal((1-psifrac)*(numentries-signal)) ;

  ws.defineSet( "yields","Nsig,Nprompt,Nbkg") ;

}

void readParameters(RooWorkspace& ws) {
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws.data("data")) ;
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws.function("pdf")) ; 
  RooArgSet* pars = pdf->getParameters(*data) ;
  
  TString filename = TString(ws.cat("decaymode")->getLabel()) + "Parameters.txt" ;

  std::cout << "Reading fit starting point from file: " << filename << std::endl ;
  pars->readFromFile( filename ) ;

  // to make the fit work, we need to scale the yields
  double Nsig = ws.var("Nsig")->getVal() ;
  double Nprompt = ws.var("Nprompt")->getVal() ;
  double Nbkg = ws.var("Nbkg")->getVal() ;
  double N = Nsig + Nprompt + Nbkg ;
  double norm = data->numEntries()/N ;
  
  ws.var("Nsig")->setVal( norm * Nsig ) ;
  ws.var("Nprompt")->setVal( norm * Nprompt ) ;
  ws.var("Nbkg")->setVal( norm * Nbkg ) ;
  ws.var("Nsig")->setError( norm * ws.var("Nsig")->getError() ) ;
  ws.var("Nprompt")->setError( norm * ws.var("Nprompt")->getError() ) ;
  ws.var("Nbkg")->setError( norm * ws.var("Nbkg")->getError() ) ;
}

void writeParameters(RooWorkspace& ws) {
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws.data("data")) ;
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws.function("pdf")) ; 
  RooArgSet* pars = pdf->getParameters(*data) ;
  TString filename = TString(ws.cat("decaymode")->getLabel()) + "ParametersOut.txt" ;
  pars->writeToFile( filename ) ;
}


// -------------------------------------------------------------------------------------------------
RooWorkspace* gMyWorkspace(0) ;

RooWorkspace* createWS(DecayMode mode)
{
  RooWorkspace* ws = new RooWorkspace("ws",kTRUE);

  if( gMyWorkspace ) delete gMyWorkspace ;
  gMyWorkspace = ws ;


  // define observables
  defineObservables(*ws) ;
  ws->cat("decaymode")->setIndex(mode) ;

  // load the data
  // Look in this directory for a tree
  TFile* inputfile = gDirectory->GetFile() ;
  TTree* tree = findFirstTree (inputfile ) ;
  if(!tree) {
    std::cout << "Couldn't find TTree ..." << std::endl ;
  } else {
    std::cout << "Found a tree: " << tree->GetName() << std::endl ;
  }

  TEventList* evtlist = dynamic_cast<TEventList*>(gDirectory->Get("evtlist")) ;
  if( evtlist ) delete evtlist ;
  tree->Draw(">>evtlist","","",2000000) ;
  evtlist = dynamic_cast<TEventList*>(gDirectory->Get("evtlist")) ;
  tree->SetEventList(evtlist) ;
  //TFile* f = new TFile("tmp.root","RECREATE") ;
  //TTree* subtree = tree->CopyTree("");
  TTree* subtree = tree ;

  //tree->Scan("m","","",5) ;
  RooDataSet* data = new RooDataSet("data","data",subtree,*ws->set("observables"),
				    "t==t && sigmat==sigmat && m==m && mdau1 == mdau1") ;
  //(const_cast<TTree*>(data->tree()))->Scan("m","","",5) ;
  //if(subtree != tree) delete subtree ;

  std::cout << "Created a dataset with " << data->numEntries() 
	    << " entries out of a tree with "
	    << tree->GetEntries() << " entries."
	    << std::endl ;

  ws->import( *data ) ;
  if( data->numEntries()==0) {
    delete ws ;
    return 0;
  }
  
  // define the pdf. since we use the sigmat projection to build the
  // sigmat pdf, we can only do this after we have imported the data
  definePDF(*ws) ;
  readParameters(*ws) ;
  
  std::cout << "Defined pdf" << std::endl ;
  {
    RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data")) ;
    RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws->function("pdf")) ; 
    RooArgSet* pars = pdf->getParameters(*data) ;
    pars->writeToFile( "parameters.txt") ;
    //pdf->printTree(std::cout) ;
  }

  inputfile->cd() ;

  ws->Print() ;

  return ws ;
}


RooWorkspace* readWorkspace(const char filename[]="fitresult.root")
{
  TFile* file = new TFile(filename) ;
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(file->Get("workspace")) ;
  delete gMyWorkspace ;
  if(ws) gMyWorkspace = ws ;
  else   gMyWorkspace = 0;
  return ws ;
}

void readResolution(const char* resofile) 
{
  if(!gMyWorkspace) return ;
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(gMyWorkspace->function("pdf")) ; 
  RooDataSet* data = dynamic_cast<RooDataSet*>(gMyWorkspace->data("data")) ;
  RooArgSet* parameters = pdf->getParameters(*data) ;
  parameters->readFromFile(resofile) ;
  parameters->writeToFile("tmp.txt") ;
}

size_t setConstant(RooWorkspace& ws, const TString& pattern, bool constant = true)
{
  size_t rc(0) ;
  RooArgList arglist(ws.allVars()) ;
  TRegexp rexp(pattern,kTRUE) ;
  for( int i=0; i<  arglist.getSize(); ++i) {
    RooAbsArg* arg = arglist.at(i) ;
        // copied from 
    if (TString(arg->GetName()).Index(rexp)>=0) {
      RooRealVar* rarg = dynamic_cast<RooRealVar*>(arg) ;
      if(rarg) {
	rarg->setConstant( constant ) ;
	++rc ;
      }
    }
  }
  return rc ;
}

void plot( RooWorkspace& ws) ;

void fit(DecayMode mode) 
{
  RooWorkspace* ws = createWS(mode) ;
  
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws->function("pdf")) ; 
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data")) ;
  //pdf->printCompactTree() ;
  //pdf->Print() ;

  RooAbsReal* nll = pdf->createNLL(*data,Extended(true)) ;
  //RooAbsReal* nll = pdf->createNLL(*data) ;
  std::cout << "nll value at starting point: " << nll->getVal() << std::endl ;

  RooFitResult *result(0) ;
  result = pdf->fitTo(*data,NumCPU(8),Extended(true),Minos(false),Save(true),Verbose(false));

  writeParameters(*ws) ;

  plot(*ws) ;

  //TFile file("fitresult.root","RECREATE") ;
  //ws->Write("workspace") ;
  //file.Write() ;
  //file.Close() ;
}


void fitMassHisto(DecayMode mode)
{
  RooWorkspace* ws = gMyWorkspace ;
  if( ws==0 ) {
    ws = createWS(mode) ;
  }
  
  //RooRealVar* m = ws->var("m") ;
  RooRealVar* mpsi = ws->var("mdau1") ;
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data")) ;
  //RooAbsPdf* masspdf = dynamic_cast<RooAbsPdf*>(ws->function("masspdf")) ;
  RooAbsPdf* psimasspdf = dynamic_cast<RooAbsPdf*>(ws->function("psimasspdf")) ;
  ws->var("Nsig")->setConstant() ;

  RooDataHist histdata("mhist","",RooArgSet(*mpsi),*data) ;
  psimasspdf->fitTo(histdata,NumCPU(8),Extended(true),Minos(false),Save(false));
  TCanvas* c1 = new TCanvas("c1","") ;
  c1->Divide(2,2) ;
  c1->cd(1) ;
  RooPlot* mpsiplot = mpsi->frame() ;
  histdata.plotOn(mpsiplot) ;
  psimasspdf->plotOn(mpsiplot) ;
  mpsiplot->Draw() ;
}

void plot( RooWorkspace& ws)
{
  gROOT->SetStyle("Plain");
  
  RooCmdArg lw = LineWidth(2.0);
  RooCmdArg xes = XErrorSize(0);
  
  TCanvas *c = new TCanvas("c","c",900,700);
  c->Divide(3,2);
  int sigcolor = kRed ;
  int psibkgcolor = kBlue ;
  int nonpsibkgcolor = kGreen ;

  //
  RooRealVar* t = ws.var("t") ;
  RooRealVar* st = ws.var("sigmat") ;
  RooRealVar* m = ws.var("m") ;
  RooRealVar* mpsi = ws.var("mdau1") ;
  RooAbsData* st_data = ws.data("data_sigt") ;
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws.data("data")) ;
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws.function("pdf")) ;

  // create some ranges
  double largeTime = 0.20;
  t->setRange("largeTime",largeTime,5);

  
  // ===========================================================================================================
  c->cd(1);
  TH2F *hist = data->createHistogram(*t,*m);
  hist->SetMarkerSize(0.3);
  hist->SetMarkerStyle(20);
  hist->SetStats(kFALSE);
  hist->GetXaxis()->SetTitle(t->getTitle(kTRUE));
  hist->GetYaxis()->SetTitle(m->getTitle(kTRUE));
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->SetTitle("");
  hist->Draw(); 

 // ===========================================================================================================
  // proper time in side bands 
  // ===========================================================================================================
  c->cd(2);
  gPad->SetLogy();

  RooPlot *_tb = t->frame(Bins(240));
  data->plotOn(_tb,Invisible());
  std::cout << "st_data: " << st_data << std::endl ;
  std::cout << "address: " << st_data->get()->find("sigmat") << " " << st << std::endl ;
  
  pdf->plotOn(_tb,ProjWData(*st_data),Components("sig_pdf"),LineColor(sigcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_tb,ProjWData(*st_data),Components("prompt_pdf"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_tb,ProjWData(*st_data),Components("bkg_pdf"),LineColor(nonpsibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_tb,ProjWData(*st_data),lw);
  data->plotOn(_tb,MarkerSize(0.5),xes);
  _tb->SetMinimum(0.1) ;
  _tb->SetTitle("");
  _tb->Draw();
  c->Update() ;

  
  // ===========================================================================================================
  c->cd(3);
  RooPlot *mf = m->frame(Bins(50));
  data->plotOn(mf,Invisible());
  pdf->plotOn(mf,ProjWData(*st_data),Components("bkg_pdf"),LineColor(nonpsibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mf,ProjWData(*st_data),Components("sig_pdf"),LineColor(sigcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mf,ProjWData(*st_data),Components("prompt_pdf"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mf,ProjWData(*st_data),lw);
  data->plotOn(mf,MarkerSize(0.7),xes);
  mf->Draw() ;
  c->Update() ;

 // ===========================================================================================================
  c->cd(4);
  
  char buf4[255];
  sprintf(buf4,"%s > %3.2f %s",(const char*)t->getTitle(),largeTime,(const char*)t->getUnit());
  RooPlot *_m = m->frame(Bins(50),Title(buf4));
  data->plotOn(_m,CutRange("largeTime"),Invisible());
  pdf->plotOn(_m,ProjWData(*st_data),ProjectionRange("largeTime"),Components("bkg_pdf"),LineColor(nonpsibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_m,ProjWData(*st_data),ProjectionRange("largeTime"),Components("sig_pdf"),LineColor(sigcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_m,ProjWData(*st_data),ProjectionRange("largeTime"),Components("prompt_pdf"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(_m,ProjWData(*st_data),ProjectionRange("largeTime"),lw);
  data->plotOn(_m,CutRange("largeTime"),MarkerSize(0.7),xes);
  _m->Draw() ;
  c->Update() ;


  // ===========================================================================================================
  // proper time in signal region
  // ===========================================================================================================
  c->cd(5);
  
  RooPlot *mpsiplot = mpsi->frame(Bins(50));
  data->plotOn(mpsiplot,Invisible());
  pdf->plotOn(mpsiplot,ProjWData(*st_data),Components("sig_pdf"),LineColor(sigcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mpsiplot,ProjWData(*st_data),Components("bkg_pdf"),LineColor(nonpsibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mpsiplot,ProjWData(*st_data),Components("prompt_pdf"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  pdf->plotOn(mpsiplot,ProjWData(*st_data),lw);
  data->plotOn(mpsiplot,MarkerSize(0.7),xes);
  mpsiplot->Draw() ;
  c->Update() ;
  
  
  c->cd(6) ;
  RooPlot *stplot = st->frame() ;
  st_data->plotOn( stplot ) ;
  stplot->Draw() ;

  //c->cd(6) ;
  //(const_cast<TTree*>(data->tree()))->Draw("chi2dof") ;

  c->Update() ;

  c->Print(TString(ws.cat("decaymode")->getLabel()) + "Fit.eps") ;

  //   std::ostringstream fileout;
  //   fileout<<"summary_"<<Bparticle<<".eps";
  //   c->Print(fileout.str().c_str());
  
}

void makeMassSPlots( DecayMode mode, bool refit = false )
{
  RooWorkspace* ws = gMyWorkspace ;
  if(ws == 0) ws = createWS(mode) ;
  
  // let's try to make an splot in the time. we need a new pdf first:
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws->function("pdf")) ;
  RooAbsPdf* masspdf = dynamic_cast<RooAbsPdf*>(ws->function("masspdf")) ;
  RooAbsPdf* t_prompt_pdf = dynamic_cast<RooAbsPdf*>(ws->function("t_prompt")) ;
  RooAbsPdf* t_sig_pdf = dynamic_cast<RooAbsPdf*>(ws->function("t_sig")) ;
  RooAbsPdf* t_bkg_pdf = dynamic_cast<RooAbsPdf*>(ws->function("t_bkg")) ;
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data")) ;
  RooRealVar* t = ws->var("t") ;
  RooRealVar* sigmat = ws->var("sigmat") ;
  sigmat->setBins(20) ;
  t->setBins(1100) ;
  
  if(masspdf==0 || data==0) return ;

  // 
  if(refit) 
    masspdf->fitTo(*data, NumCPU(4),Extended(true),Minos(false)) ;
  TCanvas* c2 = new TCanvas("c2","",700,300) ;
  c2->Divide(2,1) ;
  c2->cd(1) ;
  RooPlot* mplot = ws->var("m")->frame() ;
  data->plotOn(mplot) ;
  masspdf->plotOn(mplot) ;
  mplot->Draw() ;
  c2->cd(2) ;
  RooPlot* mpsiplot = ws->var("mdau1")->frame() ;
  data->plotOn(mpsiplot) ;
  masspdf->plotOn(mpsiplot) ;
  mpsiplot->Draw() ;
  c2->Update() ;

  //RooMsgService::instance().addStream(DEBUG,Topic(Tracing)) ; 
  RooMsgService::instance().addStream(INFO,Topic(Tracing)) ; 
  RooMsgService::instance().addStream(INFO,Topic(Contents)) ; 
  RooMsgService::instance().addStream(INFO,Topic(InputArguments)) ; 
  RooMsgService::instance().setGlobalKillBelow(DEBUG) ;
  //RooMsgService::instance().addStream(INFO,Topic(Integration)) ;

  // create an splot
  RooStats::SPlot* splot = 
    new RooStats::SPlot("splotdata","splotdata",*data,masspdf,*(ws->set("yields"))) ;

  std::cout << *(ws->set("yields")) << std::endl ;

  RooDataSet* wdata = splot->GetSDataSet() ;
  wdata->tree()->Print() ;
  //return ;

  TCanvas* c = new TCanvas("c1","",900,300) ;
  c->Divide(3,1) ;

  // release the variables that we want to fit 
  setConstant(*ws,"t_*_fll",false) ;
  setConstant(*ws,"t_*_fml",false) ;
  setConstant(*ws,"t_*_ml_tau",false) ;
  setConstant(*ws,"tres_*_*",false) ;
  setConstant(*ws,"tres_*_s3",true) ;
  pdf->getParameters( *data)->writeToFile("test.txt") ;

  RooDataSet * tmpdata ;
  RooDataHist* histdata ;
  
  c->cd(2) ;
  gPad->SetLogy() ;
  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nprompt_sw") ;
  histdata  = new RooDataHist("wtsigmat","",RooArgSet(*t,*sigmat),*tmpdata) ;
  sigmat->setVal(0.05) ;
  t_prompt_pdf->fitTo( *histdata, NumCPU(4),Minos(false), SumW2Error(kTRUE)) ;
  t_prompt_pdf->getParameters( *histdata)->writeToFile("t_prompt_pdf.txt") ;
  RooPlot* psiwtframe = t->frame(Bins(220),Range(-2,10));
  tmpdata->plotOn( psiwtframe,DataError(RooAbsData::SumW2)) ;
  t_prompt_pdf->plotOn( psiwtframe, ProjWData(*histdata) ) ;
  //pdf->plotOn(psiwtframe,ProjWData(*st_data),Components("psibkg*"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  psiwtframe->SetTitle("propertime splot for prompt psi") ;
  psiwtframe->SetMinimum(0.1) ;
  psiwtframe->Draw() ;

  c->cd(3) ;
  gPad->SetLogy() ;
  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nbkg_sw") ;
  histdata  = new RooDataHist("wtsigmat","",RooArgSet(*t,*sigmat),*tmpdata) ;
  sigmat->setVal(0.05) ;
  t_bkg_pdf->fitTo( *histdata, NumCPU(4),Minos(false), SumW2Error(kTRUE)) ;
  t_bkg_pdf->getParameters( *histdata)->writeToFile("t_bkg_pdf.txt") ;
  RooPlot* bkgwtframe = t->frame(Bins(220),Range(-2,10));
  tmpdata->plotOn( bkgwtframe,DataError(RooAbsData::SumW2)) ;
  t_bkg_pdf->plotOn( bkgwtframe, ProjWData(*histdata) ) ;
  //pdf->plotOn(bkgwtframe,ProjWData(*st_data),Components("psibkg*"),LineColor(psibkgcolor),LineStyle(kDashed),lw);
  bkgwtframe->SetTitle("propertime splot for background") ;
  bkgwtframe->SetMinimum(0.1) ;
  bkgwtframe->Draw() ;

  c->cd(1) ; 
  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nsig_sw") ;
  RooPlot* sigwtframe = t->frame(Bins(50),Range(-2,10));
  histdata = new RooDataHist("wtsigmat","",RooArgSet(*t,*sigmat),*tmpdata) ;
  //pdf->plotOn(sigwtframe,Components("sig*"),LineColor(sigcolor),LineStyle(kDashed),lw);
  //pdf->plotOn(sigwtframe,ProjWData(*st_data),LineStyle(kDashed),lw) ;
  //pdf->plotOn(sigwtframe,ProjWData(*st_data),Components("sig*"),LineColor(sigcolor),LineStyle(kDashed),lw);
  tmpdata->plotOn( sigwtframe,DataError(RooAbsData::SumW2) ) ;
  t_sig_pdf->plotOn(sigwtframe, ProjWData(*histdata) ) ;
  sigwtframe->SetTitle("propertime splot for signal") ;
  sigwtframe->Draw() ;
  
  c->Print(TString(ws->cat("decaymode")->getLabel()) + "lifetimesplots.eps") ;

  setConstant(*ws,"t_*",true) ;
  setConstant(*ws,"t_sig_tau",false) ;
  setConstant(*ws,"tres_*",true) ;
  writeParameters(*ws) ;

  /*
  if(result) {
    //RooRealVar *s_ms = data->meanVar(st); // ,0,"signalRegion");
    RooRealVar *s_sig_s1 = (RooRealVar*)result->floatParsFinal().find("sig_res_s1");
    RooRealVar *s_sig_s2 = (RooRealVar*)result->floatParsFinal().find("sig_res_s2");
    RooRealVar *s_bkg_s1 = (RooRealVar*)result->floatParsFinal().find("bkg_res_s1");
    RooRealVar *s_bkg_s2 = (RooRealVar*)result->floatParsFinal().find("bkg_res_s2");
    RooRealVar *s_sig_f2 = (RooRealVar*)result->floatParsFinal().find("sig_res_f2");
    RooRealVar *s_bkg_f2 = (RooRealVar*)result->floatParsFinal().find("bkg_res_f2");
    
    std::cout << s_sig_s1 << " " << s_sig_s2 << " " << s_sig_f2 << std::endl ;
    std::cout << s_bkg_s1 << " " << s_bkg_s2 << " " << s_bkg_f2 << std::endl ;
    
    std::cout << "s_sig_s1: " << s_sig_s1->getVal() << "+/-" << s_sig_s1->getError() << std::endl ;
    std::cout << "s_sig_s2: " << s_sig_s2->getVal() << "+/-" << s_sig_s2->getError() << std::endl ;
    std::cout << "s_sig_f2: " << s_sig_f2->getVal() << "+/-" << s_sig_f2->getError() << std::endl ;
    
    std::cout << "s_bkg_s1: " << s_bkg_s1->getVal() << "+/-" << s_bkg_s1->getError() << std::endl ;
    std::cout << "s_bkg_s2: " << s_bkg_s2->getVal() << "+/-" << s_bkg_s2->getError() << std::endl ;
    std::cout << "s_bkg_f2: " << s_bkg_f2->getVal() << "+/-" << s_bkg_f2->getError() << std::endl ;


    double avsigt = data->moment(st,1)  ;
    double ressig = avsigt * sqrt( (1-s_sig_f2->getVal()) * sqr(s_sig_s1->getVal()) +
				   s_sig_f2->getVal() * sqr(s_sig_s2->getVal() ) ) ;
    
    double resbkg = avsigt * sqrt( (1-s_bkg_f2->getVal()) * sqr(s_bkg_s1->getVal()) +
				   s_bkg_f2->getVal() * sqr(s_bkg_s2->getVal() ) ) ;
    
    
    std::cout << "Mean propertime error: " << data->moment(st,1) << std::endl ;
    std::cout << "Mean propertime error x scale factor: "     << ressig << std::endl ;
    std::cout << "Mean propertime error x scale factor bkg: " << resbkg << std::endl ;
  }

*/
}

void makeSPlots( DecayMode mode )
{
  RooWorkspace* ws = gMyWorkspace ;
  if(ws == 0) ws = createWS(mode) ;
  
  // let's try to make an splot in the time. we need a new pdf first:
  RooAbsPdf* pdf = dynamic_cast<RooAbsPdf*>(ws->function("pdf")) ;
  RooDataSet* data = dynamic_cast<RooDataSet*>(ws->data("data")) ;

  RooRealVar* pt = ws->var("pt") ;
  RooRealVar* mdau2 = ws->var("mdau2") ;
  RooRealVar* chi2 = ws->var("chi2dof") ;

  if(pdf==0 || data==0) return ;

  //RooMsgService::instance().addStream(DEBUG,Topic(Tracing)) ; 
  RooMsgService::instance().addStream(INFO,Topic(Tracing)) ; 
  RooMsgService::instance().addStream(INFO,Topic(Contents)) ; 
  RooMsgService::instance().addStream(INFO,Topic(InputArguments)) ; 
  RooMsgService::instance().setGlobalKillBelow(DEBUG) ;
  //RooMsgService::instance().addStream(INFO,Topic(Integration)) ;

  // create an splot
  RooStats::SPlot* splot = 
    new RooStats::SPlot("splotdata","splotdata",*data,pdf,*(ws->set("yields"))) ;

  RooDataSet* wdata = splot->GetSDataSet() ;
  wdata->tree()->Print() ;
  //return ;

  TCanvas* c1 = new TCanvas("c1","",900,900) ;
  c1->Divide(3,3) ;

  // release the variables that we want to fit 
  setConstant(*ws,"t_*_fll",false) ;
  setConstant(*ws,"t_*_fml",false) ;
  setConstant(*ws,"t_*_ml_tau",false) ;
  setConstant(*ws,"tres_*_*",false) ;
  setConstant(*ws,"tres_*_s3",true) ;
  pdf->getParameters( *data)->writeToFile("test.txt") ;

  RooDataSet * tmpdata ;
  RooPlot* ptframe ;
  RooPlot* chi2frame ;
  RooPlot* mdau2frame ;

  double mmin = 0 ;
  double mmax = 2000 ;
  switch( ws->cat("decaymode")->getIndex() ) {
  case JPsiPhi:
    mmin = 1000 ; mmax = 1039 ; break ;
  case JPsiKstar:
    mmin = 820 ; mmax = 980 ; break ;
  case JPsiKs:
    mmin = 465 ; mmax = 530 ; break ;
  default: break ;
  }
  
  // Prompt  
  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nsig_sw") ;
  c1->cd(1) ;
  mdau2frame = mdau2->frame(Bins(50),Range(mmin,mmax)) ;
  tmpdata->plotOn( mdau2frame,DataError(RooAbsData::SumW2) ) ;
  mdau2frame->SetTitle("mdau2 splot for signal") ;
  mdau2frame->Draw() ;


  c1->cd(4) ;
  ptframe = pt->frame(Bins(50),Range(0.,20000.));
  tmpdata->plotOn( ptframe,DataError(RooAbsData::SumW2) ) ;
  ptframe->SetTitle("pt(B) splot for signal") ;
  ptframe->Draw() ;

  c1->cd(7) ;
  chi2frame = chi2->frame(Bins(20));
  tmpdata->plotOn( chi2frame,DataError(RooAbsData::SumW2) ) ;
  chi2frame->SetTitle("chi2/dof splot for signal") ;
  chi2frame->Draw() ;
  delete tmpdata ;

  // Prompt psi
  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nprompt_sw") ;
  c1->cd(2) ;
  mdau2frame = mdau2->frame(Bins(50),Range(mmin,mmax)) ;
  tmpdata->plotOn( mdau2frame,DataError(RooAbsData::SumW2) ) ;
  mdau2frame->SetTitle("mdau2 splot for prompt psi") ;
  mdau2frame->Draw() ;

  c1->cd(5) ;
  ptframe = pt->frame(Bins(50),Range(0.,20000.));
  tmpdata->plotOn( ptframe,DataError(RooAbsData::SumW2) ) ;
  ptframe->SetTitle("pt(B) splot for prompt psi") ;
  ptframe->Draw() ;

  c1->cd(8) ;
  chi2frame = chi2->frame(Bins(20));
  tmpdata->plotOn( chi2frame,DataError(RooAbsData::SumW2) ) ;
  chi2frame->SetTitle("chi2/dof splot for prompt psi") ;
  chi2frame->Draw() ;
  delete tmpdata ;

  // Background

  tmpdata = new RooDataSet(wdata->GetName(),wdata->GetTitle(),wdata,*wdata->get(),0,"Nbkg_sw") ;
  c1->cd(3) ;
  mdau2frame = mdau2->frame(Bins(50),Range(mmin,mmax)) ;
  tmpdata->plotOn( mdau2frame,DataError(RooAbsData::SumW2) ) ;
  mdau2frame->SetTitle("mdau2 splot for background") ;
  mdau2frame->Draw() ;

  c1->cd(6) ;
  ptframe = pt->frame(Bins(50),Range(0.,20000.));
  tmpdata->plotOn( ptframe,DataError(RooAbsData::SumW2) ) ;
  ptframe->SetTitle("pt(B) splot for background") ;
  ptframe->Draw() ;

  c1->cd(9) ;
  chi2frame = chi2->frame(Bins(20));
  tmpdata->plotOn( chi2frame,DataError(RooAbsData::SumW2) ) ;
  chi2frame->SetTitle("chi2/dof splot for background") ;
  chi2frame->Draw() ;
  delete tmpdata ;

  c1->Print(TString(ws->cat("decaymode")->getLabel()) + "splots.eps") ;
}
