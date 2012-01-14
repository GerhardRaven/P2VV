#include <TH1.h>
#include <TSystem.h>
#include <TCanvas.h>

#include <RooGlobalFunc.h>
#include <RooWorkspace.h>
#include <RooDataHist.h>
#include <RooRealVar.h>
#include <RooStudyManager.h>
#include <RooGenFitStudy.h>
#include <RooDataSet.h>
#include <RooPullVar.h>
#include <RooPlot.h>
#include <RooPullVar.h>

using namespace RooFit;

RooWorkspace* toyMC() {
   gSystem->Load("libP2VV");
   RooWorkspace* w = new RooWorkspace("w");

   w->factory("t[0,10]");

   w->factory("GaussModel::tres(t,tres_m[0],tres_s[0.0005])");

   w->factory("Decay::pdf(t,tau[2,0.01,4.0],tres,SingleSided)");
   w->defineSet("observables","t");

   unsigned int nbins = 10;
   TH1F* effh1 = new TH1F( "effh1", "effh1", nbins, w->var("t")->getMin(), w->var("t")->getMax());
   for (unsigned int i = 1; i < 6; ++i) {
      effh1->SetBinContent(i, 0.1 * i);
   }
   for (unsigned int i = 6; i < nbins + 1; ++i) {
      effh1->SetBinContent(i, 1);
   }
   RooDataHist effdatahist("effdatahist", "effdatahist", RooArgList(*(w->var("t"))), effh1);

   w->import(effdatahist);

   w->factory("HistFunc::eff(t, effdatahist)");
   w->factory("EffHistProd::acc_pdf(pdf, eff)");
   w->addClassDeclImportDir("..");
   w->addClassImplImportDir("..");
   w->importClassCode();

   RooAbsPdf* acc_pdf = w->pdf("acc_pdf");

   RooDataSet* data = acc_pdf->generate(*(w->set("observables")), 10000);

   TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
   canvas->cd(1);
   RooPlot* frame = w->var("t")->frame(RooFit::Bins(30));
   data->plotOn(frame, RooFit::MarkerSize(0.5), RooFit::XErrorSize(0));
   acc_pdf->plotOn(frame, RooFit::LineWidth(2));
   frame->Draw();

   // RooGenFitStudy gfs;
   // gfs.setGenConfig("acc_pdf", "t", NumEvents(10000));
   // gfs.setFitConfig("acc_pdf", "t", PrintLevel(-1));

   // RooStudyManager mgr(*w, gfs);

   // mgr.run(1000);
   // mgr.prepareBatchInput("aap", 100, kTRUE);

   // RooDataSet* data = gfs.summaryData();
   // const RooArgSet* args = data->get();
   // RooPullVar tau_pull("tau_pull", "tau_pull", *(static_cast<RooRealVar*>(args->find("tau"))), *(w->var("t")));
   // data->addColumn(tau_pull);

   // TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 1000);
   // canvas->Divide(2, 2);
   // canvas->cd(1);

   return w;
}
