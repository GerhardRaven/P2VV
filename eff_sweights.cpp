#include <TCanvas.h>

#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include "RooStats/SPlot.h"

#include "eff_sweights.h"

using namespace RooFit;

RooDataSet* sPlotData( RooWorkspace& w, RooDataSet& data ) {

   // also build the mass-only pdf
   w.factory( "Gaussian::m_sig(m,m_sig_mean[5279,5250,5310],m_sig_sigma[10,5,20])" );
   w.factory( "Exponential::m_bkg(m,m_bkg_exp[-0.0002,-0.01,0.01])" );

   // J/psi mass pdf
   w.factory( "CBShape::mpsi_sig(mpsi,mpsi_sig_mean[3097,3085,3110],mpsi_sig_sigma[13.1,5,20],"
              "m_psi_sig_alpha[1.36,0.5,3],m_psi_sig_n[3,1,5])" );
   w.factory( "Exponential::mpsi_bkg(mpsi,mpsi_bkg_exp[-0.0002,-0.01,0.01])" );

   w.factory( "SUM::masspdf(n_sig[0,10000]*PROD( m_sig, mpsi_sig),"
              "n_prompt[0,1000000]*PROD( m_bkg, mpsi_sig ),n_bkg[0,1000000] "
              "* PROD( m_bkg, mpsi_bkg) )");
   w.defineSet( "yields","n_sig,n_prompt,n_bkg") ;

   // let's try to make an splot in the time. we need a new pdf first:
   RooAbsPdf* masspdf = dynamic_cast< RooAbsPdf* >( w.pdf( "masspdf" ) );

   masspdf->fitTo( data, NumCPU( 8 ), Extended( true ), Minos( false ) );
   TCanvas* c2 = new TCanvas( "mass_canvas", "mass_canvas", 1000, 500 );
   c2->Divide( 2, 1 );
   c2->cd( 1 );
   RooPlot* mplot = w.var( "m" )->frame();
   data.plotOn(mplot);
   masspdf->plotOn(mplot);
   mplot->Draw();

   c2->cd( 2 ) ;
   mplot = w.var( "mpsi" )->frame();
   data.plotOn( mplot );
   masspdf->plotOn( mplot );
   mplot->Draw();

   //RooMsgService::instance().addStream(DEBUG,Topic(Tracing)) ; 
   RooMsgService::instance().addStream(INFO,Topic(Tracing)) ; 
   RooMsgService::instance().addStream(INFO,Topic(Contents)) ; 
   RooMsgService::instance().addStream(INFO,Topic(InputArguments)) ; 
   RooMsgService::instance().setGlobalKillBelow(DEBUG) ;

   // create an splot
   RooStats::SPlot* splot = 
      new RooStats::SPlot( "splotdata", "splotdata", data, masspdf, *( w.set( "yields" ) ) );

   std::cout << *( w.set( "yields" ) ) << std::endl;

   RooDataSet* wdata = splot->GetSDataSet();

   return new RooDataSet( "sig_data", "sig_data", wdata, *(data.get()), "", "n_sig_sw" );
}
