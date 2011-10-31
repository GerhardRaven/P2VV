from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
import rootStyle
from ModelBuilders import _buildAngularFunction
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

from ModelBuilders import *

##############################
### Create ws, observables ###
##############################

ws = RooWorkspace('ws')

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,14], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0],m[5200,5550]}"%(-pi,pi))

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

ws.factory("{rz2[0.60,0.,1.],rperp2[0.16,0.,1.],rs2[0.,0.,1.]}")
ws.factory("RooFormulaVar::rpar2('1-@0-@1-@2',{rz2,rperp2,rs2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
ws.factory("RooFormulaVar::rs('sqrt(@0)',{rs2})")

ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f],deltas[0.,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi,-2*pi,2*pi))

#Parametrize differently
#ws.factory("{NAzAz[1],NAparApar[0.44,-1,1],NAperpAperp[0.4,-1,1],ReAparAperp[-0.114,-1,1],ReAzAperp[0.08,-1,1],ReAzApar[-0.66,-1,1],ImAparAperp[0.410,-1,1],ImAzAperp[-0.63,-1,1]}")

ws.factory("expr::ReAz   ('@0    * cos(@1)',   {rz,deltaz})")
ws.factory("expr::ImAz   ('@0    * sin(@1)',   {rz,deltaz})")
ws.factory("expr::ReApar ('@0    * cos(@1)', {rpar,deltapar})")
ws.factory("expr::ImApar ('@0    * sin(@1)', {rpar,deltapar})")
ws.factory("expr::ReAperp('@0    * cos(@1)',{rperp,deltaperp})")
ws.factory("expr::ImAperp('@0    * sin(@1)',{rperp,deltaperp})")
ws.factory("expr::ReAs('@0 * cos(@1)',{rs,deltas})")
ws.factory("expr::ImAs('@0 * sin(@1)',{rs,deltas})")

# define the relevant combinations of strong amplitudes
ws.factory("expr::NAzAz      ('( @0 * @0 + @1 * @1 )',{ReAz,    ImAz                      })")  # |A_z|^2
ws.factory("expr::NAparApar  ('( @0 * @0 + @1 * @1 )',{ReApar,  ImApar                    })")  # |A_par|^2
ws.factory("expr::NAperpAperp('( @0 * @0 + @1 * @1 )',{ReAperp, ImAperp                   })")  # |A_perp|^2
ws.factory("expr::ReAparAperp('( @0 * @2 + @1 * @3 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par||A_perp| cos(delta_perp - delta_par)
ws.factory("expr::ReAzAperp  ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   cos(delta_perp - delta_z)
ws.factory("expr::ReAzApar   ('( @0 * @2 + @1 * @3 )',{ReAz,    ImAz,    ReApar,  ImApar  })")  # |A_z||A_par|    cos(delta_par  - delta_z)
ws.factory("expr::ImAparAperp('( @0 * @3 - @1 * @2 )',{ReApar,  ImApar,  ReAperp, ImAperp })")  # |A_par|A_perp|  sin(delta_perp - delta_par)
ws.factory("expr::ImAzAperp  ('( @0 * @3 - @1 * @2 )',{ReAz,    ImAz,    ReAperp, ImAperp })")  # |A_z||A_perp|   sin(delta_perp - delta_z)

ws.factory("expr::NAsAs      ('( @0 * @0 + @1 * @1 )',{ReAs,    ImAs                      })")  # |A_s|^2
ws.factory("expr::ReAsAz     ('( @0 * @2 + @1 * @3 )',{ReAs,    ImAs,    ReAz,    ImAz    })")  # |A_s||A_z| cos(delta_z - delta_s)
ws.factory("expr::ReAsApar   ('( @0 * @2 + @1 * @3 )',{ReAs,    ImAs,    ReApar,  ImApar  })")  # |A_s||A_par| cos(delta_par - delta_s)
ws.factory("expr::ImAsAperp  ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReAperp, ImAperp })")  # |A_s||A_perp| sin(delta_perp - delta_s)
ws.factory("expr::ImAsAz     ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReAz,    ImAz    })")  # |A_s||A_z| sin(delta_z - delta_s)
ws.factory("expr::ImAsApar   ('( @0 * @3 - @1 * @2 )',{ReAs,    ImAs,    ReApar,  ImApar  })")  # |A_s||A_par| sin(delta_par - delta_s)

##########################
### physics parameters ###
##########################
ws.factory("{#Gamma[0.68,0.4,0.9]}")
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")

ws.factory("{t_sig_dG[0.060,-2,2]}")

ws.factory("{t_sig_dm[17.8,15,20]}")

#Parametrize in terms of phis:
ws.factory('{phis[-0.04,%f,%f]}'%(-2*pi,2*pi))
ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")

ws.factory("tagomega[0.,0.,0.5]")
ws.factory("expr::wtag('tagomega',tagomega)")

##############################
### Proper Time Acceptance ###
##############################
ws.factory("expr::effshape('1/(1+(a*t)**(-c))',t,a[1.45],c[2.37])")
#ws.factory("expr::effshape('(1+b*t)/(1+(a*t)**(-c))',t,a[1.45],b[-0.0157],c[2.37])")
effhist = ws['effshape'].createHistogram('effhist',ws['t'],RooFit.Binning(200,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")
ws.factory("HistPdf::effpdf2(t,effdatahist)")

#################
### PDF's     ###
#################

#Resolution Model
ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.05])")

#signal pdf
newpdf = buildJpsiphiSWave(ws,'newpdf', True,'tres')

# Signal mass
ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")

ws.factory("PROD::sig_pdf(m_sig,newpdf)")

#BKG mass
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

#BKG angles
ws.factory("Uniform::ang_bkg({trcostheta,trcospsi,trphi})")

#BKG time
#IT'S HERE!!!!!! Single no ERROR, Double does have ERROR!

#ws.factory("RooDecay::t_bkg(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")

ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
ws.factory("SUM::t_bkg(t_bkg_fll[0.3,0.,1.]*ll,ml)")

ws.factory("PROD::bkg_pdf(t_bkg,ang_bkg,m_bkg)")

#ws.factory("SUM::pdf_ext(Nsig[1186]*sig_pdf,Nbkg[568]*bkg_pdf)")

ws.factory("EffHistProd::acc_sig_pdf(sig_pdf,effpdf)")
ws.factory("EffHistProd::acc_bkg_pdf(bkg_pdf,effpdf)")
ws.factory("SUM::accpdf_ext( Nsig[1186]*acc_sig_pdf,Nbkg[568]*acc_bkg_pdf)")

#########################
### What do you want? ###
#########################
pdf = ws.pdf('accpdf_ext')
#pdf = ws.pdf('pdf_ext')
data = pdf.generate(RooArgSet(ws.var('t')),0)

signame = 'acc_sig_pdf'
bkgname = 'acc_bkg_pdf'

#signame = 'sig_pdf'
#bkgname = 'bkg_pdf'

###############################################################
### To see the fit problem for the bkg component of the pdf ###
###############################################################

lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))

sigcolor = RooFit.kGreen 
bkgcolor = RooFit.kRed

testCanvas = TCanvas('testCanvas','testCanvas') 
gPad.SetLogy()
_tb = ws.var('t').frame(RooFit.Bins(100))
data.plotOn(_tb,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(_tb,lw)
pdf.plotOn(_tb,RooFit.Components(signame),RooFit.LineColor(sigcolor),RooFit.LineStyle(kDashed),lw)
pdf.plotOn(_tb,RooFit.Components(bkgname),RooFit.LineColor(bkgcolor),RooFit.LineStyle(kDashed),lw)
#_tb.SetMinimum(0.1) 
#_tb.SetTitle("")
_tb.Draw()
testCanvas.Update()



################################################################

