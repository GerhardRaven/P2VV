from ROOT import *
gSystem.Load("libP2VV")
from math import sqrt,pi

from RooFitDecorators import *
from ModelBuilders import *
from ModelBuilders import _buildAngularFunction

ws = RooWorkspace("ws")

##################################
### Define variables and PDF's ###
##################################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,14.], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0],m[5200,5550]}"%(-pi,pi))

ws.defineSet("transversityangles","trcospsi,trcostheta,trphi")
angles = ws.set('transversityangles')

ab = abasis(ws,angles)

ws.factory("{rz2[0.601,0.4,0.7],rperp2[0.16,0.1,0.5],rs2[0.00001,0.,0.15]}")
ws.factory("RooFormulaVar::rpar2('1-@0-@1-@2',{rz2,rperp2,rs2})")
ws.factory("RooFormulaVar::rz('sqrt(@0)',{rz2})")
ws.factory("RooFormulaVar::rperp('sqrt(@0)',{rperp2})")
ws.factory("RooFormulaVar::rpar('sqrt(@0)',{rpar2})")
ws.factory("RooFormulaVar::rs('sqrt(@0)',{rs2})")

ws.factory("{deltaz[0.],deltapar[2.5,%f,%f],deltaperp[-0.17,%f,%f],deltas[0.5,%f,%f]}"%(-2*pi,2*pi,-2*pi,2*pi,-2*pi,2*pi))

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

ws.factory("{#Gamma[0.68,0.4,0.9],t_sig_dG[0.060,-2,2],t_sig_dm[17.8,15,20], phis[-0.04,%f,%f],tagomega[0.33,0.,0.501]}"%(-2*pi,2*pi))
ws.factory("expr::t_sig_tau('1/@0',{#Gamma})")
ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")
ws.factory("expr::wtag('tagomega',tagomega)")

ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.0005])")

###########
### SIG ###
###########
t_ang_pdf = buildJpsiphiSWave(ws,'t_ang_sig', True,'tres')

ws.factory("Gaussian::m_sig(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")
ws.factory("PROD::sigpdf(t_ang_sig, m_sig)")
ws.defineSet("observables","t,trcospsi,trcostheta,trphi,m")
sigpdf = ws['sigpdf']

###########
### BKG ###
###########
# Bkg mass
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

# Bkg time 
#Single exponential for biased sample!!!!
ws.factory("RooDecay::t_bkg(t,t_bkg_tau[0.2,0.01,2.0],tres,SingleSided)")

# Bkg angles
#Build angular background with RooP2VVAngleBasis
_ba = lambda name,comp : _buildAngularFunction(ws,ab,name,comp)

_ba("Bkg_0000",  [ ( 0,0,0,0,1.) ] )
c0000 = RooRealVar('c0000','c0000',1.)
#Flat background in all angles
ang_bkg = RooRealSumPdf('ang_bkg','ang_bkg',RooArgList(ws['Bkg_0000_basis']),RooArgList(c0000))
ws.put(ang_bkg)

#ws.factory("PROD::bkgpdf_tm( m_bkg, t_bkg)")
#ws.factory("PROD::bkgpdf( bkgpdf_tm, ang_bkg)")
ws.factory("PROD::bkgpdf( m_bkg, t_bkg, ang_bkg)")

sigpdf = ws['sigpdf']
bkgpdf = ws['bkgpdf']
#bkgpdf_tm = ws['bkgpdf_tm']

ws.factory("SUM::pdf(f_sig[0.71,0.,1.0]*sigpdf,bkgpdf)")

pdf = ws['pdf']
#pdf = ws['sigpdf']
#pdf = ws['bkgpdf']

#####################
### Generate data ###
#####################
print 'BEFORE GENERATING DATA!!!!!!!!!!!!!!!!!!!!!!!!'
data = pdf.generate(ws.set('observables'),10000)
print 'AFTER GENERATING DATA!!!!!!!!!!!!!!!!!!!!!!!!'

##############################################
### Define acceptance function a la Wouter ###
##############################################
ws.factory("expr::effshape('1/(1+(a*t)**(-c))',t,a[1.45],c[2.37])")
effhist = ws['effshape'].createHistogram('effhist',ws['t'],RooFit.Binning(10,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")

ws.factory("EffHistProd::accpdf(pdf,effpdf)")

accpdf = ws['accpdf']

#####################
### Generate data ###
#####################
print 'BEFORE GENERATING BIASED DATA!!!!!!!!!!!!!!!!!!!!!!!!'
accdata = accpdf.generate(ws.set('observables'),10000)
print 'AFTER GENERATING BIASED DATA!!!!!!!!!!!!!!!!!!!!!!!!'

accpdf.fitTo(accdata)
############
### Plot ###
############
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

#Make time plots
Acc = TCanvas('Acc','Acc')
Acc.Divide(2)

Acc.cd(1)
tframe = ws['t'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(tframe)
tframe.Draw()

Acc.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
ws['effpdf'].plotOn(mframe)
mframe.Draw()

timecanvas = TCanvas('timecanvas','timecanvas')
timecanvas.Divide(3,2)

timecanvas.cd(1)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

timecanvas.cd(2)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
pdf.plotOn(tframe,lw)
tframe.Draw()

timecanvas.cd(3)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(tframe,lw)
tframe.Draw()

timecanvas.cd(4)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
accdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
tframe.Draw()

timecanvas.cd(5)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
#accpdf.plotOn(tframe,lw,RooFit.ProjWData(RooArgSet(ws['m']),accdata))
accpdf.plotOn(tframe,lw,RooFit.Project(RooArgSet(ws['m'])))
#accpdf.plotOn(tframe,lw)
tframe.Draw()

timecanvas.cd(6)
#gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
accdata.plotOn(tframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(tframe,lw)
tframe.Draw()

#Make trcospsi plots
trcospsicanvas = TCanvas('trcospsicanvas','trcospsicanvas')
trcospsicanvas.Divide(3,2)

trcospsicanvas.cd(1)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
data.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
trcospsiframe.Draw()

trcospsicanvas.cd(2)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
pdf.plotOn(trcospsiframe,lw)
trcospsiframe.Draw()

trcospsicanvas.cd(3)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
data.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcospsiframe,lw)
trcospsiframe.Draw()

trcospsicanvas.cd(4)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
accdata.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
trcospsiframe.Draw()

trcospsicanvas.cd(5)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
accpdf.plotOn(trcospsiframe,lw)
trcospsiframe.Draw()

trcospsicanvas.cd(6)
trcospsiframe = ws['trcospsi'].frame(RooFit.Bins(10))
accdata.plotOn(trcospsiframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(trcospsiframe,lw)
trcospsiframe.Draw()

#Make trcostheta plots
trcosthetacanvas = TCanvas('trcosthetacanvas','trcosthetacanvas')
trcosthetacanvas.Divide(3,2)

trcosthetacanvas.cd(1)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
data.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
trcosthetaframe.Draw()

trcosthetacanvas.cd(2)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
pdf.plotOn(trcosthetaframe,lw)
trcosthetaframe.Draw()

trcosthetacanvas.cd(3)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
data.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trcosthetaframe,lw)
trcosthetaframe.Draw()

trcosthetacanvas.cd(4)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
accdata.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
trcosthetaframe.Draw()

trcosthetacanvas.cd(5)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
accpdf.plotOn(trcosthetaframe,lw)
trcosthetaframe.Draw()

trcosthetacanvas.cd(6)
trcosthetaframe = ws['trcostheta'].frame(RooFit.Bins(10))
accdata.plotOn(trcosthetaframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(trcosthetaframe,lw)
trcosthetaframe.Draw()

#Make trphi plots
trphicanvas = TCanvas('trphicanvas','trphicanvas')
trphicanvas.Divide(3,2)

trphicanvas.cd(1)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
data.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
trphiframe.Draw()

trphicanvas.cd(2)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
pdf.plotOn(trphiframe,lw)
trphiframe.Draw()

trphicanvas.cd(3)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
data.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(trphiframe,lw)
trphiframe.Draw()

trphicanvas.cd(4)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
accdata.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
trphiframe.Draw()

trphicanvas.cd(5)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
accpdf.plotOn(trphiframe,lw)
trphiframe.Draw()

trphicanvas.cd(6)
trphiframe = ws['trphi'].frame(RooFit.Bins(10))
accdata.plotOn(trphiframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(trphiframe,lw)
trphiframe.Draw()

#Make mass plots
masscanvas = TCanvas('masscanvas','masscanvas')
masscanvas.Divide(3,2)

masscanvas.cd(1)
mframe = ws['m'].frame(RooFit.Bins(10))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

masscanvas.cd(2)
mframe = ws['m'].frame(RooFit.Bins(10))
pdf.plotOn(mframe,lw)
mframe.Draw()

masscanvas.cd(3)
mframe = ws['m'].frame(RooFit.Bins(10))
data.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
pdf.plotOn(mframe,lw)
mframe.Draw()

masscanvas.cd(4)
mframe = ws['m'].frame(RooFit.Bins(10))
accdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
mframe.Draw()

masscanvas.cd(5)
mframe = ws['m'].frame(RooFit.Bins(10))
accpdf.plotOn(mframe,lw)
mframe.Draw()

masscanvas.cd(6)
mframe = ws['m'].frame(RooFit.Bins(10))
accdata.plotOn(mframe,RooFit.MarkerSize(0.5),xes)
accpdf.plotOn(mframe,lw)
mframe.Draw()
