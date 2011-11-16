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

############################
### Flags and rootfiles  ###
############################

ws = RooWorkspace("ws")

###################
### Observables ###
###################

ws.factory("{ trcospsi[-1,1], trcostheta[-1,1], trphi[%f,%f], t[0.3,14], tagdecision[Bs_Jpsiphi=+1,Bsbar_Jpsiphi=-1,untagged=0]}"%(-pi,pi))
ws.factory("{m[5200,5550],mdau1[3030,3150],mdau2[1007.46,1031.46]}")
ws.factory("{biased[Biased=+1,NotBiased=0],unbiased[Unbiased=+1,NotUnbiased=0]}")

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

#Parametrize in terms of S,D,C:
#ws.factory("{S[-0.04,-2.,2.],D[1.,-2.,2.],C[0]}")

#Parametrize in terms of phis:
ws.factory('{phis[-0.04,%f,%f]}'%(-2*pi,2*pi))
ws.factory("{expr::S('-1*sin(@0)',{phis}),expr::D('cos(@0)',{phis}),C[0]}")

#Parametrize in terms of lambda
#ws.factory('{phis[-0.04,%f,%f]}'%(-2*pi,2*pi))
#ws.factory("{relambda[1.,-10,10],imlambda[-0.04,-10,10]}")
#ws.factory("{expr::relambda('cos(@0)',{phis}),expr::imlambda('-1*sin(@0)',{phis})}")
#ws.factory("{expr::lambda2('@0*@0+@1*@1',{relambda,imlambda}),expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),expr::C('(1-@0)/(1+@0)',{lambda2})}")
#ws.factory("{expr::lambda2('@0*@0+@1*@1',{relambda,imlambda}),expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),C[0]}")
#ws.factory("{lambda2[1],expr::S('(2*@0)/(1+@1)',{imlambda,lambda2}),expr::D('(2*@0)/(1+@1)',{relambda,lambda2}),expr::C('(1-@0)/(1+@0)',{lambda2})}")

###############################
### Experimental parameters ###
###############################

#2011
ws.factory("GaussModel::tres_3(t,tres_mu[-0.0027],tres_s3[0.513],tres_SF[1.00,0.5,1.5])")
ws.factory("GaussModel::tres_2(t,tres_mu,tres_s2[0.0853],tres_SF)")
ws.factory("GaussModel::tres_1(t,tres_mu,tres_s1[0.0434],tres_SF)")
ws.factory("AddModel::tres({tres_3,tres_2,tres_1},{tres_f3[0.0017],tres_f2[0.165]})")

#ws.factory("GaussModel::tres(t,tres_mean[0.0],tres_sigma[0.05])")

# For determination of efficiency from MC, put wtag to 0, later integrate out (or set to 0.5 as I used to do before). Does that imply rebuilding the JpsiPhi pdf?
#ws.factory("tagomega[0.,0.,0.5]")
ws.factory("tagomega[0.33]")
ws.factory("expr::wtag('tagomega',tagomega)")

#Build the signal PDF
newpdf = buildJpsiphiSWave(ws,'newpdf', True,'tres')

# Make SuperCategory from (triggeredByUnbiasedHlt1AndHlt2,triggeredByBiasedHlt1AndHlt2)
TypeCat = RooSuperCategory('TypeCat','TypeCat',RooArgSet(ws['unbiased'],ws['biased']))
ws.put(TypeCat)
# Make MappedCategory from SuperCategory to split in unbiased and fullybiased
#fitcat = RooMappedCategory('fitcat','fitcat',ws['TypeCat'],'00')#'00' means NotUnbiased && NotBiased
fitcat = RooMappedCategory('fitcat','fitcat',ws['TypeCat'],'AllUnbiased')#'00' means NotUnbiased && NotBiased
fitcat.map("{Unbiased;NotBiased}","AllUnbiased")
fitcat.map("{Unbiased;Biased}","AllUnbiased") 
fitcat.map("{NotUnbiased;Biased}","FullyBiased")
fitcat = fitcat.createFundamental()
ws.put(fitcat)

##############################
### Proper Time Acceptance ###
##############################
ws.factory("expr::effshape('1/(1+(a*t)**(-c))',t,a[1.45],c[2.37])")
#ws.factory("expr::effshape('(1+b*t)/(1+(a*t)**(-c))',t,a[1.45],b[-0.0157],c[2.37])")
effhist = ws['effshape'].createHistogram('effhist',ws['t'],RooFit.Binning(10,ws['t'].getMin(),ws['t'].getMax()))
effdatahist = RooDataHist("effdatahist","effdatahist",RooArgList(ws['t']),effhist)
ws.put(effdatahist)
ws.factory("HistPdf::effpdf(t,effdatahist)")

#######################
### Build the PDF's ###
#######################
# Signal mass
ws.factory("Gaussian::m_sig_1(m,m_sig_mean[5365,5360,5370],m_sig_sigma_1[6.,0.,20.])")
ws.factory("expr::m_sig_sigma_2('2.14*@0',{m_sig_sigma_1})")
ws.factory("Gaussian::m_sig_2(m,m_sig_mean,m_sig_sigma_2)")
ws.factory("SUM::m_sig(m_bkg_f[0.83]*m_sig_1,m_sig_2)")

# Bkg mass
ws.factory("Exponential::m_bkg(m,m_bkg_exp[-0.001,-0.01,-0.0001])")

# Bkg time
#Double exponential for unbiased sample
ws.factory("RooDecay::ml(t,t_bkg_ml_tau[0.21,0.01,0.5],tres,SingleSided)")
ws.factory("RooDecay::ll(t,t_bkg_ll_tau[1.92,0.5,2.5],tres,SingleSided)")
ws.factory("SUM::t_bkg_UB(t_bkg_fll[0.3,0.,1.]*ll,ml)")

#Single exponential for biased sample
ws.factory("RooDecay::t_bkg_B(t,t_bkg_B_tau[0.2,0.01,2.0],tres,SingleSided)")

ws.factory("Uniform::ang_bkg({trcostheta,trcospsi,trphi})")

###########################
### Putting it together ###
###########################
# Mass only PDF
ws.factory("SUM::m_pdf(Nsig_all[1000,0,16000]*m_sig,Nbkg_all[1000,0,16000]*m_bkg)")

# Sig PDF

#UB
# SigPDF(t,angles,m)
ws.factory("PROD::sig_pdf_UB( m_sig, newpdf)")

#B
# Biased SigPDF(t,angles)
ws.factory("EffHistProd::acc_newpdf(newpdf,effpdf)")
# Full Biased SigPDF(t,angles,m)
ws.factory("PROD::sig_pdf_B( m_sig, acc_newpdf)")

# Bkg PDF

#UB
ws.factory("PROD::bkg_pdf_UB( m_bkg, t_bkg_UB, ang_bkg)")

#B
# Biased BkgPDF(t)
ws.factory("EffHistProd::acc_t_bkg_B(t_bkg_B,effpdf)")
# Full BkgPDF(t,angles,m)
ws.factory("PROD::bkg_pdf_B( m_bkg,acc_t_bkg_B,ang_bkg)")

ws.factory("SUM::pdf_ext_UB(Nsig_UB[7288]*sig_pdf_UB,Nbkg_UB[3739]*bkg_pdf_UB)")
ws.factory("SUM::pdf_UB(f_sig_UB[0.71,0.,1.0]*sig_pdf_UB,bkg_pdf_UB)")

ws.factory("SUM::pdf_ext_B(Nsig_B[1186]*sig_pdf_B,Nbkg_B[568]*bkg_pdf_B)")
ws.factory("SUM::pdf_B(f_sig_B[0.71,0.,1.0]*sig_pdf_B,bkg_pdf_B)")

ws.factory("Simultaneous::simpdf(fitcat)")
ws['simpdf'].addPdf(ws['pdf_ext_UB'],'AllUnbiased')
ws['simpdf'].addPdf(ws['pdf_ext_B'],'FullyBiased')

wsfile = TFile('ToySimWS.root','RECREATE')
ws.Write()
wsfile.Close()

pdf =  ws['simpdf']
ws.defineSet("observables","t,m,fitcat")
ras = ws.set('observables')
data = pdf.generate(ras,0)

ws['rs2'].setVal(0.04)

#sw = TStopwatch()
#sw.Start()
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'),ws.pdf('tres_SFconstraint'))))
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True))
#sw.Stop()
#result.SaveAs('result_%s.root'%(name))
#result.writepars(name,False)
#result.writecorr(name)

######PLOT######
lw = RooCmdArg(RooFit.LineWidth(2))
xes = RooCmdArg(RooFit.XErrorSize(0))
err = RooCmdArg(RooFit.DrawOption('E'))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))

signame_UB = 'sig_pdf_UB'
bkgname_UB = 'bkg_pdf_UB'
signame_B = 'sig_pdf_B'
bkgname_B = 'bkg_pdf_B'

sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
nonpsicolor = RooCmdArg(RooFit.LineColor(RooFit.kOrange))

msigmin = 5345.
msigmax = 5387.
mmin = 5200.
mmax = 5550.
ws['m'].setRange('sigRegion',msigmin,msigmax)
ws['m'].setRange('leftSideband',mmin,msigmin)
ws['m'].setRange('rightSideband',msigmax,mmax)

projectdata = RooDataSet('projectdata','projectdata',data,RooArgSet(ws['fitcat']))

#Make Sanity plots
masscanvas = TCanvas('masscanvas','masscanvas')
masscanvas.Divide(1,2)

masscanvas.cd(1)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.Cut("fitcat==fitcat::AllUnbiased"))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(signame_UB),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(bkgname_UB),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
mframe.Draw()

masscanvas.cd(2)
mframe = ws['m'].frame(RooFit.Bins(30))
data.plotOn(mframe,RooFit.Cut("fitcat==fitcat::FullyBiased"))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(signame_B),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(mframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(bkgname_B),RooFit.ProjWData(projectdata),dashed,bkgcolor)
mframe.Draw()

masscanvas.SaveAs('plaatjemass.eps')

timecanvas = TCanvas('timecanvas','timecanvas')
timecanvas.Divide(1,2)

timecanvas.cd(1)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.Cut("fitcat==fitcat::AllUnbiased"))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(signame_UB),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"AllUnbiased"),RooFit.Components(bkgname_UB),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
tframe.Draw()

timecanvas.cd(2)
gPad.SetLogy()
tframe = ws['t'].frame(RooFit.Bins(30))
data.plotOn(tframe,RooFit.Cut("fitcat==fitcat::FullyBiased"))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.ProjWData(projectdata))
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(signame_B),RooFit.ProjWData(projectdata),dashed,sigcolor)
pdf.plotOn(tframe,RooFit.Slice(ws['fitcat'],"FullyBiased"),RooFit.Components(bkgname_B),RooFit.ProjWData(projectdata),dashed,bkgcolor) 
tframe.Draw()

timecanvas.SaveAs('plaatjetime.eps')
