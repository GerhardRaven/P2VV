from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *
import rootStyle
from ModelBuilders import _buildAngularFunction
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

##############
### Get WS ###
##############

wsfile = TFile('UntaggedWS.root')
ws = wsfile.Get('ws')

angcorrpdf = ws.pdf('pdf_ext_angcorrpdf')
data = ws.data('data')

################
### Blinding ###
################
blinded = False

if blinded:
   #Building blinded parameters
   ws.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
   ws.factory("expr::t_sig_tau_blind('1/@0',{#Gamma_unblind})")
   ws.factory("RooUnblindUniform::t_sig_dG_blind('BsHobbes',0.2,t_sig_dG)")

   customizer = RooCustomizer(angcorrpdf,'blinded')

   tau_unblind = ws.function('t_sig_tau')
   tau_blind = ws.function('t_sig_tau_blind')

   dG_unblind = ws.var('t_sig_dG')
   dG_blind = ws.function('t_sig_dG_blind')

   customizer.replaceArg( tau_unblind, tau_blind )
   customizer.replaceArg( dG_unblind, dG_blind )

   blindedpdf = customizer.build()
   ws.put(blindedpdf)

if blinded:
    pdf = ws.pdf('pdf_ext_angcorrpdf_blinded')
else:
    pdf = ws.pdf('pdf_ext_angcorrpdf')

#Proper time acceptance 
#ptacc = RooFormulaVar('ptacc','ptacc','1-0.025*@0',RooArgList(t))
#finalpdf = RooEffProd('finalpdf','finalpdf',pdf,ptacc)

###########
### Fit ###
###########

#Set back some values before the fit!

ws['wtag'].setVal(0.5)
ws['phis'].setVal(0)

result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

paramlist = pdf.getParameters(data)
writeFitParamsLatex(paramlist)

dict = writeCorrMatrixLatex(result)

################
### Profiles ###
################

gamma = ws.var('#Gamma')
deltaGamma = ws.var('t_sig_dG')
phis = ws.var('phis')

#setting back values
ws.var('#Gamma').setVal(0.68)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(0.0)

#MakeProfile('ProfiledGamma_Gamma',data,pdf,12,gamma,0.55,0.85,deltaGamma,-0.35,0.45)

#setting back values
ws.var('#Gamma').setVal(0.68)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(0.0)
ws.var('phis').setConstant(kFALSE)
#With phis unconstrained we now also have sensitivity to deltaperp!!!! 
ws.var('deltaperp').setMin(-2*pi)
ws.var('deltaperp').setMax(2*pi)
ws.var('deltaperp').setConstant(kFALSE)

#MakeProfile('ProfiledGamma_phis_untagged',data,pdf,15,phis,-pi,pi,deltaGamma,-0.7,0.7)


#We might still need this to see if the fits are fine actually, I remember seeing something fits hitting borders.....

## param1 = deltaGamma
## param2 = phis

## param1.setMin(-0.7)
## param1.setMax(0.7)
## param2.setMin(-pi)
## param2.setMax(pi)

## param1.setConstant(kFALSE)
## param2.setConstant(kFALSE)

## param1_min = param1.getMin()
## param1_max = param1.getMax()
## param1_int = param1_min-param1_max

## #get range for param2
## param2_min = param2.getMin()
## param2_max = param2.getMax()
## param2_int = param2_max - param2_min

## print '**************************************************'
## print 'Minimizing NLL'
## print '**************************************************'
## nll = pdf.createNLL(data,RooFit.NumCPU(8))
## m = RooMinuit(nll)
## #m.setVerbose(kTRUE)
## m.migrad()
## pdf.getParameters(data).Print("s")

## assert False

## for i in range(1,npoints+1):
##     param1_val = param1_min + (i-1)*(param1_int/npoints)
##     for j in range(1,npoints+1):
##         param2_val = param2_min + (j-1)*(param2_int/npoints)
##         print '***************************************************************************'
##         print 'At gridpoint i = %i from %i and j = %i from %i'%(i,npoints,j,npoints)
##         print '%s at current gridpoint ='%(param1.GetName()), param1_val
##         print '%s at current gridpoint ='%(param2.GetName()), param2_val
##         print '***************************************************************************'
##         param1.setVal(param1_val)
##         #param1.setConstant(kTRUE)
##         param2.setVal(param2_val)
##         #param2.setConstant(kTRUE)
##         #result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
##         #nll = pdf.createNLL(data,RooFit.NumCPU(8))
##         #nllval = nll.getVal()
##         #ProfileLikelihood.SetBinContent(i,j,nllval)
##         ProfileLikelihood.SetBinContent(i,j,prof.getVal())
##         #Heights.SetBinContent(i,j,2*(-1*(nllMINval)+nllval))
     
## gStyle.SetPalette(1)
## gStyle.SetOptStat(0)
## Canvas = TCanvas('Canvas','Canvas')

## ProfileLikelihood.Draw()

## tfile = TFile('profilelikelihood.root','RECREATE')
## ProfileLikelihood.Write()
## tfile.Close()

#################
### Single LL ###
#################
LL = False
if LL:
    CanL = TCanvas('CanL','CanL')
    CanL.Divide(2,1)

    #Construct unbinned likelihood
    nll = pdf.createNLL(data,RooFit.NumCPU(8)) 

    #Minimize likelihood w.r.t all parameters before making plots
    minuit = RooMinuit(nll)
    minuit.migrad()

    #Plot likelihood scan frac 
    CanL.cd(1)
    gammaframe = gamma.frame(RooFit.Bins(10),RooFit.Title("LL in #Gamma"))#Range(0.01,0.95)
    gammaframe.SetXTitle('#Gamma')
    nll.plotOn(gammaframe,RooFit.ShiftToZero())
    gammaframe.Draw()

    #Plot likelihood scan in sigma_g2
    CanL.cd(2)
    deltaGammaframe = deltaGamma.frame(RooFit.Bins(10),RooFit.Title("LL in #Delta #Gamma"))#Range(3.3,5.0),
    deltaGammaframe.SetXTitle('#Delta #Gamma')
    nll.plotOn(deltaGammaframe,RooFit.ShiftToZero())
    deltaGammaframe.Draw()
    
##################
###   Plotting ###
##################

# Watch out, the names of the signal and background components are hard-coded!!!!

plotting = False

if plotting:
    msigmin = 5345.
    msigmax = 5387.
    mmin = 5200.
    mmax = 5550
    ws.var('m').setRange('sigRegion',msigmin,msigmax)
    ws.var('m').setRange('leftSideband',mmin,msigmin)
    ws.var('m').setRange('rightSideband',msigmax,mmax)

    canvas = TCanvas('MassTime','MassTime')
    canvas.Divide(3,2)

    canvas.cd(1)
    myline1=TLine(tmin,msigmin,tmax,msigmin)
    myline1.SetLineColor(3)
    myline1.SetLineWidth(2)
    myline1.SetLineStyle(1)
    myline2=TLine(tmin,msigmax,tmax,msigmax)
    myline2.SetLineColor(3)
    myline2.SetLineWidth(2)
    myline2.SetLineStyle(1)

    hist = data.createHistogram(ws.var('t'),ws.var('m'))
    hist.SetMarkerSize(0.3)
    hist.SetMarkerStyle(20)
    hist.SetStats(kFALSE)
    hist.GetXaxis().SetTitle(str(ws.var('t').getTitle()))
    hist.GetYaxis().SetTitle(str(ws.var('m').getTitle()))
    hist.SetTitle('B_{s} mass vs. proper time')
    hist.Draw()
    myline1.Draw('same')
    myline2.Draw('same')
    canvas.Update()


    _c2 = plot( canvas.cd(2),ws.var('m'),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('B_{s} mass'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0) )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = False
                )

    _c3 = plot( canvas.cd(3),ws.var('t'),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('proper time full mass range'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0) )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = True
                )

    _c4 = plot( canvas.cd(4),ws.var('t'),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title("proper time signal region") )
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('sigRegion') )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = True
                )

    _c5 = plot( canvas.cd(5),ws.var('t'),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title("proper time signal region") )
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('leftSideband') )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('leftSideband')) 
                , logy = True
                )

    _c6 = plot( canvas.cd(6),ws.var('t'),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title("proper time signal region") )
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('rightSideband') )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('rightSideband')) 
                , logy = True
                )

    anglescanvas = TCanvas('Angles','Angles')
    anglescanvas.Divide(3,4)

    anglesnamelist = angles.nameList()

    _ac1 = plot( anglescanvas.cd(1),ws.var(anglesnamelist[0]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#psi)'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = False
                )

    _ac2 = plot( anglescanvas.cd(2),ws.var(anglesnamelist[1]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#theta)'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = False
                )

    _ac3 = plot( anglescanvas.cd(3),ws.var(anglesnamelist[2]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('#phi'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  )
                , pdfOpts = ( RooFit.LineWidth(2), ) 
                , logy = False
                )

    _ac4 = plot( anglescanvas.cd(4),ws.var(anglesnamelist[0]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#psi) in signal region'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('sigRegion')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('sigRegion')) 
                , logy = False
                )

    _ac5 = plot( anglescanvas.cd(5),ws.var(anglesnamelist[1]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#theta) in signal region'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('sigRegion')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('sigRegion')) 
                , logy = False
                )

    _ac6 = plot( anglescanvas.cd(6),ws.var(anglesnamelist[2]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('#phi in signal region'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  RooFit.CutRange('sigRegion'))
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('sigRegion')) 
                , logy = False
                )

    _ac7 = plot( anglescanvas.cd(7),ws.var(anglesnamelist[0]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#psi) in left sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('leftSideband')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('leftSideband')) 
                , logy = False
                )

    _ac8 = plot( anglescanvas.cd(8),ws.var(anglesnamelist[1]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#theta) in left sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('leftSideband')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('leftSideband')) 
                , logy = False
                )

    _ac9 = plot( anglescanvas.cd(9),ws.var(anglesnamelist[2]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('#phi in left sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  RooFit.CutRange('leftSideband'))
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('leftSideband')) 
                , logy = False
                )

    _ac10 = plot( anglescanvas.cd(10),ws.var(anglesnamelist[0]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#psi) in right sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('rightSideband')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('rightSideband')) 
                , logy = False
                )

    _ac11 = plot( anglescanvas.cd(11),ws.var(anglesnamelist[1]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('cos(#theta) in right sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0), RooFit.CutRange('rightSideband')  )
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('rightSideband')) 
                , logy = False
                )

    _ac12 = plot( anglescanvas.cd(12),ws.var(anglesnamelist[2]),data,pdf
                , { 'sig_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) )
                    , 'bkg_pdf_angcorrpdf'  : ( RooFit.LineColor(RooFit.kRed),RooFit.LineStyle(kDashed) )
                    }
                , frameOpts = ( RooFit.Bins(30), RooFit.Title('#phi in right sideband'))
                , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0),  RooFit.CutRange('rightSideband'))
                , pdfOpts = ( RooFit.LineWidth(2), RooFit.ProjectionRange('rightSideband')) 
                , logy = False
                )

