import sys
from math import sqrt, pi, cos, sin
from RooFitWrappers import *
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')
from ROOT import RooMsgService
from ROOT import kDashed, kRed, kGreen, kBlack, kOrange, TCanvas, TLatex, TFile
from P2VVGeneralUtils import plot
#RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)

import RootStyle
from ROOT import (gROOT,gStyle,TStyle)
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

from P2VVGeneralUtils import numCPU
fitOpts = dict(  Timer=1
                #, NumCPU = numCPU()
                #, NumCPU = 1
                , Save = True
                , Optimize = 1
                #, Verbose = True
                #, Minimizer = ('Minuit2','minimize')
                )

tmincut = 0.3
timeacceptance = True
fit = False

# define observables
m    = RealVar('B_s0_MM',  Title = 'm(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  50
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('J_psi_1S_MM', Title = 'm(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('phi_1020_MM', Title = 'm(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  32 )
t    = RealVar('B_s0_LOKI_DTF_CTAU',  Title = 'Decay time',    Unit = 'ps',  Observable = True, MinMax = (tmincut, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integration when making plots. Checked the data: no events below 0.007.
st   = RealVar('B_s0_LOKI_DTF_CTAU',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.007, 0.12),  nBins =  20 )

## eta_os  = RealVar('tagomega_os',      Title = '#eta_{OS}',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
## #The peak at 0.5 seems to be shifted to -2 in the SS eta!
## eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
## iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
## #The peak at 0 seems to be shifted to -1000 in the SS tagdecision
## iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
## sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
## triggerdec = Category( 'triggerDecision',            Title = 'triggerdec',                 Observable = True, States = { 'triggered': +1 } )
## from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
## #angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
## angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
##                     , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
##                     , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
##                     )

from P2VVGeneralUtils import readData

#Read data
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
                 , dataSetName = 'DecayTree'
                 , NTuple = True
                 , observables = [ m, mpsi, mphi, t, st, eta_os, eta_ss, iTag_os, iTag_ss, sel, triggerdec, angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]
                 , Rename = 'Data_DecayTree'
                 )

#TODO: Fix bug: When Rename is on, data = nil, but the dataset is imported in the ws anyway, so this is a quick fix:
#data = obj.ws().data('Data_DecayTree')

print 'Number of events', data.numEntries()

data.table(iTag_os).Print('v')
data.table(iTag_ss).Print('v')

########
# PLOT #
########
from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack, RooArgSet
from P2VVGeneralUtils import plot, Sorter

#PhiMass:
phimassplot = True
#phimassplot = False

if phimassplot:
    canvas_phimass = TCanvas('phimassplot','phimassplot',972,600)
    canvas_phimass.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts = dict()
    plot( canvas_phimass
          , mphi
          , data
          #, pdf
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.1
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          #, pdfOpts = dict( LineWidth = 3
          #                  , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
          #                  ,**pdfOpts
          #                  )
          #, plotResidHist = True
          #, logy = False
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.6,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=70
                              ) 
          )
    canvas_phimass.SaveAs('canvas_phimass.eps')

#JpsiMass:
jpsimassplot = True
#jpsimassplot = False

if jpsimassplot:
    canvas_jpsimass = TCanvas('jpsimassplot','jpsimassplot',972,600)
    canvas_jpsimass.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts = dict()
    plot( canvas_jpsimass
          , mpsi
          , data
          #, pdf
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.1
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          #, pdfOpts = dict( LineWidth = 3
          #                  , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
          #                  ,**pdfOpts
          #                  )
          #, plotResidHist = True
          #, logy = False
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.2,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=70
                              ) 
          )
    canvas_jpsimass.SaveAs('canvas_jpsimass.eps')
