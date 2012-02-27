from math import sqrt, pi, cos, sin
from RooFitWrappers import *
indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')
from ROOT import RooMsgService
#RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
#RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)

from P2VVGeneralUtils import numCPU
fitOpts = dict(
    NumCPU = 1
    #NumCPU = numCPU()
    , Timer=1
    , Save = True
    , Verbose = True
    , Minimizer = ('Minuit2','minimize')
    , Optimize = 2
    )

tmincut = 0.3

# define observables
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (tmincut, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integration when making plots. Checked the data: no events below 0.007.


st   = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps',  Observable = True, MinMax = (0.007, 0.12),  nBins =  50)
# st   = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = False, Value = 0.1, MinMax = (0.007, 0.12), nBins = 50, Constant = True)

eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#The peak at 0 seems to be shifted to -1000 in the SS tagdecision
iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
triggerdec = Category( 'triggerDecision',            Title = 'triggerdec',                 Observable = True, States = { 'triggered': +1 } )
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

from P2VVGeneralUtils import readData

#Read data
data = readData( '/stuff/PhD/p2vv/data/Bs2JpsiPhi_ntupleB_for_fitting_20120120.root'
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

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

#Signal: per-event error
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
sigtres = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st, Cache = False)

#Background
bkgtres = sigtres

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = bkgtres.model()
                         , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.21 )
                         , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.06, MinMax = (0.5,2.5) )
                         , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.15, MinMax = (0.01,0.5) )
                         )


##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
MoriondAcceptance = Moriond2012_TimeAcceptance( time = t, Input = '/stuff/PhD/p2vv/data/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')

from ROOT import TH1F
flat_hist = TH1F('flat', 'flat', 20, 0.3, 14)
for i in range(1, flat_hist.GetNbinsX() + 1):
    flat_hist.SetBinContent(i, 1.)
UnityAcceptance = HistFunc(Name = 'unity_time_acceptance', Observables = [t], Histogram = flat_hist)


##########################################
### Plot proper time acceptance curves ###
##########################################

from ROOT import TCanvas
unitycanvas = TCanvas()
tframe = t.frame()
UnityAcceptance.plotOn(tframe)
tframe.Draw()

moriondcanvas = TCanvas()
tframe = t.frame()
MoriondAcceptance._acceptance.plotOn(tframe)
tframe.Draw()

############
# BKG COMP #
############
sidebanddata =      data.reduce(CutRange = 'leftsideband' )
sidebanddata.append(data.reduce(CutRange = 'rightsideband'))

nbkg = 10500

######################################################
### Turn these three on/off to see the problem !!! ###
######################################################
#background = Component('bkg'   , ( MoriondAcceptance*bkg_t.pdf(),), Yield = ( nbkg, 0.9, 1.1*nbkg) )
#background = Component('bkg'   , ( UnityAcceptance*bkg_t.pdf(),), Yield = ( nbkg, 0.9, 1.1*nbkg) )
background = Component('bkg'   , ( bkg_t.pdf(),), Yield = ( nbkg, 0.9, 1.1*nbkg) )

###############
### BKG fit ###
###############
bkgpdf = buildPdf((background,), Observables = (t,), Name='bkgpdf')
bkgpdf.Print()

## from ROOT import RooArgSet
## nset = RooArgSet(t._target_())
## print bkgpdf.createIntegral(nset).getVal()

## RooMsgService.instance().addStream(RooFit.INFO, RooFit.Topic(RooFit.Integration))
## RooAbsPdf.verboseEval(2)

## print bkgpdf_flat.createIntegral(nset).getVal()

## boundaries = [0.3, 0.985, 1.67, 2.355, 3.04, 3.725, 4.41, 5.095, 5.78, 6.465, 7.15, 7.835, 8.52, 9.205, 9.89, 10.575, 11.26, 11.945, 12.63, 13.315, 14]
## ints = []
## for i in range(len(boundaries) - 1):
##     r = 'test_range_%d' % i
##     t.setRange(r, (boundaries[i], boundaries[i + 1]))
##     ints.append(bkgpdf.createIntegral(nset, r))

## s = 0.
## for i in ints:
##     s += i.getVal()
## print s


bkgfitresult = bkgpdf.fitTo(sidebanddata, **fitOpts)

projWDataSet = []
projWDataSet += [ st ]
projWData     = dict( ProjWData = ( sidebanddata.reduce(  ArgSet = projWDataSet ), True ) )

bkgtimeplot = True
if bkgtimeplot:
    from ROOT import kDashed, kRed, kGreen, TCanvas, TLatex
    from P2VVGeneralUtils import plot
    canvas = TCanvas()
    plot( canvas
          , t
          , sidebanddata
          , bkgpdf
          #, components = {'bkg_m' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed   )
                          #'sig_m' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen )
          #                }
          , pdfOpts = dict( list( projWData.items() ), LineWidth = 3 )
          #, pdfOpts = dict( LineWidth = 3 )
          , plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = 'B_{s}#rightarrow J/#psi#phi'
                              , TitleOffset = (1.2,'y')
                              , Object = ( TLatex(0.55,.8,"#splitline{LHCb preliminary}{#sqrt{s} = 7 TeV, L = 1.03 fb^{-1}}", NDC = True), )
                              , Bins=70 ) 
          )

