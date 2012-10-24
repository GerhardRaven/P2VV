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
m    = RealVar('mass',  Title = 'm(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  50
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'm(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'm(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'Decay time',    Unit = 'ps',  Observable = True, MinMax = (tmincut, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integration when making plots. Checked the data: no events below 0.007.
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.007, 0.12),  nBins =  20 )
eta_os  = RealVar('tagomega_os',      Title = '#eta_{OS}',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  20)
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

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372) ) )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

#Time Resolution Model
#Three Gaussians
#from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as DataTimeResolution
#sigtres = DataTimeResolution( time = t, timeResSFConstraint = True ) 

#Signal: per-event error
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
sigtres = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st)

#Background
#from P2VVParameterizations.TimeResolution import Truth_TimeResolution as BkgTimeResolution
#from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as BkgTimeResolution
#from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as BkgTimeResolution
#bkgtres = BkgTimeResolution( time = t)
bkgtres = sigtres

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , dGamma = dict( Name = 'dGamma'
                                                            , Value = 0.060
                                                            #, Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
                                                            )
                                       , dM = dict( Value = 17.58, MinMax = (16.5,18.5), Constant = False)
                                       #, dMConstraint = True
                                      )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os, p0Constraint = True, p1Constraint = True )
tagging.addConditional(iTag_os)

from P2VVParameterizations.CPVParams import LambdaArg_CPParam
CP = LambdaArg_CPParam( phiCP      = dict( Name = 'phi_s'
                                              , Value = -0.04
                                              , MinMax = (-pi,pi)
                                              #, Blind =  ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
                                              )
                           )

#from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
#CP = LambdaSqArg_CPParam(  phiCP      = dict( Name = 'phi_s'
#                                              , Value = -0.04
#                                              , MinMax = (-pi,pi)
#                                              #Can't have kwarg Constant when blinded???
#                                              #, Constant = False
#                                              #, Blind =  ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 ))
#                         , lambdaCPSq = dict( Value = 1., Constant = True )
#                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.52, A0Phase = 0
                                              , AperpMag2 = 0.25, AperpPhase = 2.77 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                              , AparPhase = 3.2
                                              , f_S = dict( Value = 0.02, Constant = False )
                                              , ASPhase = dict( Value = 2.7, Constant = False )
                                             )

# need to specify order in which to traverse...
from RooFitWrappers import RealCategory
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes, CP, ['A0','Apar','Aperp','AS'] )

from RooFitWrappers import BTagDecay
sig_t_angles_iTag = BTagDecay(  Name                   = 'sig_t_angles_iTag'
                              , time                   = t
                              , iTag                   = iTag_os
                              , dm                     = lifetimeParams['dM']
                              , tau                    = lifetimeParams['MeanLifetime']
                              , dGamma                 = lifetimeParams['dGamma']
                              , resolutionModel        = sigtres['model']
                              , coshCoef               = basisCoefficients['cosh']
                              , cosCoef                = basisCoefficients['cos']
                              , sinhCoef               = basisCoefficients['sinh']
                              , sinCoef                = basisCoefficients['sin']
                              , dilution               = tagging['dilution']
                              , ADilWTag               = tagging['ADilWTag']
                              , avgCEven               = tagging['avgCEven']
                              , avgCOdd                = tagging['avgCOdd']
                              , ConditionalObservables = sigtres.conditionalObservables() + tagging.conditionalObservables()
                              , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                                         + sigtres.externalConstraints()\
                                                         + tagging.externalConstraints()
                             )

from P2VVParameterizations.FlavourTagging import Trivial_TagPdf
sig_itagospdf = Trivial_TagPdf( iTag_os, NamePrefix = 'os_tag_sig' ).pdf()
bkg_itagospdf = Trivial_TagPdf( iTag_os, NamePrefix = 'os_tag_bkg' ).pdf()
# TODO: count # of iTag_os = +1,0,-1 in sweighted data set in order to set the additonal
#       parameters to their eventual fit values...

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = bkgtres.model()
                         , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.21 )
                         , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.06, MinMax = (0.5,2.5) )
                         , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.15, MinMax = (0.01,0.5) )
                         )

print "*****************************"
print [ i.GetName() for i in sig_t_angles_iTag.ConditionalObservables()  ]
print "*****************************"

#####################################
### Angular acceptance correction ###
#####################################
from P2VVGeneralUtils import RealMomentsBuilder
#nset = angles.angles.values()

#canomoms = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, MCpdf,nset) for v in angles.functions.itervalues() for i in v if i ) )
#canomoms.compute(MCdata)
#canomoms.Print(Scales = [1./(16.*sqrt(pi)),1./(16.*sqrt(pi)),1./(16.*sqrt(pi))])

from itertools import chain
#momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))
#These are the relevant terms as found with MinSignificance>3
momindices = [(0,0,0),(0,2,0),(0,2,2),(2,0,0)]

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/data/bfys/dveijk/DataJpsiPhi/2012/effmoments_tcut_%s.txt'%(str(tmincut)))

#Build Angular acceptance corrected PDF
sig_t_angles_iTag = eff * sig_t_angles_iTag

##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
if timeacceptance:
    sig_t_angles_iTag = acceptance * sig_t_angles_iTag
sig_t_angles_iTag.setAttribute("NOCacheAndTrack")

#######################
### S weighted data ###
#######################
nsig = 21000
nbkg = 10500
from P2VVGeneralUtils import createSData
from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
sig_m_comp = Component('sig_m_comp', ( sig_m.pdf(),), Yield = ( nsig, 0.9, 1.1*nsig) )
bkg_m_comp = Component('bkg_m_comp', ( bkg_m.pdf(),), Yield = ( nbkg, 0.9, 1.1*nbkg) )
masspdf = buildPdf((sig_m_comp, bkg_m_comp), Observables = (m,), Name='masspdf')
masspdf.fitTo(data,**fitOpts)
from P2VVGeneralUtils import SData
splot_m = SData(Pdf = masspdf, Data = data, Name = 'splot_m')
sigdata = splot_m.data('sig_m_comp')
bkgdata = splot_m.data('bkg_m_comp')

massplot = False

if massplot:
    canvas_m = TCanvas('massplot','massplot',972,600)
    canvas_m.SetFixedAspectRatio(True)
    plot( canvas_m
          , m
          , data
          , masspdf
          , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , yTitleOffset = 1.2
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack)
          , pdfOpts = dict( LineWidth = 3 )
          #, plotResidHist = True
          , logy = False
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=70
                              ) 
          )

########################
### Build omega PDFs ###
########################
sig_omegaospdf = HistPdf( Name = 'sig_omegaospdf', Observables = (eta_os,), Data = sigdata )
bkg_omegaospdf = HistPdf( Name = 'bkg_omegaospdf', Observables = (eta_os,), Data = bkgdata )

#from RooFitWrappers import UniformPdf
#sig_omegaospdf =  UniformPdf('sig_omegaospdf', Arguments = ( eta_os, ) )
#bkg_omegaospdf =  UniformPdf('bkg_omegaospdf', Arguments = ( eta_os, ) )

#########################
### Build sigmat PDFs ###
#########################
sig_sigmatpdf = HistPdf( Name = 'sig_sigmatpdf', Observables = (st,), Data = sigdata )
bkg_sigmatpdf = HistPdf( Name = 'bkg_sigmatpdf', Observables = (st,), Data = bkgdata )

#from RooFitWrappers import UniformPdf
#sig_sigmatpdf =  UniformPdf('sig_sigmatpdf', Arguments = ( st, ) )
#bkg_sigmatpdf =  UniformPdf('bkg_sigmatpdf', Arguments = ( st, ) )

############
# BKG COMP #
############
background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf(), bkg_itagospdf, bkg_omegaospdf, bkg_sigmatpdf), Yield = ( nbkg, 0.9, 1.1*nbkg) )

### BKG angles 

# create PDF for angular background
if False :
    # make sweighted dataset using J/psi phi mass
    from P2VVGeneralUtils import createSData
    from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
    splot_m = createSData( Name = 'Mass' , Components =  (signal,background), Observables = (m,), FitOpts = fitOpts, Data = data )
    mompdfBuilder = SPlot_Moment_Angles( angles.angles , splot_m )
    background += mompdfBuilder.pdf( Component = background.GetName()
                                       , Indices = [ i for i in indices(3,4) ]
                                       , Name = 'bkg_angles'
                                       , MinSignificance = 0.5
                                       , Scale = sqrt(50.) )
elif False:
    for i in angles.angles.itervalues(): background[i]=None
else :
    background += HistPdf( Name = 'bkg_angles'
                             , Observables = angles.angles.itervalues()
                             , Binning =  { angles.angles['cpsi']   : 5
                                          , angles.angles['ctheta'] : 7
                                          , angles.angles['phi' ]   : 9
                                          }
                             , Data  = bkgdata
                             )

############
# SIG COMP #
############
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles_iTag, sig_itagospdf, sig_omegaospdf, sig_sigmatpdf), Yield = ( nsig, 0.9, 1.1*nsig) )

#######
# FIT #
#######

#pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os)+tuple(angles.angles.itervalues()), Name='fullpdf')
#pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name='fullpdf')
pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os,eta_os,st)+tuple(angles.angles.itervalues()), Name='fullpdf')
#pdf   = buildPdf((signal,background), Observables = (m,t)+tuple(angles.angles.itervalues()), Name='fullpdf')
#pdf   = buildPdf((signal,background), Observables = (m,t)+tuple(angles.angles.itervalues()), Name='fullpdf')

def search(fname,path) :
    import os
    for f in ( os.path.join(p,fname) for p in os.path.split(os.pathsep) ) :
        if os.path.exists(f) : return f
    return None

import os

if timeacceptance:
    fitname = 'cfit_acc_floatdeltams'
else:
    fitname = 'cfit_noacc_floatdeltams'

paramfile = search(fitname+'params.txt',os.pathsep.join(['.','FitScripts']) )
if paramfile :
    print 'Reading fit result from %s' % paramfile
    fitset = pdf.getParameters(data)
    fitset.readFromFile(paramfile)

if fit or not paramfile:
    cfitresult = pdf.fitTo(data, **fitOpts)
    f = TFile(fitname+"_result","RECREATE")
    cfitresult.Write(fitname+"_result") 
    f.Close()
    cfitresult.writepars(fitname+'result',False)
    cfitresult.Print()
    fitset = pdf.getParameters(data)
    fitset.writeToFile(fitname+"params.txt")

print 120 * '='
print 'DvEFit: parameters in NLL:'
for par in pdf.getVariables() : par.Print()
print 120 * '='

assert False

def correctgamma():
    gamma = obj.ws()['Gamma']
    uncval = (1.)/(gamma.getVal())
    alpha = 0.0112
    from math import sqrt
    corrval = ( (1.+alpha*uncval)/(4.*alpha) ) - ( (1.)/(4.*alpha) ) * sqrt(alpha**2*uncval**2-6*alpha*uncval+1.)
    print 'The fitted value of Gamma is %f +/- %f and is corrected for time acceptance as %f +/- %f'%(gamma.getVal(),gamma.getError(),1./corrval,gamma.getError())
    return

#Correct for Gamma!
#correctgamma()

def checkambiguity():
##     fitresult = search(fitname+'_result',os.pathsep.join(['.','FitScripts']) )
##     if fitresult:
##         f = TFile(fitname+"_result")
##         fitresultpresent = f.Get(fitname+"_result")
    print 'PDF value is %f'%(pdf.getVal())
    ll = pdf._var.createNLL(data)
    print 'NLL is %f'%(ll.getVal())
    print "Changing the parameters to the discrete ambiguity:"
    obj.ws()['dGamma'].setVal(-1*obj.ws()['dGamma'].getVal())
    obj.ws()['phi_s'].setVal(pi-obj.ws()['phi_s'].getVal())
    obj.ws()['ASPhase'].setVal(-1*obj.ws()['ASPhase'].getVal())
    obj.ws()['AparPhase'].setVal(-1*obj.ws()['AparPhase'].getVal())
    obj.ws()['AperpPhase'].setVal(pi-obj.ws()['AperpPhase'].getVal())
    print 'PDF value in ambiguity is %f'%(pdf.getVal())
    lltwo = pdf._var.createNLL(data)
    print 'Value of ll in ambiguity is %f'%(lltwo.getVal())

########
# PLOT #
########
from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack, RooArgSet
from P2VVGeneralUtils import plot, Sorter

####################
### Full region: ###
####################
#Mass:
#massplot = True
massplot = False

if massplot:
    canvas_m = TCanvas('massplot','massplot',972,600)
    canvas_m.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts = dict()
    plot( canvas_m
          , m
          , data
          , pdf
          , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , yTitleOffset = 1.2
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = False
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=70
                              ) 
          )
    canvas_m.SaveAs('canvas_m_%s.eps'%(fitname))

#Angles:
#angleplot = True
angleplot = False

if angleplot:
    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]
    canvas_angles = TCanvas('anglesplot','anglesplot',972,600)
    canvas_angles.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts = dict()
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas_angles.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from ROOT import RooArgSet
        plot( p, o, data, pdf
              , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                              #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                              ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                              }
              , yTitleOffset = 1.5
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,)
              , pdfOpts  = dict( LineWidth = 2
                                 , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                                 , **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.2,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  , Bins = 25
                                  )
              , logy = ( o == t )
              )
    canvas_angles.SaveAs('canvas_angles_%s.eps'%(fitname))

#Time:
#timeplot = True
timeplot = False

if timeplot:
    canvas_t = TCanvas('timeplot','timeplot',972,600)
    canvas_t.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts = dict()
    plot( canvas_t
          , t
          , data
          , pdf
          , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , yTitleOffset = 1.2
          , xTitleOffset = 1.1
          , yScale = (0.1,10000)
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137
                              ) 
          )
    canvas_t.SaveAs('canvas_t_%s.eps'%(fitname))

#####################
### Signal Region ###
#####################
#Angles:
#anglesplot_sig = True
anglesplot_sig = False

if anglesplot_sig:
    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]

    canvas_angles_sig = TCanvas('anglesplot_sig','anglesplot_sig',972,600)
    canvas_angles_sig.SetFixedAspectRatio(True)
    dataOpts = dict( CutRange =        'signal' )
    pdfOpts  = dict( ProjectionRange = 'signal' ) 
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas_angles_sig.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from ROOT import RooArgSet
        plot( p, o, data, pdf
              , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                              #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                              ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                              }
              , yTitleOffset = 1.5
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,**dataOpts)
              , pdfOpts  = dict( LineWidth = 2
                                 , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                                 , **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.25,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  ,Bins = 25
                                  )
              , logy = ( o == t )
              )
    canvas_angles_sig.SaveAs('canvas_angles_sig_%s.eps'%(fitname))

#Time:
#timeplot_sig = True
timeplot_sig = False

if timeplot_sig:
    canvas_t_sig = TCanvas('timeplot_sig','timeplot_sig',972,600)
    canvas_t_sig.SetFixedAspectRatio(True)
    dataOpts = dict( CutRange =        'signal' )
    pdfOpts  = dict( ProjectionRange = 'signal' ) 
    plot( canvas_t_sig
          , t
          , data
          , pdf
          , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , yTitleOffset = 1.2
          , xTitleOffset = 1.1
          , yScale = (0.1,10000)
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137
                              ) 
          )
    canvas_t_sig.SaveAs('canvas_t_sig_%s.eps'%(fitname))

#########################
### Background Region ###
#########################
#Angles:
#anglesplot_bkg = True
anglesplot_bkg = False

if anglesplot_bkg:
    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]

    canvas_angles_bkg = TCanvas('anglesplot_bkg','anglesplot_bkg',972,600)
    canvas_angles_bkg.SetFixedAspectRatio(True)
    dataOpts = dict( CutRange =        'leftsideband,rightsideband' )
    pdfOpts  = dict( ProjectionRange = 'leftsideband,rightsideband' ) 
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas_angles_bkg.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from ROOT import RooArgSet
        plot( p, o, data, pdf
              , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                              #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                              ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                              }
              , yTitleOffset = 1.5
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,**dataOpts)
              , pdfOpts  = dict( LineWidth = 2
                                 , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                                 , **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.25,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  , Bins = 25
                                  )
              , logy = ( o == t )
              )
    canvas_angles_bkg.SaveAs('canvas_angles_bkg_%s.eps'%(fitname))

#Time:
#timeplot_bkg = True
timeplot_bkg = False

if timeplot_bkg:
    canvas_t_bkg = TCanvas('timeplot_bkg','timeplot_bkg',972,600)
    canvas_t_bkg.SetFixedAspectRatio(True)
    dataOpts = dict( CutRange =        'leftsideband,rightsideband' )
    pdfOpts  = dict( ProjectionRange = 'leftsideband,rightsideband' ) 
    plot( canvas_t_bkg
          , t
          , data
          , pdf
          , components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                          #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
                          ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
                          }
          , yTitleOffset = 1.2
          , yScale = (0.1,10000)
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( data.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137
                              ) 
          )
    canvas_t_bkg.SaveAs('canvas_t_bkg_%s.eps'%(fitname))

###########################
### S Plot based signal ###
###########################
if timeacceptance:
    signalcomponent = obj._ws['signal_BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins_X_sig_t_angles_iTag_x_Eff_X_sig_m_X_sig_omegaospdf_X_sig_sigmatpdf_X_os_tag_sig_Trivial_TagPdf']
else:
    signalcomponent = obj._ws['signal_sig_t_angles_iTag_x_Eff_X_sig_m_X_sig_omegaospdf_X_sig_sigmatpdf_X_os_tag_sig_Trivial_TagPdf']

#Time:
#timeplot_sig_s = True
timeplot_sig_s = False

if timeplot_sig_s:
    canvas_t_sig_s = TCanvas('timeplot_sig_s','timeplot_sig_s',972,600)
    canvas_t_sig_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    plot( canvas_t_sig_s
          , t
          , sigdata
          , signalcomponent
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
          #                #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.2
          , yScale = (0.1,10000)
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( sigdata.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137
                              ) 
          )
    canvas_t_sig_s.SaveAs('canvas_t_sig_s_%s.eps'%(fitname))

#Angles:
#anglesplot_sig_s = True
anglesplot_sig_s = False

if anglesplot_sig_s:
    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]

    canvas_angles_sig_s = TCanvas('anglesplot_sig_s','anglesplot_sig_s',972,600)
    canvas_angles_sig_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas_angles_sig_s.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from ROOT import RooArgSet
        plot( p, o, sigdata, signalcomponent
              #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                              #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
              #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
              #                }
              , yTitleOffset = 1.5
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8
                                 , MarkerColor = kBlack
                                 ,**dataOpts)
              , pdfOpts  = dict( LineWidth = 2
                                 , ProjWData = ( sigdata.reduce(RooArgSet(st,eta_os)), True )
                                 , **pdfOpts
                                 )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.25,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  , Bins = 25
                                  )
              , logy = ( o == t )
              )
    canvas_angles_sig_s.SaveAs('canvas_angles_sig_s_%s.eps'%(fitname))

#Eta:
#etaplot_sig_s = True
etaplot_sig_s = False

if etaplot_sig_s:
    canvas_eta_sig_s = TCanvas('etaplot_sig_s','etaplot_sig_s',972,600)
    canvas_eta_sig_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    plot( canvas_eta_sig_s
          , eta_os
          , sigdata
          , signalcomponent
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
          #                #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.2
          #, yScale = (0.1,100000)
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( bkgdata.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          #, logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=137
                              ) 
          )
    canvas_eta_sig_s.SaveAs('canvas_eta_sig_s_%s.eps'%(fitname))

###############################
### S Plot based background ###
###############################
bkgcomponent = obj._ws['bkg_bkg_angles_X_bkg_m_X_bkg_omegaospdf_X_bkg_sigmatpdf_X_bkg_t_X_os_tag_bkg_Trivial_TagPdf']

#Time:
#timeplot_bkg_s = True
timeplot_bkg_s = False

if timeplot_bkg_s:
    canvas_t_bkg_s = TCanvas('timeplot_bkg_s','timeplot_bkg_s',972,600)
    canvas_t_bkg_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    plot( canvas_t_bkg_s
          , t
          , bkgdata
          , bkgcomponent
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
          #                #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.2
          , yScale = (0.1,10000)
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( bkgdata.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          , logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              , Bins=137
                              ) 
          )
    canvas_t_bkg_s.SaveAs('canvas_t_bkg_s_%s.eps'%(fitname))

#Angles:
#anglesplot_bkg_s = True
anglesplot_bkg_s = False

if anglesplot_bkg_s:
    orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )
    obs = [angles.angles['cpsi'] ,angles.angles['ctheta'] ,angles.angles['phi']]
    canvas_angles_bkg_s = TCanvas('anglesplot_bkg_s','anglesplot_bkg_s',972,600)
    canvas_angles_bkg_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    obs =  filter( lambda x : hasattr(x,'frame'), obs ) 
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas_angles_bkg_s.pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        from ROOT import RooArgSet
        plot( p, o, bkgdata, bkgcomponent
              #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
                              #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
              #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
              #                }
              , yTitleOffset = 1.5
              , xTitleOffset = 1.1
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack,**dataOpts)
              , pdfOpts  = dict( LineWidth = 2
                                 , ProjWData = ( bkgdata.reduce(RooArgSet(st,eta_os)), True )
                                 , **pdfOpts )
              , plotResidHist = False
              , frameOpts = dict( Object = ( TLatex(0.30,.25,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                                  , Title = ""
                                  , Bins = 25
                                  )
              , logy = ( o == t )
              )
    canvas_angles_bkg_s.SaveAs('canvas_angles_bkg_s_%s.eps'%(fitname))

#Eta:
etaplot_bkg_s = True
#etaplot_bkg_s = False

if etaplot_bkg_s:
    canvas_eta_bkg_s = TCanvas('etaplot_bkg_s','etaplot_bkg_s',972,600)
    canvas_eta_bkg_s.SetFixedAspectRatio(True)
    dataOpts = dict()
    pdfOpts  = dict()
    plot( canvas_eta_bkg_s
          , eta_os
          , bkgdata
          , bkgcomponent
          #, components = {'sig*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen   )
          #                #,'psi*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kOrange )
          #                ,'bkg*' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed )
          #                }
          , yTitleOffset = 1.2
          , yScale = (-0.1,7000)
          , xTitleOffset = 1.1
          , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts)
          , pdfOpts = dict( LineWidth = 3
                            , ProjWData = ( bkgdata.reduce(RooArgSet(st,eta_os)), True )
                            ,**pdfOpts
                            )
          #, plotResidHist = True
          #, logy = True
          , frameOpts = dict( Title = ''
                              , Object = ( TLatex(0.55,.8,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}", NDC = True), )
                              #, Bins=137
                              ) 
          )
    canvas_eta_bkg_s.SaveAs('canvas_eta_bkg_s_%s.eps'%(fitname))
