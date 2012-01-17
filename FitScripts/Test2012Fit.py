from math import sqrt, pi, cos, sin
from RooFitWrappers import *

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

from P2VVGeneralUtils import numCPU
from ROOTDecorators import  ROOTversion as Rv
fitOpts = dict( NumCPU = numCPU() 
              , Timer=1
              , Save = True
              , Verbose = False
              , Optimize = True
#              , Minimizer = ('Minuit2','minimize')
              )

# define observables
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.3, 14),    nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (-2.0001,0.50001),  nBins =  25)
iTag_os = Category( 'tagdecision_os', Title = 'initial state flavour tag OS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#The peak at 0 seems to be shifted to -1000 in the SS tagdecision
iTag_ss = Category( 'tagdecision_ss', Title = 'initial state flavour tag SS', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : -1000 } )
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
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
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

#####################
#Use this formalism later for adding the SS tagger
####################

## # Make SuperCategory from (iTag_os and iTag_ss)
## TypeCat = SuperCategory('TypeCat',[biased,unbiased])
## data.table(TypeCat).Print('v')
## fitcat = MappedCategory('fitcat',TypeCat,{"AllUnbiased":["{NotBiased;Unbiased}","{Biased;Unbiased}"],"FullyBiased":["{Biased;NotUnbiased}"]})
## data.table(fitcat).Print('v')

## fitcat = data.addColumn(fitcat._var)#Whaa, hackie?
## fitcat.SetName('fitcat')
## fitcat.setRange("unbiased","AllUnbiased")
## fitcat.setRange("fullybiased","FullyBiased")

## unbiaseddata = data.reduce(RooFit.CutRange('unbiased'))
## fullybiaseddata = data.reduce(RooFit.CutRange('fullybiased'))

## unbiaseddata.SetTitle('unbiaseddata')
## unbiaseddata.SetName('unbiaseddata')
## fullybiaseddata.SetTitle('fullybiaseddata')
## fullybiaseddata.SetName('fullybiaseddata')

## fulldata = data.reduce(RooFit.CutRange('unbiased'))
## fulldata.SetTitle('fulldata')
## fulldata.SetName('fulldata')
## fulldata.append(fullybiaseddata)

# B mass pdf
from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as Signal_BMass, LP2011_Background_Mass as Background_BMass
sig_m = Signal_BMass(     Name = 'sig_m', mass = m, m_sig_mean = dict( Value = 5365, MinMax = (5363,5372) ) )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m, m_bkg_exp  = dict( Name = 'm_bkg_exp' ) )

#Time Resolution Model
#Data
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as DataTimeResolution
tresdata = DataTimeResolution( time = t, timeResSFConstraint = True ) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
externalConstraints = list()
externalConstraints += tresdata.externalConstraints()

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , deltaGamma = dict( Name = 'dGamma'
                                                            , Value = 0.060
                                                            , Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
                                                            )
                                       , deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = False) 
                                       , deltaMConstraint = True
                                      )
externalConstraints += lifetimeParams.externalConstraints()

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os, p0Constraint = True, p1Constraint = True ) # Constant = False, Constrain = True ) TODO!!!
externalConstraints += tagging.externalConstraints()

# WARNING: we don't try to describe wtag, so when plotting you must use ProjWData for eta_os !!!
#Need this, because eta_os is conditional observable in signal PDF, the actual shape doesn't matter for fitting and plotting purposes
#eta_os_pdf = UniformPdf( Name = 'eta__os_pdf', Arguments = (eta_os,) )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam(  phiCP      = dict( Name = 'phi_s'
                                              , Value = -0.04
                                              , MinMax = (-pi,pi)
                                              #Can't have kwarg Constant when blinded???
                                              #, Constant = False
                                              , Blind =  ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 ))
                         , lambdaCPSq = dict( Value = 1., Constant = True )
                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
#from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
#amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.60, A0Phase = 0
#                                              , AperpMag2 = 0.16, AperpPhase = -0.17 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
#                                              , AparPhase = 2.5
#                                              , f_S = dict( Value = 0.0, Constant = False )
#                                              , ASPhase = dict( Value = 0.0, Constant = False )
#                                             )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.60, A0Phase = 0
                                              , AperpMag2 = 0.16, AperpPhase = -0.17 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                              , AparPhase = 2.5
                                              , f_S_Re = dict( Value = 0.10 / ( 1. + 0.10 ) * cos(2.20), Constant = False )
                                              , f_S_Im = dict( Value = 0.10 / ( 1. + 0.10 ) * sin(2.20), Constant = False )
                                             )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
#from P2VVParameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet
#amplitudes = JpsiVPolar_AmplitudeSet(  A0Mag2 = 0.60, A0Phase = 0
#                                     , AperpMag2 = 0.16, AperpPhase = -0.17 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
#                                     , AparPhase = 2.5
#                                     , ASMag2 = dict( Value = 0.01, Constant = False )
#                                     , ASPhase = dict( Value = 0.5, Constant = False )
#                                     , PWaveNorm = False
#                                    )

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions
                                                  , amplitudes
                                                  , CP
                                                  , Product('tag',(iTag_os,tagging['dilution']))
                                                  , ['A0','Apar','Aperp','AS'] ) 

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tresdata.model()
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.3 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.92, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.21, MinMax = (0.01,0.5) ) )

sig_t_angles = BDecay( Name      = 'sig_t_angles'
                       , time      = t
                       , dm        = lifetimeParams['deltaM'] 
                       , tau       = lifetimeParams['MeanLifetime']
                       , dGamma    = lifetimeParams['deltaGamma'] 
                       , resolutionModel = tresdata.model()
                       , coshCoef  = basisCoefficients['cosh']
                       , cosCoef   = basisCoefficients['cos']
                       , sinhCoef  = basisCoefficients['sinh']
                       , sinCoef   = basisCoefficients['sin']
#                       , ConditionalObservables = ( eta_os, )
#                       , ConditionalObservables = ( eta_os, iTag_os, )
                       )

#####################################
### Angular acceptance correction ###
#####################################
from P2VVGeneralUtils import RealMomentsBuilder
#nset = angles.angles.values()

#canomoms = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, MCpdf,nset) for v in angles.functions.itervalues() for i in v if i ) )
#canomoms.compute(MCdata)
#canomoms.Print(Scales = [1./(16.*sqrt(pi)),1./(16.*sqrt(pi)),1./(16.*sqrt(pi))])

from itertools import chain
#momindices = indices(3,3)
momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/effmoments.txt')
eff.Print()

#Build Angular acceptance corrected PDF
#sig_t_angles = eff * sig_t_angles

##############################
### Proper time acceptance ###
##############################
a = RealVar('a', Title = 'a', Value = 1.45, MinMax = (1, 2))
c = RealVar('c', Title = 'c', Value = -2.37, MinMax = (-3, 2))
eff = FormulaVar('eff_shape', "(@0 > 0.) ? (1 / (1 + (@1 * @0) ** (@2))) : 0.0001", [t, a, c])

from P2VVBinningBuilders import build1DVerticalBinning
binning, eff_func = build1DVerticalBinning('time_binning', eff, t, 0.05, 1.)

acceptance = BinnedPdf(Name = 'time_acceptance', Observables = [t], Function = eff, Binning = binning)

#Build proper time acceptance corrected PDF
#sig_t_angles = acceptance * sig_t_angles

##################
### Build PDFs ###
##################

############
# BKG ONLY #
############
sidebanddata = data.reduce(CutRange = 'leftsideband')
rightsidebanddata = data.reduce(CutRange = 'rightsideband')
sidebanddata.append(rightsidebanddata)

nbkg = 20000
#background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf()), Yield = ( nbkg, 0, 2.0*nbkg) )
#The following doesn't make a difference, indeed!
background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf()), Yield = ( nbkg, 0, 2.0*nbkg) )
background[eta_os]=None
background[iTag_os]=None

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
    for i in angles.angles.itervalues():
        background[i]=None
    #background += UniformPdf( Name = 'bkg_angles', Arguments = angles.angles.itervalues() )
else :
    background += HistPdf( Name = 'bkg_angles'
                             , Observables = angles.angles.itervalues()
                             , Binning =  { angles.angles['cpsi'] : 7
                                          , angles.angles['ctheta'] : 5
                                          , angles.angles['phi' ] : 9
                                          }
                             , Data  = sidebanddata
                             )

bkgpdf = buildPdf((background,), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name = 'bkgpdf')

#bkgresult = bkgpdf.fitTo(sidebanddata,**fitOpts)
#bkgresult .writepars('bkgresult',False)

############
# SIG ONLY #
############
sigdata = data.reduce(CutRange = 'signal')
sigdata.Print()

nsig = 20000
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles ), Yield = ( nsig, 0, 2.0*nsig) )
#sig_t_angles depends on eta_os, and is NOT conditional on eta_os, so also ask for eta_os here....
sigpdf = buildPdf((signal,), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name = 'sigpdf')

#for i in ['ASPhase','f_S'] :
#for i in ['ASPhase','ASMag2'] :
#    amplitudes[i].setVal(0.0)
#    amplitudes[i].setConstant(True)

#amplitudes['ASPhase'].setVal(2.5)
#amplitudes['ASPhase'].setConstant(False)
#amplitudes['f_S'].setVal(0.01)
#amplitudes['f_S'].setConstant(False)

#amplitudes['ASPhase'].setVal(2.5)
#amplitudes['ASPhase'].setConstant(False)
#amplitudes['ASMag2'].setVal(0.01)
#amplitudes['ASMag2'].setConstant(False)

sigpdf.Print()
sigresult = sigpdf.fitTo(sigdata, ExternalConstraints = externalConstraints, **fitOpts)
sigresult.writepars('signalfitresult_As2Free_AngAcc_NoTimeAcc', False)

assert False
#############
# MASS ONLY #
#############
masspdf = buildPdf((signal,background), Observables = (m,), Name = 'masspdf')
#masspdf.fitTo(data,**fitOpts)

###########
# SIG+BKG #
###########
#sig_t_angles depends on eta_os, and is NOT conditional on eta_os, so also ask for eta_os here....
pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name='fullpdf')
pdf.Print()

###################
### CLASSIC FIT ###
###################

amplitudes['ASPhase'].setVal(2.5)
amplitudes['ASPhase'].setConstant(False)
amplitudes['f_S'].setVal(0.01)
amplitudes['f_S'].setConstant(False)

#Don't add externalconstraints to fitOpts, otherwise fits for splots might go wrong, you don't want to constrain mass fits!
#classicfitresult = pdf.fitTo(data,ExternalConstraints = externalConstraints, **fitOpts)
#classicfitresult.writepars('classicfitresult',False)

############
### SFIT ###
############
# make sweighted dataset. TODO: using J/psi phi mass
from P2VVGeneralUtils import SData, splot
from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
splot_m = SData(Pdf = masspdf, Data = data, Name = 'MassSplot')

S_sigdata = splot_m.data('signal')
S_bkgdata = splot_m.data('bkg')

## canvas2 = TCanvas()
## canvas2.Divide(2)
## canvas2.cd(1)
## mframe = m.frame()
## S_sigdata.plotOn(mframe)
## mframe.Draw()
## canvas2.cd(2)
## mframe = m.frame()
## S_bkgdata.plotOn(mframe)
## mframe.Draw()

sfitsignal = Component('sfitsignal', ( sig_t_angles, ), Yield = ( nsig, 0, 2.0*nsig) )
sfitpdf = buildPdf((sfitsignal,), Observables = (t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name='sfitpdf')

#Don't add externalconstraints to fitOpts, otherwise fits for splots might go wrong, you don't want to constrain mass fits!
#sfitresult = sfitpdf.fitTo(S_sigdata,ExternalConstraints = externalConstraints, **fitOpts)
#sfitresult.writepars('sfitresult',False)

################
# FROM GERHARD #
################

# fit & fix iTag_os parameters
#pdf_itag   = buildPdf((signal,background), Observables = (m,iTag_os), Name='pdf_itag')
#pdf_itag.fitTo( data,**fitOpts)
#for p in pdf_itag.Parameters() : p.setConstant( not p.getAttribute('Yield') )

from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
from P2VVGeneralUtils import plot

#Old bkganglecanvas
#for (p,o) in zip(canvas.pads(3),angles.angles.itervalues()):
    #p.cd()
    #f = o.frame()
    #sidebanddata.plotOn(f)
    ##background[angles.angles.itervalues()].dataHist().plotOn(f)
    #bkganglepdf.plotOn(f)
    #f.Draw()

bkganglecanvas = TCanvas()
for (p,o) in zip(bkganglecanvas.pads(3),angles.angles.itervalues()):
    plot(p,o,sidebanddata,bkganglepdf)

canvas = dict()
for rng in ( None, 'signal','leftsideband,rightsideband' ) :
    canvas[rng] = TCanvas('%s'%rng)
    dataOpts = dict( CutRange =        rng ) if rng else dict()
    pdfOpts  = dict( ProjectionRange = rng ) if rng else dict()
    from ROOT import RooArgSet
    # TODO: grab the data in the relevant range... data.reduce( **dataOpts ) 
    #       and remove the spike at 0.5 to take into account its correlation to iTag_os = 0!!!
    pdfOpts['ProjWData'] = ( RooArgSet( eta_os._var ),  data, True ) 
    obs =  [ o for o in pdf.Observables() if hasattr(o,'frame') ]
    for (p,o) in zip( canvas[rng].pads(len(obs)), obs ) :
        plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                            , 'bkg*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                            }
                             , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
                             , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
                             , logy = ( o == t )
                             )
