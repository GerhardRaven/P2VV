from math import sqrt, pi, cos, sin
from RooFitWrappers import *

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

from P2VVGeneralUtils import numCPU
fitOpts = dict( NumCPU = numCPU() 
              , Timer=1
              , Save = True
              , Verbose = False
              , Minimizer = ('Minuit2','minimize')
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
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.12),  nBins =  50 )
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
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20120120.root'
#                 , '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
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
#Three Gaussians
#from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as DataTimeResolution
#tresdata = DataTimeResolution( time = t, timeResSFConstraint = True ) 
#Per event error
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
tresdata = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st)

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , deltaGamma = dict( Name = 'dGamma'
                                                            , Value = 0.060
                                                            , Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
                                                            )
                                       , deltaM = dict( Value = 17.58, MinMax = (16.5,18.5), Constant = False) 
                                       , deltaMConstraint = True
                                      )

# define tagging parameter 
# WARNING: we don't try to describe wtag, so when plotting you must use ProjWData for eta_os !!!
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams( estWTag = eta_os, p0Constraint = True, p1Constraint = True )

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
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.52, A0Phase = 0
                                              , AperpMag2 = 0.25, AperpPhase = 2.77 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                              , AparPhase = 3.2
                                              , f_S = dict( Value = 0.02, Constant = False )
                                              , ASPhase = dict( Value = 2.7, Constant = False )
                                             )

# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions
                                                  , amplitudes
                                                  , CP
                                                  , Product('tag',(iTag_os,tagging['dilution']))
                                                  , ['A0','Apar','Aperp','AS'] ) 
basisCoefficients.externalConstraints = tagging.externalConstraints()

from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tresdata.model()
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.21 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 1.06, MinMax = (0.5,2.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.15, MinMax = (0.01,0.5) ) )

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
                       , ConditionalObservables = tresdata.model().ConditionalObservables()
                       , ExternalConstraints = lifetimeParams.externalConstraints() + tresdata.externalConstraints() + basisCoefficients.externalConstraints
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
#momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))
#These are the relevant terms as found with MinSignificance>3
momindices = [(0,0,0),(0,2,0),(0,2,2),(2,0,0)]

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/data/bfys/dveijk/DataJpsiPhi/2012/effmoments_tcut_%s.txt'%(str(tmincut)))
#Check the significant terms
#eff.read('/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/effmoments_tcut_%s.txt'%(str(tmincut)),MinSignificance = 3)
eff.Print()

#Build Angular acceptance corrected PDF
sig_t_angles = eff * sig_t_angles

##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_20bins')
#acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/data/bfys/dveijk/DataJpsiPhi/2012/propertimeacceptance.root', Histogram = 'timeacceptancehisto')
sig_t_angles = acceptance * sig_t_angles

##################
### Build PDFs ###
##################

############
# BKG COMP #
############
sidebanddata =      data.reduce(CutRange = 'leftsideband' )
sidebanddata.append(data.reduce(CutRange = 'rightsideband'))

nbkg = 10500
background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf(), { eta_os: None, iTag_os : None }), Yield = ( nbkg, 0.9, 1.1*nbkg) )

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
                             , Data  = sidebanddata
                             )

############
# SIG COMP #
############
nsig = 22000
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles ), Yield = ( nsig, 0.9, 1.1*nsig) )

#############
# MASS ONLY #
#############
#masspdf = buildPdf((signal,background), Observables = (m,), Name = 'masspdf')
#masspdf.fitTo(data,**fitOpts)

pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os,eta_os,st)+tuple(angles.angles.itervalues()), Name='fullpdf')
pdf.Print()
classicfitresult = pdf.fitTo(data,**fitOpts)

classicfitresult.writepars('classicfitresult',False)
