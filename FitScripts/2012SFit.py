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
#tresdata = DataTimeResolution( time = t, timeResSFConstraint = True ) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
#Per event error
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
tresdata = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')

externalConstraints = list()
externalConstraints += tresdata.externalConstraints()

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                     , deltaGamma = dict( Name = 'dGamma'
                                                        , Value = 0.060
                                                        #, Blind = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
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
                                              #, Blind =  ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 ))
                                              )
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
#amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.52, A0Phase = 0
#                                              , AperpMag2 = 0.25, AperpPhase = 2.77 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
#                                              , AparPhase = 3.2
#                                              , f_S_Re = dict( Value = 0.02 * cos(2.7), Constant = False )
#                                              , f_S_Im = dict( Value = 0.02 * sin(2.7), Constant = False )
#                                             )

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
from itertools import chain
momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/user/dveijk/LHCb/P2VV/FreshStart/p2vv/FitScripts/effmoments_tcut_%s.txt'%(tmincut))
eff.Print()

#Build Angular acceptance corrected PDF
sig_t_angles = eff * sig_t_angles

##############################
### Proper time acceptance ###
##############################
from P2VVParameterizations.TimeAcceptance import LP2011_TimeAcceptance
acceptance = LP2011_TimeAcceptance( time = t )
sig_t_angles = acceptance.acceptance() * sig_t_angles

##################
### Build PDFs ###
##################

nsig = 20000
nbkg = 10000
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles ), Yield = ( nsig, 0, 2.0*nsig) )
background     = Component('bkg',    ( bkg_m.pdf(), ),              Yield = ( nbkg, 0, 2.0*nbkg) )

############
### SFIT ###
############
# make sweighted dataset. TODO: using J/psi phi mass
from P2VVGeneralUtils import SData, splot

masspdf = buildPdf((signal,background), Observables = (m,), Name = 'masspdf')
masspdf.fitTo(data,**fitOpts)
for p in masspdf.Parameters() : p.setConstant( not p.getAttribute('Yield') )
splot_m = SData(Pdf = masspdf, Data = data, Name = 'MassSplot')

sfitsignal = Component('sfitsignal', ( sig_t_angles, ), Yield = ( nsig, 0, 2.0*nsig) )

pereventerror = True
if pereventerror:
    sfitpdf = buildPdf((sfitsignal,), Observables = (t,iTag_os,eta_os,st)+tuple(angles.angles.itervalues()), Name='sfitpdf')
    sfitpdf.Print()
    #Don't add externalconstraints to fitOpts, otherwise fits for splots might go wrong, you don't want to constrain mass fits!
    sfitresult = sfitpdf.fitTo( splot_m.data('signal'), ExternalConstraints = externalConstraints, SumW2Error = True, ConditionalObservables = [st], **fitOpts)

else:
    sfitpdf = buildPdf((sfitsignal,), Observables = (t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name='sfitpdf')
    sfitpdf.Print()
    #Don't add externalconstraints to fitOpts, otherwise fits for splots might go wrong, you don't want to constrain mass fits!
    sfitresult = sfitpdf.fitTo(S_sigdata,ExternalConstraints = externalConstraints, SumW2Error = True, **fitOpts)

sfitresult.writepars('sfitresult_NOTimeAcc',False)

