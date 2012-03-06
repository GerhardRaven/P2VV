from math import sqrt, pi, cos, sin
from RooFitWrappers import *
from ROOT import RooMsgService
#RooMsgService.instance().addStream(RooFit.DEBUG, RooFit.Topic(RooFit.Integration))
RooMsgService.instance().getStream(1).removeTopic(RooFit.Caching)
RooMsgService.instance().getStream(1).removeTopic(RooFit.Eval)


indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

from P2VVGeneralUtils import numCPU

fitOpts = dict(
    NumCPU = 1
    #NumCPU = numCPU()
    , Timer=1
    , Save = True
    #, Verbose = True
    #, Minimizer = ('Minuit2','minimize')
    )

# define observables
m    = RealVar('mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV/c^{2}', Observable = True, Value = 5368, MinMax = (5200, 5550), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1008,1032), nBins =  16 )
t    = RealVar('time',  Title = 'decay time',    Unit = 'ps',  Observable = True, Value = 0.5, MinMax = (0.3, 14),    nBins =  54 )
#Set the left boundary of sigmat to non-zero to prevent problems with integrations (eg when making plots). Checked the data: no events below 0.007.
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, Value = 0.007, MinMax = (0.007, 0.12),  nBins =  100 )
eta_os  = RealVar('tagomega_os',      Title = 'estimated mistag OS', Observable = True, Value = 0.25, MinMax = (0,0.50001),  nBins =  100)
#The peak at 0.5 seems to be shifted to -2 in the SS eta!
eta_ss  = RealVar('tagomega_ss',      Title = 'estimated mistag SS',          Observable = True, MinMax = (0,0.50001),  nBins =  100)
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

# add 20 bins for caching the normalization integral
for i in [ eta_os, st ] : i.setBins( 100 , 'cache' )

#Read data
from P2VVGeneralUtils import readData
data = readData( '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
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
sig_m = Signal_BMass(     Name = 'sig_m', mass = m )
bkg_m = Background_BMass( Name = 'bkg_m', mass = m )

#Time Resolution Model
from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as DataTimeResolution
tresdata = DataTimeResolution( time = t, timeResSFConstraint = True, sigmat = st, timeResSF = dict( Value = 1.45, MinMax = ( 1., 2. ) ) )

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.67
                                       , dGamma = dict( Name = 'dGamma' , Value = 0.1)
                                       , dM = 17.6
                                       , dMConstraint = True
                                     )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
tagging = TaggingParams(  estWTag = eta_os, p0 = dict( Name = 'wTagP0' ), p1 = dict( Name = 'wTagP1' )
                        , p0Constraint = True, p1Constraint = True )

from P2VVParameterizations.CPVParams import LambdaArg_CPParam, LambdaSqArg_CPParam
CP = LambdaArg_CPParam(phiCP = 0.)

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0 and fs = As2/(1+As2)
from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet
amplitudes = JpsiVPolarSWaveFrac_AmplitudeSet(  A0Mag2 = 0.50, A0Phase = 0
                                              , AperpMag2 = 0.25, AperpPhase = 2.9
                                              , AparPhase = 3.3
                                              , f_S = 0.02
                                              , ASPhase = 2.9
                                             )

from RooFitWrappers import RealCategory
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions
                                                  , amplitudes
                                                  , CP
                                                  , iTag_os
                                                  , tagging['dilution']
                                                  , ['A0','Apar','Aperp','AS'] ) 

basisCoefficients.externalConstraints = tagging.externalConstraints()

sig_t_angles = BDecay( Name      = 'sig_t_angles'
                       , time      = t
                       , dm        = lifetimeParams['dM']
                       , tau       = lifetimeParams['MeanLifetime']
                       , dGamma    = lifetimeParams['dGamma']
                       , resolutionModel = tresdata.model()
                       , coshCoef  = basisCoefficients['cosh']
                       , cosCoef   = basisCoefficients['cos']
                       , sinhCoef  = basisCoefficients['sinh']
                       , sinCoef   = basisCoefficients['sin']
                       , ConditionalObservables = set(tresdata.model().ConditionalObservables()).union( set( [eta_os,] ) )
                       #, ConditionalObservables = set(tresdata.model().ConditionalObservables()).union( set( [iTag_os,eta_os]))
                       , ExternalConstraints = lifetimeParams.externalConstraints() + tresdata.externalConstraints() + basisCoefficients.externalConstraints
                       )

print "*****************************"
print [ i.GetName() for i in sig_t_angles.ConditionalObservables()  ]
print "*****************************"

#####################################
### Angular acceptance correction ###
#####################################
from P2VVGeneralUtils import RealMomentsBuilder
from itertools import chain
#momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,10) for j0 in [1,-2]))
#These are the relevant terms as found with MinSignificance>3
momindices = [(0,0,0),(0,2,0),(0,2,2),(2,0,0)]

eff = RealMomentsBuilder()
#Don't specify pdf and normset here, we're gonna read moments and not calculate any.
eff.appendPYList( angles.angles, momindices)
eff.read('/project/bfys/jleerdam/softDevel/Erasmus/P2VV2/release/test/effmoments_tcut_0.3_Feb.txt')
eff.Print()

#Build Angular acceptance corrected PDF
sig_t_angles = eff * sig_t_angles

##############################
### Proper time acceptance ###
##############################
#from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance
#acceptance = Moriond2012_TimeAcceptance( time = t, Input = '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root', Histogram = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
#sig_t_angles = acceptance * sig_t_angles

####################
### Compose PDFs ###
####################

nEvents = data.numEntries()
signal     = Component('signal',     ( sig_m.pdf(), sig_t_angles, { eta_os: None, st : None } ), Yield = ( 0.67 * nEvents, 0., nEvents) )
background = Component('background', ( bkg_m.pdf(),                                           ), Yield = ( 0.33 * nEvents, 0., nEvents) )

############
### SFIT ###
############
# make sweighted dataset. TODO: use mumu mass as well...
from P2VVGeneralUtils import SData, splot

masspdf = buildPdf((signal,background), Observables = (m,), Name = 'masspdf')
masspdf.fitTo(data,**fitOpts)

massplot = True
if massplot:
    from ROOT import kDashed, kRed, kGreen, TCanvas, TLatex
    from P2VVGeneralUtils import plot
    canvas = TCanvas()
    plot( canvas, m, data, masspdf, components = { 'sig_m' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kGreen )
                                                 , 'bkg_m' : dict( LineStyle = kDashed, LineWidth=3, LineColor = kRed   ) 
                                                 }
                                  , pdfOpts = dict( LineWidth = 3 )
                                  , frameOpts = dict( Title = 'B_{s}#rightarrow J/#psi#phi'
                                                    , TitleOffset = (1.2,'y')
                                                    , Object = ( TLatex(0.55,.8,"#splitline{LHCb preliminary}{#sqrt{s} = 7 TeV, L = 1.03 fb^{-1}}", NDC = True), )
                                                    , Bins=70 ) 
                                  )

splot_m = SData(Pdf = masspdf, Data = data, Name = 'MassSplot')

pdf = buildPdf((signal,), Observables = (t,iTag_os)+tuple(angles.angles.itervalues()), Name='pdf')
#pdf = buildPdf((signal,), Observables = (t,)+tuple(angles.angles.itervalues()), Name='pdf')

fitData = splot_m.data('signal')
assert False
sfitresult = pdf.fitTo( fitData, SumW2Error = False, **fitOpts)
sfitresult.writepars('sfitresult_NoTimeAcc',False)

assert False

#fitset = pdf._var.getParameters(data)
#fitset = pdf.getParameters(data)
#fitset.writeToFile("nominalsfitresult.txt")
#sfitresult.Print()
#fitset.readFromFile("nominalsfitresult.txt")

#Turn this on when fit is fast with NumCPU = numCPU() working!!!
#######################
# Profile likelihoods #
#######################

pllvar = lifetimeParams._dM._var

from ROOT import RooMinuit
#Need to implement conditionalobservables and externalconstraints here

nll = pdf.createNLL(splot_m.data('signal'), NumCPU = numCPU()
                    , Verbose = False
                    #, ConditionalObservables =
                    #, ExternalConstraints = 
                    )
minuit = RooMinuit(nll)
minuit.setProfile(1) #Timer =1
minimize = minuit.migrad() #Save = True
sfitresult = minuit.save()
pll = nll.createProfile(RooArgSet(pllvar)) #This creates a RooProfileLL object!
from ROOT import TCanvas
canvas = TCanvas()
frame = pllvar.frame()
pll.plotOn(frame)
frame.Draw()

########
# PLOT #
########

from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
from P2VVGeneralUtils import plot
orderdict = dict( (i[1].GetName(), i[0]) for i in enumerate([m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]) )

canvas = dict()
for rng in ( None, ) :
    canvas[rng] = TCanvas('%s'%rng)
    obs =  [ o for o in pdf.Observables() if hasattr(o,'frame') ]
    from P2VVGeneralUtils import Sorter
    for (p,o) in zip( canvas[rng].pads(len(obs)), sorted(obs, key = Sorter(orderdict)) ) :
        dataOpts = dict( CutRange =        rng ) if rng else dict()
        pdfOpts  = dict( ProjectionRange = rng ) if rng else dict()
        from P2VVGeneralUtils import plot
        from ROOT import RooArgSet
        pdfOpts[ 'ProjWData' ] = ( RooArgSet(st._var, eta_os._var), data, True )
        plot( p, o, splot_m.data('signal'), pdf
              , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
              , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
#Error for events with negative weight for negative log
#              , logy = ( o == t )
              )
