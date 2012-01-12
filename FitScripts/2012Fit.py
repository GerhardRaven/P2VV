from math import sqrt, pi
from RooFitWrappers import *

indices = lambda i,l : ( ( _i, _l, _m ) for _i in range(i) for _l in range(l) for _m in range( -_l, _l + 1 )  )
obj  = RooObject( workspace = 'workspace')

from P2VVGeneralUtils import numCPU
from ROOTDecorators import  ROOTversion as Rv
fitOpts = dict( NumCPU = numCPU() 
              , Timer=1
              , Save = True
              , Verbose = False
              , Optimize = True if Rv[1]<32 else 0 # do NOT optimize in 5.32 or later... ( Optimize = 1 only works on a single CPU, 2 doesn't work at all )
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

#For MC dataset only
bkgcat = Category( 'bkgcat',            Title = 'bkgcat',                 Observable = True, States = { 'bkgcat0': 0, 'bkgcat10': 10 } )

from P2VVGeneralUtils import readData
#Read MC data
MCdata = readData('/data/bfys/dveijk/MC/2012/Bs2JpsiPhi_MC11a_ntupleB_for_fitting_20120109.root'
                  , dataSetName = 'DecayTree'
                  , NTuple = True
                  , observables = [ t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], iTag_os,eta_os, triggerdec,sel,bkgcat,mphi]
                  )

#Read data
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
                 , dataSetName = 'DecayTree'
                 , NTuple = True
                 , observables = [ m, mpsi, mphi, t, st, eta_os, eta_ss, iTag_os, iTag_ss, sel, triggerdec, angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi']]
                 , Rename = 'DecayTree_Renamed'
                 )

#TODO: Fix bug: When Rename is on, data = nil, but the dataset is imported in the ws anyway, so this is quick fix:
data = obj.ws().data('DecayTree_Renamed')

print 'Number of events', data.numEntries()
print 'Number of events', MCdata.numEntries()

data.table(iTag_os).Print('v')
data.table(iTag_ss).Print('v')

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
tresdata = DataTimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tresdata.setConstant('.*')

#Time Resolution Model for MC
from P2VVParameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')
#TODO:
#externalConstraints = list()
#externalConstraints += tres.ExternalConstraints()

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.679
                                       , deltaGamma = 0.060
                                       , deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = True) 
                                       )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import WTag_TaggingParams as TaggingParams
tagging = TaggingParams( wTag = eta_os ) # Constant = False, Constrain = True )
# TODO: add external constraint terms for p0 and p1... (and make p0,p1 non-constant ;-)
#externalConstraints += tagging.ExternalConstraints()

# WARNING: we don't try to describe wtag, so when plotting you must use ProjWData for eta_os !!!
#Need this, because eta_os is conditional observable in signal PDF, the actual shape doesn't matter for fitting and plotting purposes
#eta_os_pdf = UniformPdf( Name = 'eta__os_pdf', Arguments = (eta_os,) )

## # Uniform bkg itag distribution
## from P2VVParameterizations.FlavourTagging import Uniform_Background_Tag
## bkg_tag = Uniform_Background_Tag( Name = 'bkg_tag'
##                                   , tagdecision   = iTag_os
##                                   )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
from P2VVParameterizations.CPVParams import ArgOnly_CPParam
CP = ArgOnly_CPParam( phiCP      = dict( Name = 'phi_s', Value = -0.04, MinMax = (-pi,pi), Constant = False ))

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiPhiAmplitudesWinter2012
amplitudes = JpsiPhiAmplitudesWinter2012( A0Mag2 = 0.60, A0Phase = 0
                                          , AperpMag2 = 0.16, AperpPhase = -0.17 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                          , AparPhase = 2.5
                                          , f_S = dict( Value = 0.0, Constant = False )
                                          , ASPhase = dict( Value = 0.0, Constant = False )
                                          )

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

from RooFitWrappers import BDecay
MC_sig_t_angles = BDecay( Name      = 'MC_sig_t_angles'
                       , time      = t
                       , dm        = lifetimeParams['deltaM'] 
                       , tau       = lifetimeParams['MeanLifetime']
                       , dGamma    = lifetimeParams['deltaGamma'] 
                       , resolutionModel = tres.model()
                       , coshCoef  = basisCoefficients['cosh']
                       , cosCoef   = basisCoefficients['cos']
                       , sinhCoef  = basisCoefficients['sinh']
                       , sinCoef   = basisCoefficients['sin']
#                       , ConditionalObservables = ( eta_os, )
                       #                     , ConditionalObservables = ( eta_os, iTag_os, )
                       )

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
                       #                     , ConditionalObservables = ( eta_os, iTag_os, )
                       )

MCpdf = MC_sig_t_angles

print 'Number of MC events', MCdata.numEntries()
allObs = MCpdf.getObservables( MCdata.get() )
print 'MCobservables:', [ i.GetName() for i in allObs ]
o = MCpdf.getObservables(MCdata.get() )

from P2VVGeneralUtils import RealMomentsBuilder
nset = angles.angles.values()

canomoms = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, MCpdf,nset) for v in angles.functions.itervalues() for i in v if i ) )
#canomoms.compute(MCdata)
canomoms.Print(Scales = [1./(16.*sqrt(pi)),1./(16.*sqrt(pi)),1./(16.*sqrt(pi))])

#nset = MCpdf.getObservables( MCdata.get() )
#for a in angles.angles.itervalues() : nset.remove( a._var )
bfun = lambda i,l,m : P2VVAngleBasis( angles.angles, (i,0,l,m))
from itertools import chain
#momindices = indices(3,3)
momindices = chain(indices(3,3),((i0,2,j0) for i0 in range(3,20) for j0 in [1,-2]))
moments = ( RealEffMoment( bfun(*ind), float(2*ind[0]+1)/2, MCpdf, nset ) for ind in momindices )
eff = RealMomentsBuilder( Moments = moments )
#eff.compute(MCdata)
eff.Print()

#BKG ONLY
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

#bkgpdf.fitTo(sidebanddata,**fitOpts)

#SIG ONLY
sigdata = data.reduce(CutRange = 'signal')
sigdata.Print()

nsig = 20000
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles ), Yield = ( nsig, 0, 2.0*nsig) )
#sig_t_angles depends on eta_os, and is NOT conditional on eta_os, so also ask for eta_os here....
sigpdf = buildPdf((signal,), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name = 'sigpdf')

for i in ['ASPhase','f_S'] :
    amplitudes[i].setVal(0.01)
    #amplitudes[i].setConstant(False)

#sigpdf.fitTo(sigdata,**fitOpts)

#MASS ONLY
signalmass = Component('signalmass',(sig_m.pdf(),), Yield = (nsig,0,2.0*nsig))
bkgmass = Component('bkgmass',(bkg_m.pdf(),), Yield = (nbkg,0,2.0*nbkg))
masspdf = buildPdf((signalmass,bkgmass), Observables = (m,), Name = 'masspdf')

#masspdf.fitTo(data,**fitOpts)

#SIG+BKG
for i in ['ASPhase','f_S'] :
    amplitudes[i].setVal(0.01)
    #amplitudes[i].setConstant(False)

#sig_t_angles depends on eta_os, and is NOT conditional on eta_os, so also ask for eta_os here....
pdf   = buildPdf((signal,background), Observables = (m,t,iTag_os,eta_os)+tuple(angles.angles.itervalues()), Name='fullpdf')

pdf.Print()

pdf.fitTo(data,**fitOpts)

assert False


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
