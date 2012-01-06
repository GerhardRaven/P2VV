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
m    = RealVar('m',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True, MinMax = (5200, 5550), nBins =  48
                     ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                  , 'signal'        : ( 5330, 5410 )
                                  , 'rightsideband' : ( 5410, None ) 
                                  } )
mpsi = RealVar('mdau1', Title = 'M(#mu#mu)',     Unit = 'MeV', Observable = True, MinMax = (3030, 3150), nBins =  32 )
mphi = RealVar('mdau2', Title = 'M(KK)',         Unit = 'MeV', Observable = True, MinMax = (1007.46, 1031.46), nBins =  16 )
t    = RealVar('t',  Title = 'decay time',    Unit = 'ps',  Observable = True, MinMax = (0.3, 14),    nBins =  54 )
st   = RealVar('sigmat',Title = '#sigma(t)',     Unit = 'ps',  Observable = True, MinMax = (0.0, 0.15),  nBins =  50 )
eta  = RealVar('tagomega',      Title = 'estimated mistag',          Observable = True, MinMax = (0,0.50001),  nBins =  25)
iTag = Category( 'tagdecision', Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1, 'untagged' : 0 } )
#sel  = Category( 'sel',            Title = 'selection',                 Observable = True, States = { 'good': +1 } )
#unbiased  = Category( 'unbiased',            Title = 'unbiased',                 Observable = True, States = { 'Unbiased': +1, 'NotUnbiased': 0 } )
unbiased  = Category( 'unbiased',            Title = 'unbiased',                 Observable = True, States = { 'Unbiased': +1})
#biased  = Category( 'biased',            Title = 'biased',                 Observable = True, States = { 'Biased': +1, 'NotBiased': 0 } )

from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi   = dict( Name = 'trcospsi',   Title = 'cos(#psi)',        nBins = 24 )
                    , ctheta = dict( Name = 'trcostheta', Title = 'cos(#theta_{tr})', nBins = 24 )
                    , phi    = dict( Name = 'trphi',      Title = '#phi_{tr}',        nBins = 24 ) 
                    )

from P2VVGeneralUtils import readData
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2011/Pass3Version2.root'
                 , dataSetName = 'MyTree'
                 , NTuple = True
#                 , observables = [ iTag,eta, mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], st, sel ]
#                 , observables = [ iTag,eta,mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], unbiased, biased]
                 , observables = [ iTag,eta,mpsi,mphi,m,t,angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'], unbiased]
                 )


## print 'Number of events', data.numEntries()

## data.table(iTag).Print('v')

## print 'TAGDECISION FOR UNBIASED EVENTS ONLY'
## data.table(iTag,'unbiased == 1').Print('v')

## print 'TRIGGER SUMMARY FOR ALL EVENTS'
## data.table(biased).Print('v')
## data.table(unbiased).Print('v')

## # Make SuperCategory from (triggeredByUnbiasedHlt1AndHlt2,triggeredByBiasedHlt1AndHlt2)
## TypeCat = SuperCategory('TypeCat',[biased,unbiased])
## data.table(TypeCat).Print('v')
## fitcat = MappedCategory('fitcat',TypeCat,{"AllUnbiased":["{NotBiased;Unbiased}","{Biased;Unbiased}"],"FullyBiased":["{Biased;NotUnbiased}"]})
## data.table(fitcat).Print('v')

## fitcat = data.addColumn(fitcat._var)#Whaa, hackie...
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

# Decay time pdf
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
tres = TimeResolution(time = t) # TODO: extend _util_parse_mixin so that we can add: , Constant = '.*')
tres.setConstant('.*')
#externalConstraints = list()
#externalConstraints += tres.ExternalConstraints()

from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.676
                                     , deltaGamma = 0.107
                                     , deltaM = dict( Value = 17.8, MinMax = (16.5,18.5), Constant = True) 
                                     )

# define tagging parameter 
from P2VVParameterizations.FlavourTagging import Trivial_PerEventTaggingParams as TaggingParams
tagging = TaggingParams( wTag = eta ) # Constant = False, Constrain = True )
# TODO: add external constraint terms for p0 and p1... (and make p0,p1 non-constant ;-)
#externalConstraints += tagging.ExternalConstraints()

# WARNING: we don't try to describe wtag, so when plotting you must use ProjWData for eta !!!
#Need this, because eta is conditional observable in signal PDF, the actual shape doesn't matter for fitting and plotting purposes
eta_pdf = UniformPdf( Name = 'eta_pdf', Arguments = (eta,) )

# Uniform bkg itag distribution
from P2VVParameterizations.FlavourTagging import Uniform_Background_Tag
bkg_tag = Uniform_Background_Tag( Name = 'bkg_tag'
                                , tagdecision   = iTag
                                )

# itag distribution
## bkg_tag = Trivial_Background_Tag( Name = 'bkg_tag'
##                                 , tagdecision   = iTag
##                                 , bkg_tag_eps   = dict( Name = 'tag_bkg_eff', Value = 0.27 )
##                                 , bkg_tag_delta = dict( Name = 'tag_bkg_delta', Value = 0 , MinMax = (-0.1,0.1) ) 
##                                 )
## sig_tag = Trivial_Background_Tag( Name = 'sig_tag'
##                                 , tagdecision   = iTag
##                                 , bkg_tag_eps   = dict( Name = 'tag_sig_eff', Value = 0.35 )
##                                 , bkg_tag_delta = dict( Name = 'tag_sig_delta', Value = 0, MinMax = (-0.1,0.1) ) 
##                                 )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP      = dict( Name = 'phi_s', Value = 0.139, MinMax = (-pi,pi), Constant = False )
                        , lambdaCPSq = ConstVar( Name = 'one', Value = 1 ) 
                        )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.484, A0Phase = 0
                                    , AperpMag2 = 0.259, AperpPhase = 2.647 # , Constant = True ) # untagged with zero CP has no sensitivity to this phase
                                    , AparPhase = 3.263
#                                    , ASMag2 = dict( Value = 0, Constant = True ) , ASPhase = dict( Value = 0, Constant = True )
                                    , ASMag2 = 0.011 , ASPhase = 2.963 
                                    )
# need to specify order in which to traverse...
from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions
                                                  , amplitudes
                                                  , CP
                                                  , Product('tag',(iTag,tagging['dilution']))
                                                  , ['A0','Apar','Aperp','AS'] ) 

from RooFitWrappers import BDecay
sig_t_angles = BDecay( Name      = 'sig_t_angles'
                       , time      = t
                       , dm        = lifetimeParams['deltaM'] 
                       , tau       = lifetimeParams['MeanLifetime']
                       , dGamma    = lifetimeParams['deltaGamma'] 
                       , resolutionModel = tres.model()
                       , coshCoef  = basisCoefficients['cosh']
                       , cosCoef   = basisCoefficients['cos']
                       , sinhCoef  = basisCoefficients['sinh']
                       , sinCoef   = basisCoefficients['sin']
                       , ConditionalObservables = ( eta, )
                       #                     , ConditionalObservables = ( eta, iTag, )
                       )

### TODO: fit MC version of sig_t_angles on MC data and get efficiency,
###       then multiply sig_t_angles with the efficiency...
###       

#mcpdf = BDecay( ... )
#from P2VVGeneralUtils import RealMomentsBuilder
#eff = RealMomentsBuilder()
#eff.appendPYList( angles.angles, indices(4,2) , PDF = mcpdf, NormSet = mcpdf.getObservables( mcdata.get() ) )
#eff.compute(mcdata)
#
## replace signal pdf with efficiency corrected signal pdf
#sig_t_angles = eff * sig_t_angles
# TODO: verify that after multiplication the new PDF still has the same ConditionalObservables!!!!!!!

#####
from P2VVParameterizations.TimePDFs import LP2011_Background_Time as Background_Time
bkg_t = Background_Time( Name = 'bkg_t', time = t, resolutionModel = tres.model()
                       , t_bkg_fll    = dict( Name = 't_bkg_fll',    Value = 0.243 )
                       , t_bkg_ll_tau = dict( Name = 't_bkg_ll_tau', Value = 0.997, MinMax = (0.5,3.5) )
                       , t_bkg_ml_tau = dict( Name = 't_bkg_ml_tau', Value = 0.144 ) )

# Create components
(ntot,fsig) = ( data.numEntries(), 0.65 )
(nsig,nbkg) = ( fsig*ntot, (1-fsig)*ntot) 
#signal         = Component('signal', ( sig_m.pdf(), sig_t_angles, eta_pdf, sig_tag.pdf()), Yield = ( nsig, 0, 1.1*nsig) )
#bkg_background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf(),  eta_pdf, bkg_tag.pdf()), Yield = ( nbkg, 0, 2.0*nbkg) )
signal         = Component('signal', ( sig_m.pdf(), sig_t_angles, eta_pdf, ), Yield = ( nsig, 0, 1.1*nsig) )
bkg_background = Component('bkg'   , ( bkg_m.pdf(), bkg_t.pdf(),  eta_pdf, bkg_tag), Yield = ( nbkg, 0, 2.0*nbkg) )

sidebanddata = data.reduce(CutRange = 'leftsideband')
rightsidebanddata = data.reduce(CutRange = 'rightsideband')
sidebanddata.append(rightsidebanddata)

# create PDF for angular background
if False :
    # make sweighted dataset using J/psi phi mass
    from P2VVGeneralUtils import createSData
    from P2VVParameterizations.AngularPDFs import SPlot_Moment_Angles
    splot_m = createSData( Name = 'Mass' , Components =  (signal,bkg_background), Observables = (m,), FitOpts = fitOpts, Data = data )
    mompdfBuilder = SPlot_Moment_Angles( angles.angles , splot_m )
    bkg_background += mompdfBuilder.pdf( Component = bkg_background.GetName()
                                       , Indices = [ i for i in indices(3,4) ]
                                       , Name = 'bkg_angles'
                                       , MinSignificance = 0.5
                                       , Scale = sqrt(50.) )
elif True:
    bkg_background += UniformPdf( Name = 'bkg_angles', Arguments = angles.angles.itervalues() )
else :
    bkg_background += HistPdf( Name = 'bkg_angles'
                             , Observables = angles.angles.itervalues()
                             , Binning =  { angles.angles['cpsi'] : 7
                                          , angles.angles['ctheta'] : 5
                                          , angles.angles['phi' ] : 9
                                          }
                             , Data  = sidebanddata
                             )

bkganglepdf = buildPdf((bkg_background,), Observables = tuple(angles.angles.itervalues()), Name = 'bkganglepdf')

# fit & fix iTag parameters
#pdf_itag   = buildPdf((signal,bkg_background), Observables = (m,iTag), Name='pdf_itag')
#pdf_itag.fitTo( data,**fitOpts)
#for p in pdf_itag.Parameters() : p.setConstant( not p.getAttribute('Yield') )

#fitdata = unbiaseddata
fitdata = data
fitdata.numEntries()

masspdf = buildPdf((signal,bkg_background), Observables = (m,), Name = 'm_pdf')
m_result = masspdf.fitTo(fitdata, **fitOpts)
m_result.writepars('massfit',False)

# tagged fit
pdf   = buildPdf((signal,bkg_background ), Observables = (m,t,eta,iTag)+tuple(angles.angles.itervalues()), Name='pdf')
# fitOpts[ 'ExternalConstraints' ] = externalConstraints
c_pdf = pdf.fitTo(fitdata, **fitOpts)



from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
from P2VVGeneralUtils import plot

#Old bkganglecanvas
#for (p,o) in zip(canvas.pads(3),angles.angles.itervalues()):
    #p.cd()
    #f = o.frame()
    #sidebanddata.plotOn(f)
    ##bkg_background[angles.angles.itervalues()].dataHist().plotOn(f)
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
    #       and remove the spike at 0.5 to take into account its correlation to iTag = 0!!!
    pdfOpts['ProjWData'] = ( RooArgSet( eta._var ),  data, True ) 
    obs =  [ o for o in pdf.Observables() if hasattr(o,'frame') ]
    for (p,o) in zip( canvas[rng].pads(len(obs)), obs ) :
        plot( p, o, data, pdf, components = { 'signal*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                            , 'bkg*'     : dict( LineColor = kBlue,  LineStyle = kDashed )
                                            }
                             , dataOpts = dict( MarkerSize = 0.8, MarkerColor = kBlack, **dataOpts )
                             , pdfOpts  = dict( LineWidth = 2, **pdfOpts )
                             , logy = ( o == t )
                             )
