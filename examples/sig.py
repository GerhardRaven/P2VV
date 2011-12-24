from RooFitWrappers import *

ws = RooObject( workspace = 'myws' )

from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t         = RealVar(  't', Title = 'decay time', Unit='ps',                  Observable = True,  MinMax=(0,14)  )
iTag      = Category( 'tagdecision' , Title = 'initial state flavour tag',   Observable = True,  States = { 'B': +1, 'Bbar': -1 } ) # , 'untagged' : 0 } )
eta       = RealVar(   'eta', Title = 'estimated mis tag', Observable = False, Constant = True, Value = 0.3, MinMax=(0,0.5) )
mass = RealVar('m',Observable=True,Unit='MeV/c^2',MinMax=(5200,5550))
observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag, mass ]

for i in angles.angles.itervalues() : i.setBins(16)
t.setBins(48)
mass.setBins(48)
mass.setRange('leftsideband', (mass.getMin(),5330) )
mass.setRange('signal',(5330,5410) )
mass.setRange('rightsideband',(5410,mass.getMax()) )

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.4, 'MinMax': (-3.2,3.2) }, lambdaCPSq = ConstVar(Name='one',Value=1) )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = { 'Value' : 0, 'Constant': True} , ASPhase = { 'Value': 0, 'Constant':True } )
#amplitudes.setConstant('.*AS.*',True)

#### Package from here until the "BTagDecay('name', args)" into a dedicated class/function...
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

#from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
#basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, iTag,  ['A0','Apar','Aperp','AS'] ) 

from P2VVParameterizations.FlavourTagging import Trivial_TaggingParams
taggingParams = Trivial_TaggingParams( wTag = eta ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )

# now build the actual signal PDF...
from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.68, deltaGamma = 0.05, deltaM = dict( Value = 17.8, MinMax = (16,19), Constant = True) )

from P2VVParameterizations.TimeResolution import Truth_TimeResolution
args = { 'time'      : t
       , 'iTag'      : iTag
       , 'dm'        : lifetimeParams['deltaM'] 
       , 'tau'       : lifetimeParams['MeanLifetime']
       , 'dGamma'    : lifetimeParams['deltaGamma'] 
       , 'resolutionModel' : Truth_TimeResolution(time = t)['model']
       , 'coshCoef'  : basisCoefficients['cosh']
       , 'cosCoef'   : basisCoefficients['cos']
       , 'sinhCoef'  : basisCoefficients['sinh']
       , 'sinCoef'   : basisCoefficients['sin']
       , 'avgCEven'  : taggingParams['avgCEven'] 
       , 'avgCOdd'   : taggingParams['avgCOdd']
       , 'dilution'  : taggingParams['dilution']
       , 'ADilWTag'  : taggingParams['ADilWTag']
       }

# TODO: should be able to write BTagDecay('mypdf', **lifetimeParams.BTagDecay() + **basisCoefficients.BTagDecay() + **taggingParams.BTagDecay() )
mcpdf = BTagDecay( 'mc_pdf', **args  )
#mcpdf = BDecay( 'mc_pdf',  args )

# update resolution model, and build again...
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
args[ 'resolutionModel' ]  = LP2011_TimeResolution(time = t)['model']

#sigpdf = BTagDecay( 'sig_pdf', args )

pdf = mcpdf

from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass
signal = Component('sig',(  LP2011_Signal_Mass( mass = mass ).pdf(),  pdf ), Yield = (1000,0,15000) )

from P2VVParameterizations.MassPDFs import LP2011_Background_Mass
from P2VVParameterizations.TimePDFs import LP2011_Background_Time
from P2VVParameterizations.FlavourTagging import Trivial_Background_Tag
from P2VVParameterizations.AngularPDFs import Uniform_Angles
bkg  = Component('bkg',(  LP2011_Background_Mass( mass = mass ).pdf()
                       ,  LP2011_Background_Time( time = t , resolutionModel = LP2011_TimeResolution(time = t)['model']).pdf()
                       ,  Trivial_Background_Tag( tagdecision = iTag, bkg_tag_delta = 0.3 ).pdf()
                       ,  Uniform_Angles( angles = angles.angles ).pdf()
                       ), Yield = (4000,1000,15000) )

pdf = buildPdf( (signal,bkg), Observables = observables,  Name = 'jointpdf' )



if True : 
    print 'generating data'
    data = pdf.generate( observables ) # , NumEvents = 10000 )
    print 'generated %s events' % data.numEntries()
else  :
    mcfilename =  '/data/bfys/dveijk/MC/2011/MC2011_UB.root'
    mcfilename =  '/data/bfys/graven/ntupleB_MC10Bs2JpsiPhi_Reco10_UpDown_simple_with_MCtime_angles_tag.root'
    mcfilename =  '/data/bfys/graven/aladaan.root'
    mcfilename =  '/tmp/aladaan.root'
    from P2VVGeneralUtils import readData
    data = readData( mcfilename, 'MyTree', True, observables )
    print 'got dataset with %s entries' % data.numEntries()


if False :
    print 'computing efficiency moments'
    from P2VVGeneralUtils import RealMomentsBuilder
    # eff = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, pdf, angles.angles.itervalues() ) for v in angles.functions.itervalues() for i in v if i ) )
    eff = RealMomentsBuilder()
    indices  = [ ( i, l, m ) for i in range(3)
                             for l in range(3)
                             for m in range(-l,l+1) ]
    indices += [ ( i, 2, m ) for i in range(3,10)  # 3,20)
                             for m in [-2,1] ] # these are for the 'infinite' series in the signal PDF
    eff.appendPYList( angles.angles, indices, PDF = pdf, NormSet = angles.angles.itervalues() )
    eff.compute(data)

    from math import sqrt,pi
    #eff.Print( MinSignificance = 0., Names = '.*_ang_.*',   Scales = ( 1. / (16*sqrt(pi)), 1. / (16*sqrt(pi)), 1. ) )
    eff.Print( MinSignificance = 0., Names = 'p2vvab.*',    Scales = ( 1. / ( 2*sqrt(pi)), 1. / ( 2*sqrt(pi)), 1. ) )

    #pdf.Print("T")
    pdf2 = eff * pdf
    #pdf2.Print("T")

    # create generic PDF to describe the angular distribution
    moms = RealMomentsBuilder()
    moms.appendPYList( angles.angles, indices )
    moms.compute(data)
    mom_pdf = moms.createPDF( Name = 'mom_pdf' )

from P2VVGeneralUtils import numCPU
pdf.fitTo(data, NumCPU = numCPU(), Timer = 1 , Minimizer = ('Minuit2','minimize'), Optimize = 0)

from ROOT import TCanvas
plots = [ o for o in observables if hasattr(o,'frame') ]

# TODO: plot signal and sideband seperately!!!
c = dict()
for rng in ( None, 'signal', 'leftsideband','rightsideband','leftsideband,rightsideband' ) :
    c[rng] = TCanvas('%s'%rng)
    for (cc,(lab,ind)) in zip(c[rng].pads(1,len(iTag.states().items())),iTag.states().items()) :
        dataCuts = dict( Cut = '%s == %s' % ( iTag.GetName(), ind ), **(dict( CutRange = rng ) if rng else {}) )
        pdfCuts  = dict( Slice = ( iTag, lab ), **(dict( ProjectionRange = rng ) if rng else {}) )
        for (ccc,o) in zip(cc.pads(len(plots)),plots) :
            from P2VVGeneralUtils import plot
            plot( ccc, o, data, pdf, components = { 'sig*' : { 'LineColor' : RooFit.kGreen, 'LineStyle' : RooFit.kDashed }
                                                  , 'bkg*' : { 'LineColor' : RooFit.kRed, 'LineStyle' : RooFit.kDashed }
                                                  }
                                   , dataOpts = dict( MarkerSize = 0.8, MarkerColor = RooFit.kBlack, **dataCuts )
                                   , pdfOpts  = dict( LineWidth = 2 , **pdfCuts )
                                   , logy = o == t
                                   )
        #if o in mom_pdf.getObservables( data )  : mom_pdf.plotOn( f, LineColor = RooFit.kGreen)
