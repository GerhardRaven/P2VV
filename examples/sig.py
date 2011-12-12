from RooFitWrappers import *

#from ROOT import RooMsgService
#RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject( workspace = 'myws' )

from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t         = RealVar(  't', Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(0,14)  )
iTag      = Category( 'tagdecision' , Title = 'initial state flavour tag',           Observable = True,  States = { 'B': +1, 'Bbar': -1 } ) # TODO: , 'untagged' : 0 } )

observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag ]

from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.04, 'MinMax': (0,9.9) }, lambdaCPSq = ConstVar('one',Value=1) )
#CP._phiCP.Print("V")
#CP._lambdaCPSq.Print("V")

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VVParameterizations.DecayAmplitudes import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = 0, ASPhase = 0 )
#amplitudes.setConstant('.*AS.*',True)

#### Package from here until the "BTagDecay('name', args)" into a dedicated class/function...
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

from P2VVParameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
#basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, iTag,  ['A0','Apar','Aperp','AS'] ) 


from P2VVParameterizations.BBbarAsymmetries import  ProdTagNorm_CEvenOdd, Trivial_CEvenOdd
ANuissance = Trivial_CEvenOdd()
#minus = ConstVar('minus',  Value = -1  )
#ANuissance = ProdTagNorm_CEvenOdd( AProd   = RealVar(    'AProd',    Title = 'production asymmetry',          Value = 0 )
#                                 , ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry',  Value = 0 )
#                                 , ANorm   = Product(    'ANorm',   [minus,CP['C']],  Title = 'normalization asymmetry' )
#                                 )


# now build the actual signal PDF...
from P2VVParameterizations.TimeResolution import Truth_TimeResolution
args = { 'dm'        : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}',  Value = 17.8 )  
       , 'tau'       : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',       Value =  1.0/0.681,  MinMax = ( 1.3, 1.8) )
       , 'dGamma'    : RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}',  Value =  0.060    ,  MinMax = (-0.3, 0.3) )
       , 'resolutionModel' : Truth_TimeResolution(time = t)['model']
       , 'time'      : t
       , 'coshCoef'  : basisCoefficients['cosh']
       , 'cosCoef'   : basisCoefficients['cos']
       , 'sinhCoef'  : basisCoefficients['sinh']
       , 'sinCoef'   : basisCoefficients['sin']
       , 'avgCEven'  : ANuissance['avgCEven'] 
       , 'avgCOdd'   : ANuissance['avgCOdd']
       , 'iTag'      : iTag
       , 'dilution'  : ConstVar('one',Value = 1) # FormulaVar('tagdilution', '1-2*@0', [ tagOmega ] )
       , 'ADilWTag'  : RealVar( 'ADilWTag',    Title = 'dilution/wrong tag asymmetry',  Value = 0 )
       } 

mcpdf = BTagDecay( 'mc_pdf',  args )
#mcpdf = BDecay( 'mc_pdf',  args )

if False :
   ws.ws().importClassCode()
   ws.ws().defineSet('observables',','.join( i.GetName() for i in observables ) )
   ws.ws().writeToFile('/tmp/pdf.root')
   from sys import exit
   exit(0)

# update resolution model, and build again...
from P2VVParameterizations.TimeResolution import LP2011_TimeResolution
args[ 'resolutionModel' ]  = LP2011_TimeResolution(time = t)['model']

#sigpdf = BTagDecay( 'sig_pdf', args )

pdf = mcpdf
#pdf.Print("t")



#l = RooArgSet()
#pdf.branchNodeServerList( l )
#for i in l  : 
#    for m in 'ReReRe','ReImIm','ImReIm','ImImRe','P2VVAngle', 'Re','Im' : # , 'a_':
#        if m in i.GetName() : i.setAttribute( "CacheAndTrack" )
#    if i.GetName()[:2] != 'a_' : i.setAttribute( "CacheAndTrack" ) 
#    if not i.getAttribute('CacheAndTrack')  : print i.GetName()
#    #print i.GetName(), i.getAttribute('CacheAndTrack') 
#


if True : 
    print 'generating data'
    data = pdf.generate( observables , NumEvents = 10000 )
    print 'generated %s events' % data.numEntries()

else  :
    mcfilename =  '/data/bfys/dveijk/MC/2011/MC2011_UB.root'
    mcfilename =  '/data/bfys/graven/ntupleB_MC10Bs2JpsiPhi_Reco10_UpDown_simple_with_MCtime_angles_tag.root'
    mcfilename =  '/data/bfys/graven/aladaan.root'
    mcfilename =  '/tmp/aladaan.root'
    from ROOT import TFile
    MCfile = TFile(mcfilename)
    MCtuple = MCfile.Get('MyTree')
    assert MCtuple
    from ROOT import RooDataSet
    noNAN  = ' && '.join( '%s==%s' % (i.GetName(),i.GetName()) for i in observables )
    data = RooDataSet('MCdata','MCdata',MCtuple, ( i._target_() for i in observables  ),noNAN ) # ' && '.join([ noNAN, 'sel>0.5', '( (triggeredByUnbiasedHlt1AndHlt2>0.5) || (triggeredByBiasedHlt1AndHlt2>0.5) )' ]))
    print 'got dataset with %s entries' % data.numEntries()


print 'computing efficiency moments'
from P2VVGeneralUtils import RealMomentsBuilder
eff = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, pdf, angles.angles.itervalues() ) for v in angles.functions.itervalues() for i in v if i ) )
moms2Indices  = [ ( i, l, m ) for i in range(3)
                              for l in range(3)
                              for m in range(-l,l+1) ]
moms2Indices += [ ( i, 2, m ) for i in range(3,10)  # 3,20)
                              for m in [-2,1] ] # these are for the 'infinite' series in the signal PDF
eff.appendPYList( angles.angles, moms2Indices, PDF = pdf, NormSet = angles.angles.itervalues() )

eff.compute(data)

from math import sqrt,pi
eff.Print( MinSignificance = 0., Names = '.*_ang_.*',        Scale = ( 1. / (16*sqrt(pi)), 1. / (16*sqrt(pi)), 1. ) )
eff.Print( MinSignificance = 3., Names = 'p2vvab.*', Scale = ( 1. / ( 2*sqrt(pi)), 1. / ( 2*sqrt(pi)), 1. ) )

pdf.Print("T")
pdf2 = eff * pdf
pdf2.Print("T")

if False :
    print 'fitting data including efficiency'
    pdf.fitTo(data, NumCPU = 4, Timer = 1 , Minimizer = ('Minuit2','minimize'))


from ROOT import TCanvas
c = TCanvas()
from itertools import chain
for (cc,a) in zip(c.pads(4),chain(angles.angles.itervalues(),[t])) :
    f = a.frame( Bins = 24 )
    data.plotOn(f, MarkerSize = 0.8, MarkerColor = RooFit.kRed )
    pdf.plotOn( f , LineColor = RooFit.kBlack)
    pdf2.plotOn( f , LineColor = RooFit.kBlue)
    f.Draw( pad = cc)
