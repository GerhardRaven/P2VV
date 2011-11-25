from RooFitWrappers import *

#from ROOT import RooMsgService
#RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject( workspace = 'myws' )

from parameterizations import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t         = RealVar(  't', Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(0,14)  )
iTag      = Category( 'tagdecision' , Title = 'initial state flavour tag',           Observable = True,  States = { 'B': +1, 'Bbar': -1 } ) # TODO: , 'untagged' : 0 } )

observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag ]

from parameterizations import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.04, 'MinMax': (0,999) }, lambdaCPSq = ConstVar('one',Value=1) )
CP._phiCP.Print("V")
CP._lambdaCPSq.Print("V")

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from parameterizations import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = 0, ASPhase = 0 )
#amplitudes.setConstant('.*AS.*',True)

#### Package from here until the "BTagDecay('name', args)" into a dedicated class/function...
from parameterizations import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

from parameterizations import JpsiphiBDecayBasisCoefficients
#basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, iTag,  ['A0','Apar','Aperp','AS'] ) 


from parameterizations import  ProdTagNorm_CEvenOdd, Trivial_CEvenOdd
ANuissance = Trivial_CEvenOdd()
#minus = ConstVar('minus',  Value = -1  )
#ANuissance = ProdTagNorm_CEvenOdd( AProd   = RealVar(    'AProd',    Title = 'production asymmetry',          Value = 0 )
#                                 , ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry',  Value = 0 )
#                                 , ANorm   = Product(    'ANorm',   [minus,CP.C],  Title = 'normalization asymmetry' )
#                                 )

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
args = { 'dm'        : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}',  Value = 17.8 )  
       , 'tau'       : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',       Value =  1.0/0.681,  MinMax = ( 1.3, 1.8) )
       , 'dGamma'    : RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}',  Value =  0.060    ,  MinMax = (-0.3, 0.3) )
       , 'resolutionModel' : ResolutionModel( 'resModel', Type = TruthModel, Observables = [ t ] )
       , 'decayType' : 'SingleSided' 
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
from parameterizations import ResolutionModelLP2011
args[ 'resolutionModel' ]  = ResolutionModelLP2011( t ).Model

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


if False : 
    print 'generating data'
    data = pdf.generate( observables , 10000 )

else  :
    mcfilename =  '/data/bfys/dveijk/MC/2011/MC2011_UB.root'
    mcfilename =  '/data/bfys/graven/ntupleB_MC10Bs2JpsiPhi_Reco10_UpDown_simple_with_MCtime_angles_tag.root'
    mcfilename =  '/data/bfys/graven/aladaan.root'
    from ROOT import TFile
    MCfile = TFile(mcfilename)
    MCtuple = MCfile.Get('MyTree')
    from ROOT import RooDataSet
    _obs = RooArgSet()
    cutvar = [] # [ RealVar(i,MinMax=(-1,2)) for i in [ 'sel','triggeredByUnbiasedHlt1AndHlt2','triggeredByBiasedHlt1AndHlt2'] ]
    for i in observables + cutvar : _obs +=  i._target_() 
    noNAN  = ' && '.join( '%s==%s' % (i.GetName(),i.GetName()) for i in _obs )
    data = RooDataSet('MCdata','MCdata',MCtuple,_obs,noNAN ) # ' && '.join([ noNAN, 'sel>0.5', '( (triggeredByUnbiasedHlt1AndHlt2>0.5) || (triggeredByBiasedHlt1AndHlt2>0.5) )' ]))
    print 'got dataset with %s entries' % data.numEntries()


print 'computing efficiency moments'
moms = [ EffMoment( i, 1, pdf, angles.angles.itervalues() ) for v in angles.functions.itervalues() for i in v if i ] 

_bm = lambda i,l,m : EffMoment( P2VVAngleBasis(angles.angles, i,0,l,m,1. ), float(2*l+1)/2, pdf, angles.angles.itervalues() )
moms2  = [ _bm(i,l,m) for i in range(3) 
                      for l in range(3) 
                      for m in range(-l,l+1) ]
moms2 += [ _bm(i,2,m) for i in range(3,10)  # 3,20)
                      for m in [-2,1] ] # these are for the 'infinite' series in the signal PDF 
moms2 = [ ]

computeMoments( data, moms + moms2 )
from pprint import pprint
from math import sqrt,pi
stsp = 16*sqrt(pi)
pprint( [ (m.GetName(), m.coefficient()/stsp, sqrt(m.variance())/stsp, m.significance() ) for m in moms  ] )
pprint( [ (m.GetName(), m.coefficient()/stsp, sqrt(m.variance())/stsp, m.significance() ) for m in moms2 ] )

### TODO: multiply signal PDF with moms2....

#from parameterizations import buildEff_x_PDF
# [ ( m.basis() , m.coefficient() ) for m in moments if m.significance()>signif]
#eff_pdf = buildEff_x_PDF('eff_pdf',pdf._var,[ ( m.basis(), m.coefficient())  for m in moms2 ] )
#pdf = eff_pdf

if True :
    print 'fitting data'
    from ROOT import RooCmdArg
    NumCPU = RooCmdArg( RooFit.NumCPU(7) )
    pdf.fitTo(data, NumCPU, RooFit.Timer(1)) # , RooFit.Minimizer('Minuit2','minimize'))
