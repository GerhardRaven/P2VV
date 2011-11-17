from RooFitWrappers import *

#from ROOT import RooMsgService
#RooMsgService.instance().addStream(RooFit.INFO,RooFit.Topic(RooFit.Optimization))

ws = RooObject( workspace = 'myws' )

from parameterizations import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcpsi', ctheta = 'trctheta', phi = 'trphi' )
t         = RealVar(  't',          Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(-5,14)  )
iTag      = Category( 'tagdecision',Title = 'initial state flavour tag',           Observable = True,  States = { 'B': +1, 'Bbar': -1 } )

observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag ]

from parameterizations import LambdaSqArg_CPParam
from math import pi
#from ROOT import RooUnblindUniform as UnblindUniform
CP = LambdaSqArg_CPParam( lambdaSq  = RealVar( 'lambda^2', Title = 'CP violation param |lambda|^2',  Value = 1)
                        , lambdaArg = RealVar( 'phiCP',    Title = 'CP violation param. phi_s',      Value = -0.2,  MinMax = (-2*pi,2*pi) ) # , Blind = (UnblindUniform,'BsRooBarb',0.2) )
                        )

from parameterizations import  ProdTagNorm_CEvenOdd, Trivial_CEvenOdd
# ANuissance = Trivial_CEvenOdd()
minus = ConstVar('minus',  Value = -1  )
ANuissance = ProdTagNorm_CEvenOdd( AProd   = RealVar(    'AProd',    Title = 'production asymmetry',          Value = 0 )
                                 , ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry',  Value = 0 )
                                 , ANorm   = Product(    'ANorm',   [minus,CP.C],  Title = 'normalization asymmetry' )
                                 )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from parameterizations import JpsiphiAmplitudesLP2011
amplitudes = JpsiphiAmplitudesLP2011()

#### Package from here until the "BTagDecay('name', args)" into a dedicated function...
from parameterizations import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
args = { 'dm'        : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}',  Value = 17.8 )  
       , 'tau'       : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',       Value =  1.5,  MinMax = (1.3, 1.8) )
       , 'dGamma'    : RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}',  Value =  0.05, MinMax = (-0.3,0.3) )
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
       , 'dilution'  : RealVar( 'tagDilution', Title = 'Average Tagging Dilution',      Value = 1 )
       , 'ADilWTag'  : RealVar( 'ADilWTag',    Title = 'dilution/wrong tag asymmetry',  Value = 0 )
       } 

mcpdf = BTagDecay( 'mc_pdf',  args )

if False :
   ws.ws().importClassCode()
   ws.ws().defineSet('observables',','.join( i.GetName() for i in observables ) )
   ws.ws().writeToFile('/tmp/pdf.root')
   from sys import exit
   exit(0)

# update resolution model, and build again...
from parameterizations import ResolutionModelLP2011
args[ 'resolutionModel' ]  = ResolutionModelLP2011( t ).Model

sigpdf = BTagDecay( 'sig_pdf', args )

pdf = mcpdf
#pdf.Print("t")



#l = RooArgSet()
#pdf.branchNodeServerList( l )
#for i in l : i.setAttribute( "CacheAndTrack" )


print 'generating data'
data = pdf.generate( observables , 10000 )

#mcfilename =  '/data/bfys/dveijk/MC/2011/MC2011_UB.root'
#from ROOT import TFile
#MCfile = TFile(mcfilename)
#MCtuple = MCfile.Get('MyTree')
#from ROOT import RooDataSet
#_obs = RooArgSet()
#for i in observables : _obs +=  i._target_() 
#noNAN  = ' && '.join( '%s==%s' % (i.GetName(),i.GetName()) for i in _obs ) 
#print noNAN
#data = RooDataSet('MCdata','MCdata',MCtuple,_obs,noNAN)
#print data.numEntries()


print 'computing efficiency moments'
moms = [ EffMoment( i, 1, pdf, angles.angles.itervalues() ) for v in angles.functions.itervalues() for i in v if i ] 

_bm = lambda i,l,m : EffMoment( P2VVAngleBasis(angles.angles, i,0,l,m,1. ), float(2*l+1)/2, pdf, angles.angles.itervalues() )
moms2  = [ _bm(i,l,m) for i in range(3) for l in range(3) for m in range(-l,l+1) ]
moms2 += [ _bm(i,2,m) for i in range(3,20) for m in [-2,1] ] # these are for the 'infinite' terms in the signal PDF 


computeMoments( data, moms + moms2 )
from pprint import pprint
pprint( [ (m.GetName(), m.coefficient()) for m in moms ] )
pprint( [ (m.GetName(), m.coefficient()) for m in moms2 ] )

#print 'fitting data'
#from ROOT import RooCmdArg
#NumCPU = RooCmdArg( RooFit.NumCPU(8) )
#pdf.fitTo(data, NumCPU, RooFit.Timer(1), RooFit.Minimizer('Minuit2','minimize'))
