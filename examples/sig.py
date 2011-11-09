from RooFitWrappers import *
from ROOT import gSystem
gSystem.Load('libP2VV.so')

ws = RooObject()
ws.setWorkspace( RooWorkspace('myworkspace') )

from math import pi
cpsiAng   = RealVar(  'helcthetaK', Title = 'cosine of kaon polarization angle',   Observable = True,  MinMax=(-1., 1.))
cthetaAng = RealVar(  'helcthetaL', Title = 'cosine of lepton polarization angle', Observable = True,  MinMax=(-1., 1.))
phiAng    = RealVar(  'helphi',     Title = 'angle between decay planes',          Observable = True,  MinMax=(-pi, pi))
t         = RealVar(  't',          Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(-5,14)  )
iTag      = Category( 'iTag',       Title = 'initial state flavour tag',           Observable = True,  States = { 'B': +1, 'Bbar': -1 } )
helicityAngles = [ cpsiAng,cthetaAng,phiAng ] # WARNING: order counts!!
observables = helicityAngles + [ t,iTag ]

from parameterizations import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( lambdaSq = RealVar( 'lambda^2', Title = 'CP violation param |lambda|^2',  Value = 1)
                        , arg      = RealVar( 'phiCP',    Title = 'CP violation param. phi_s',      Value = -0.2,  MinMax = (-2*pi,2*pi) )
                        )

from parameterizations import  ProdTagNorm_CEvenOdd, Trivial_CEvenOdd
# ANuissance = Trivial_CEvenOdd()
minus = ConstVar('minus',  Value = -1  )
ANuissance = ProdTagNorm_CEvenOdd( AProd   = RealVar(    'AProd',    Title = 'production asymmetry',          Value = 0 )
                                 , ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry',  Value = 0 )
                                 , ANorm   = Product(    'ANorm',   [minus,CP.C],  Title = 'normalization asymmetry' )
                                 )

# polar transversity amplitudes -- this is 'internal only'
_A0Mag2    = RealVar('A0Mag2',    Title = '|A0|^2',      Value = 0.556,   MinMax = (0., 1.))
_A0Ph      = RealVar('delta0',    Title = 'delta_0',     Value = 0. )
_AperpMag2 = RealVar('AperpMag2', Title = '|A_perp|^2',  Value = 0.233,   MinMax = ( 0., 1.))
_AperpPh   = RealVar('deltaPerp', Title = 'delta_perp',  Value = pi-2.91, MinMax=( -2. * pi, 2. * pi))
_AparMag2  = FormulaVar('AparMag2', '1. - @0 - @1', [_A0Mag2, _AperpMag2],  Title = '|A_par|^2' )
_AparPh    = RealVar('deltaPar',  Title = 'delta_par',   Value = 2.93,    MinMax = ( -2. * pi, 2. * pi))
_ASMag2    = RealVar('ASMag2',    Title = '|A_S|^2',     Value = 0.05,    MinMax=( 0., 1.))
_ASPh      = RealVar('deltaS',    Title = 'delta_S',     Value = 2.2,     MinMax=( -2. * pi, 2. * pi))

from parameterizations import Polar2_Amplitude
Amplitudes = { 'A0'    : Polar2_Amplitude( 'A0',    _A0Mag2,    _A0Ph,    +1 )
             , 'Apar'  : Polar2_Amplitude( 'Apar',  _AparMag2,  _AparPh,  +1 )
             , 'Aperp' : Polar2_Amplitude( 'Aperp', _AperpMag2, _AperpPh, -1 )
             , 'AS'    : Polar2_Amplitude( 'AS',    _ASMag2,    _ASPh,    -1 )
             }

# using transversity amplitudes and helicity angles
from parameterizations import JpsiphiTransversityAmplitudesHelicityAngles
angFuncs = JpsiphiTransversityAmplitudesHelicityAngles( cpsi = cpsiAng, ctheta = cthetaAng, phi = phiAng )

from parameterizations import JpsiphiBTagDecayBasisCoefficients
# need to specify order in which to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angFuncs, Amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

# now build the actual signal PDF...
from ROOT import RooTruthModel as TruthModel
pdf = BTagDecay( 'sig_pdf',  { 'dm'        : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}',  Value = 17.8 )  
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
               )
#ws.ws().Print('V')
#ws.ws().writeToFile("pdf.root")

print 'generating data'
data = pdf.generate( observables , 10000 )


print 'computing efficiency moments'
moms = [ EffMoment( i, 1, pdf, helicityAngles ) for v in angFuncs.itervalues() for i in v if i ] 

_bm = lambda i,l,m : EffMoment( P2VVAngleBasis(helicityAngles, i,0,l,m,1. ), float(2*l+1)/2, pdf, helicityAngles )
moms2  = [ _bm(i,l,m) for i in range(3) for l in range(3) for m in range(-l,l+1) ]
moms2 += [ _bm(i,2,m) for i in range(3,20) for m in [-2,1] ] # these are for the 'infinite' terms in the signal PDF 


computeMoments( data, moms + moms2 )

from pprint import pprint
pprint( [ (m.GetName(), m.coefficient()) for m in moms ] )
pprint( [ (m.GetName(), m.coefficient()) for m in moms2 ] )

print 'fitting data'
from ROOT import RooCmdArg
NumCPU = RooCmdArg( RooFit.NumCPU(8) )
pdf.fitTo(data, NumCPU)
