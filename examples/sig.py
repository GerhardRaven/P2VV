from RooFitWrappers import *
from ROOT import gSystem
gSystem.Load('libP2VV.so')

ws = RooObject()
ws.setWorkspace( RooWorkspace('myworkspace') )

plus  = ConstVar('plus', Value = 1)
minus = ConstVar('minus',  Value = -1  )
one   = plus

from math import pi
cpsiAng   = RealVar(  'helcthetaK', Title = 'cosine of kaon polarization angle',   Observable = True,  MinMax=(-1., 1.))
cthetaAng = RealVar(  'helcthetaL', Title = 'cosine of lepton polarization angle', Observable = True,  MinMax=(-1., 1.))
phiAng    = RealVar(  'helphi',     Title = 'angle between decay planes',          Observable = True,  MinMax=(-pi, pi))
t         = RealVar(  't',          Title = 'decay time', Unit='ps',               Observable = True,  MinMax=(-5,14)  )

unbiased  = Category( 'unbiased',   Title = 'Unbiased trigger?',         Observable = True, States = { 'yes': 1 ,   'no' : 0 } )
decay     = Category( 'decaytype',  Title = 'J/psiX decay type',         Observable = True, States = { 'JpsiKplus' : 10, 'JpsiKmin' : 11, 'JpsiKstar0' : 20, 'JpsiKstar0bar': 21, 'Jpsiphi': 40 } )
iTag      = Category( 'iTag',       Title = 'initial state flavour tag', Observable = True, States = { 'B': +1, 'Bbar': -1 } )

from parameterizations import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( lambdaSq = RealVar( 'lambda^2', Title = 'CP violation param |lambda|^2', Observable = False, Value = 1)
                        , arg      = RealVar( 'phiCP',    Title = 'CP violation param. phi_s',     Observable = False, Value = -0.2,  MinMax = (-2*pi,2*pi) )
                        )

from parameterizations import  ProdTagNorm_CEvenOdd, Trivial_CEvenOdd
# ANuissance = Trivial_CEvenOdd()
ANuissance = ProdTagNorm_CEvenOdd( AProd   = RealVar(    'AProd',    Title = 'production asymmetry',         Observable = False, Value = 0 )
                                 , ATagEff = RealVar(    'ATagEff',  Title = 'tagging efficiency asymmetry', Observable = False, Value = 0 )
                                 , ANorm   = Product(    'ANorm',   [minus,CP.C],  Title = 'normalization asymmetry' )
                                 )

# polar transversity amplitudes -- this is 'internal only'
_A0Mag2    = RealVar('A0Mag2',    Title = '|A0|^2',     Observable = False, Value = 0.556, MinMax = (0., 1.))
_A0Ph      = RealVar('delta0',    Title = 'delta_0',    Observable = False, Value = 0. )
_AperpMag2 = RealVar('AperpMag2', Title = '|A_perp|^2', Observable = False, Value = 0.233, MinMax = ( 0., 1.))
_AperpPh   = RealVar('deltaPerp', Title = 'delta_perp', Observable = False, Value = pi-2.91, MinMax=( -2. * pi, 2. * pi))
_AparMag2  = FormulaVar('AparMag2', '1. - @0 - @1', [_A0Mag2, _AperpMag2],  Title = '|A_par|^2' )
_AparPh    = RealVar('deltaPar',  Title = 'delta_par',  Observable = False, Value = 2.93, MinMax = ( -2. * pi, 2. * pi))
_ASMag2    = RealVar('ASMag2',    Title = '|A_S|^2',    Observable = False, Value = 0.05, MinMax=( 0., 1.))
_ASPh      = RealVar('deltaS',    Title = 'delta_S',    Observable = False, Value = 2.2, MinMax=( -2. * pi, 2. * pi))

# in python 2.7 and later, could use collections.OrderedDict so we can traverse without having to think...
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
args = dict()
args.update( { 'dm'     : RealVar( 'dm',        Title = 'delta m',       Unit = 'ps^{-1}', Observable = False, Value = 17.8 )  
             , 'tau'    : RealVar( 't_sig_tau', Title = 'mean lifetime', Unit = 'ps',      Observable = False, Value =  1.5, MinMax = (1.3, 1.8) )
             , 'dGamma' : RealVar( 'dGamma',    Title = 'dGamma',        Unit = 'ps^{-1}', Observable = False, Value = 0.05, MinMax = (-0.3,0.3) )
             , 'resolutionModel' : ResolutionModel( 'resModel', Type = TruthModel, Observables = [ t ] )
             , 'decayType' : 'SingleSided' } )
args.update( { 'avgCEven' : ANuissance['avgCEven'] 
             , 'avgCOdd'  : ANuissance['avgCOdd'] } )
args.update( { 'time'     : t
             , 'iTag'     : iTag
             , 'dilution' : RealVar( 'tagDilution', Title = 'Average Tagging Dilution',     Observable = False, Value = 1 )
             , 'ADilWTag' : RealVar( 'ADilWTag',    Title = 'dilution/wrong tag asymmetry', Observable = False, Value = 0 )
             } )
for i in ['cosh','sinh','cos','sin' ] : args.update( { '%sCoef'%i: basisCoefficients[i] } )
pdf = BTagDecay( 'sig_pdf', args )
ws.ws().Print('V')

print 'generating data '
data = pdf.generate( [ cpsiAng,cthetaAng,phiAng,t,iTag, ] , 10000 )
print 'fitting data '
pdf.fitTo(data)
