from math import pi, sin, cos, sqrt

# job parameters
generateData = False
nEvents = 10000

dataSetName = 'JpsiphiData'
dataSetFile = 'fitJpsiV.root'
#dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/testSample.root'
NTuple = False

# values of transversity amplitudes
A0Mag2Val    = 0.4
A0PhVal      = 0.
AparMag2Val  = 0.3
AparPhVal    = 2.4
AperpMag2Val = 0.3
AperpPhVal   = -0.79
ASMag2Val    = 0.3
ASPhVal      = 2.4

# values of CP violation parameters
carthLambdaCP = True
phiCPVal      = -pi / 4.
lambdaCPSqVal = 0.6

# values of B lifetime parameters
GammaVal  = 0.68
dGammaVal = 0.05
dmVal     = 17.8


###########################################################################################################################################

# import RooFit wrappers and load P2VV library
from RooFitWrappers import *
from P2VV import loadP2VVLib, setRooFitOutput
loadP2VVLib()
setRooFitOutput()

# workspace
ws = RooObject(workspace = 'ws')

# constants
zero  = ConstVar('zero',  Value =  0.)
one   = ConstVar('one',   Value = +1.)
minus = ConstVar('minus', Value = -1.)

# variables
cthetak   = RealVar('hel_cthetak', Title = 'cosine of kaon polarization angle',   Observable = True, Value = 0., MinMax=(-1., 1.))
cthetal   = RealVar('hel_cthetal', Title = 'cosine of lepton polarization angle', Observable = True, Value = 0., MinMax=(-1., 1.))
phiHel    = RealVar('hel_phi',     Title = 'angle between decay planes',          Observable = True, Value = 0., MinMax=(-pi, pi))
time      = RealVar('t',           Title = 'decay time', Unit = 'ps',             Observable = True, Value = 0., MinMax=(-0.5, 5.))
iTag      = Category('tagInitial', Title = 'initial state flavour tag',           Observable = True, States = {'B': +1, 'Bbar': -1})

helAngles = [cthetak, cthetal, phiHel]
observables = helAngles + [time, iTag]

# angular functions
from parameterizations import JpsiphiTransversityAmplitudesHelicityAngles
angFuncs = JpsiphiTransversityAmplitudesHelicityAngles(cpsi = cthetak, ctheta = cthetal, phi = phiHel)

# transversity amplitudes
ReA0    = RealVar('ReA0',    Title = 'Re(A_0)',    Value = 1.)
ImA0    = RealVar('ImA0',    Title = 'Im(A_0)',    Value = 0.)
ReApar  = RealVar('ReApar',  Title = 'Re(A_par)',  Value = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal),  MinMax = (-1., 1.))
ImApar  = RealVar('ImApar',  Title = 'Im(A_par)',  Value = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal),  MinMax = (-1., 1.))
ReAperp = RealVar('ReAperp', Title = 'Re(A_perp)', Value = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal), MinMax = (-1., 1.))
ImAperp = RealVar('ImAperp', Title = 'Im(A_perp)', Value = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal), MinMax = (-1., 1.))
ReAS    = RealVar('ReAS',    Title = 'Re(A_S)',    Value = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal),    MinMax = (-1., 1.))
ImAS    = RealVar('ImAS',    Title = 'Im(A_S)',    Value = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal),    MinMax = (-1., 1.))

from parameterizations import Carthesian_Amplitude
transAmps = {
    'A0'    : Carthesian_Amplitude('A0',    ReA0,    ImA0,    +1)
  , 'Apar'  : Carthesian_Amplitude('Apar',  ReApar,  ImApar,  +1)
  , 'Aperp' : Carthesian_Amplitude('Aperp', ReAperp, ImAperp, -1)
  , 'AS'    : Carthesian_Amplitude('AS',    ReAS,    ImAS,    -1)
}

# B lifetime
Gamma  = RealVar('Gamma',  Title = 'Gamma',       Unit = 'ps^{-1}', Value = GammaVal,  MinMax = (  0.4,  0.9))
dGamma = RealVar('dGamma', Title = 'delta Gamma', Unit = 'ps^{-1}', Value = dGammaVal, MinMax = (- 0.3,  0.3))
dm     = RealVar('dm',     Title = 'delta m',     Unit = 'ps^{-1}', Value = dmVal,     MinMax = ( 13.,  23.))  

from ROOT import RooGaussModel as GaussModel
timeError = RealVar('BLifetimeError', Title = 'B lifetime error reslution model', Unit = 'ps', Value = 0.05)
resModel  = ResolutionModel('resModel', Type = GaussModel, Observables = [time], Parameters = [zero, timeError])

# CP violation parameters
if carthLambdaCP :
  # carthesian lambda
  from parameterizations import LambdaCarth_CPParam
  lambdaCP = LambdaCarth_CPParam(
      ReLambda = RealVar('ReLambdaCP', Title = 'CP violation param. Re(lambda)',
                         Value = sqrt(lambdaCPSqVal) * cos(-phiCPVal), MinMax = (-2., 2.))
    , ImLambda = RealVar('ImLambdaCP', Title = 'CP violation param. Im(lambda)',
                         Value = sqrt(lambdaCPSqVal) * sin(-phiCPVal), MinMax = (-2., 2.))
  )

else :
  # polar lambda
  from parameterizations import LambdaSqArg_CPParam
  lambdaCP = LambdaSqArg_CPParam(
      lambdaSq  = RealVar('lambdaCPSq', Title = 'CP violation param. lambda^2',
                         Value = lambdaCPSqVal, MinMax = (0., 5.))
    , lambdaArg = RealVar('lambdaCPArg', Title = 'CP violation param. arg(lambda)',
                         Value = -phiCPVal, MinMax = (-2. * pi, 2. * pi))
  )

# tagging
wTag    = RealVar('wTag',    Title = 'wrong tag fraction B',      Value = 0.1, MinMax = (0., 0.5))
wTagBar = RealVar('wTagBar', Title = 'wrong tag fraction anti-B', Value = 0.2, MinMax = (0., 0.5))

# nuissance asymmetries
from parameterizations import ProdTagNorm_CEvenOdd
ANuissance = ProdTagNorm_CEvenOdd(
    AProd   = RealVar('AProd',   Title = 'production asymmetry',         Value =  0.4, MinMax = (-1., 1.))
  , ATagEff = RealVar('ATagEff', Title = 'tagging efficiency asymmetry', Value = -0.5, MinMax = (-1., 1.))
  , ANorm   = Product('ANorm',   [minus, lambdaCP.C],  Title = 'normalization asymmetry')
)

# coefficients for time functions
from parameterizations import JpsiphiBTagDecayBasisCoefficients
timeBasisCoefs = JpsiphiBTagDecayBasisCoefficients(angFuncs, transAmps, lambdaCP, ['A0','Apar','Aperp','AS']) 

# build the B_s -> J/psi phi signal PDF
args = {
    'time'            : time
  , 'iTag'            : iTag
  , 'tau'             : FormulaVar('BMeanLife', '1. / @0', [Gamma], Title = 'mean lifetime')
  , 'dGamma'          : dGamma
  , 'dm'              : dm
  , 'dilution'        : FormulaVar('tagDilution', '1. - @0 - @1',               [wTag, wTagBar], Title = 'average tagging dilution')
  , 'ADilWTag'        : FormulaVar('ADilWTag',    '(@0 - @1) / (1. - @0 - @1)', [wTag, wTagBar], Title = 'dilution/wrong tag asymmetry')
  , 'avgCEven'        : ANuissance['avgCEven'] 
  , 'avgCOdd'         : ANuissance['avgCOdd']
  , 'coshCoef'        : timeBasisCoefs['cosh']
  , 'sinhCoef'        : timeBasisCoefs['sinh']
  , 'cosCoef'         : timeBasisCoefs['cos']
  , 'sinCoef'         : timeBasisCoefs['sin']
  , 'resolutionModel' : resModel
  , 'decayType'       : 'SingleSided' 
}

pdf = BTagDecay('JpsiphiPDF', args)

###########################################################################################################################################

# generate data
if generateData :
  print 'fitJpsiV: generating %d events' % nEvents
  data = pdf.generate(observables, nEvents)

  from P2VV import writeData
  writeData(dataSetFile, dataSetName, data, NTuple)

else :
  from P2VV import readData
  data = readData(dataSetFile, dataSetName, NTuple)

# fit data
print 'fitJpsiV: fitting %d events' % data.numEntries()
pdf.fitTo(data, RooFit.NumCPU(12), RooFit.Timer(1))

