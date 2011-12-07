from math import pi, sin, cos, sqrt

# job parameters
generateData = False
nEvents = 10000

dataSetName = 'JpsiphiData'
dataSetFile = 'JvLFit.root'
#dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/testSample.root'
NTuple = False

# values of transversity amplitudes
A0Mag2Val    =  0.4
A0PhVal      =  0.
AparMag2Val  =  0.3
AparPhVal    =  2.4
AperpMag2Val =  0.3
AperpPhVal   = -0.79
ASMag2Val    =  0.3
ASPhVal      =  2.4

# values of CP violation parameters
carthLambdaCP = False
phiCPVal      = -pi / 4.
lambdaCPSqVal = 0.6

# values of B lifetime parameters
GammaVal  = 0.68
dGammaVal = 0.05
dmVal     = 17.8

# values of asymmetries
AProdVal       =  0.4
ANormVal       = -(1. - lambdaCPSqVal) / (1. + lambdaCPSqVal)
ATagEffVal     = -0.5
nuissanceAsyms = False


###########################################################################################################################################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# constants
zero = ConstVar('zero', Value = 0.)

# angular functions
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles
angleFuncs = JpsiphiHelicityAngles(cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi')

# variables in PDF
time = RealVar('t',           Title = 'Decay time', Unit = 'ps',   Observable = True, Value = 0., MinMax = (-0.5, 5.))
iTag = Category('tagInitial', Title = 'Initial state flavour tag', Observable = True, States = {'B': +1, 'Bbar': -1})
observables = [angle for angle in angleFuncs.angles.itervalues()] + [time, iTag]

# tagging
wTag    = RealVar('wTag',    Title = 'Wrong tag fraction B',      Value = 0.1, MinMax = (0., 0.5))
wTagBar = RealVar('wTagBar', Title = 'Wrong tag fraction anti-B', Value = 0.2, MinMax = (0., 0.5))

# transversity amplitudes
from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesianAmplitudes
transAmps = JpsiVCarthesianAmplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                      , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                      , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                      , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                      , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                      , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                     )

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
  from P2VVParameterizations.CPVParams import LambdaCarth_CPParam
  lambdaCP = LambdaCarth_CPParam( ReLambdaCP = sqrt(lambdaCPSqVal) * cos(-phiCPVal), ImLambdaCP = sqrt(lambdaCPSqVal) * sin(-phiCPVal) )

else :
  # polar lambda
  from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam
  lambdaCP = LambdaSqArg_CPParam( lambdaCPSq = lambdaCPSqVal, phiCP = phiCPVal )

# nuissance asymmetries
if nuissanceAsyms :
  # use nuissance asymmetries directly
  from P2VVParameterizations.BBbarAsymmetries import ProdTagNorm_CEvenOdd
  ANuissance = ProdTagNorm_CEvenOdd(AProd = AProdVal, ATagEff = ATagEffVal, CPParam = lambdaCP)

else :
  # use average and even coefficients
  from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
  ANuissance = Coefficients_CEvenOdd(  avgCEven = 1.
                                     , avgCOdd = (AProdVal + ANormVal + ATagEffVal + AProdVal * ANormVal * ATagEffVal)
                                                  / (1. + AProdVal * ANormVal + AProdVal * ATagEffVal + ANormVal * ATagEffVal)
                                    )

# coefficients for time functions
from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
timeBasisCoefs = JpsiphiBTagDecayBasisCoefficients(angleFuncs.functions, transAmps, lambdaCP, ['A0','Apar','Aperp','AS']) 

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
from P2VVLoad import RooFitOutput
if generateData :
  print 'fitJpsiV: generating %d events' % nEvents
  data = pdf.generate(observables, nEvents)

  from P2VVGeneralUtils import writeData
  writeData(dataSetFile, dataSetName, data, NTuple)

else :
  from P2VVGeneralUtils import readData
  data = readData(dataSetFile, dataSetName, NTuple)

# fix values of some parameters
#ANuissance.setConstant('avgCOdd')
#lambdaCP.setConstant('phiCP')
wTag.setConstant()
wTagBar.setConstant()

# fit data
pdf.fitTo(data, ConditionalObservables = [iTag], NumCPU = 12, Timer = 1)

