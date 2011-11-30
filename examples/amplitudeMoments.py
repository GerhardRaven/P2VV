from math import pi, sin, cos, sqrt

# job parameters
generateData = True
nEvents = 10000

dataSetName = 'JpsiKstarData'
dataSetFile = 'amplitudeMoments.root'
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


###########################################################################################################################################

# import RooFit wrappers and load P2VV library
from RooFitWrappers import *
from P2VV import loadP2VVLib, setRooFitOutput
loadP2VVLib()
setRooFitOutput()

# workspace
ws = RooObject(workspace = 'ws')

# variables
from parameterizations import JpsiphiHelicityAngles
angles = JpsiphiHelicityAngles(cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi')

observables = [angle for angle in angles.angles.itervalues()]

# transversity amplitudes
from parameterizations import JpsiVCarthesianAmplitudes
transAmps = JpsiVCarthesianAmplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                      , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                      , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                      , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                      , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                      , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                     )

# build angular PDF
from parameterizations import Amplitudes_AngularPdfTerms
pdfTerms = Amplitudes_AngularPdfTerms(AmpNames = [ 'A0', 'Apar', 'Aperp' ], Amplitudes = transAmps, AngFunctions = angles.functions)
pdf = pdfTerms.buildSumPdf('AngularPDF')


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
pdf.fitTo(data, NumCPU = 2, Timer = 1)

