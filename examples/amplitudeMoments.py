###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData = False
nEvents = 10000

# data parameters
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
## build the angular PDF ##
###########################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# variables
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles
angles = JpsiphiHelicityAngles(cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi')

observables = [angle for angle in angles.angles.itervalues()]

# transversity amplitudes
from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesianAmplitudes
transAmps = JpsiVCarthesianAmplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                      , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                      , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                      , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                      , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                      , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                     )

# build angular PDF
from P2VVParameterizations.AngularPDFs import Amplitudes_AngularPdfTerms
pdfTerms = Amplitudes_AngularPdfTerms(AmpNames = [ 'A0', 'Apar', 'Aperp' ], Amplitudes = transAmps, AngFunctions = angles.functions)
pdf = pdfTerms.buildSumPdf('AngularPDF')


###########################################################################################################################################
## generate/read data, fit, calculate moments, ... ##
#####################################################

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

# fit data
print 'fitJpsiV: fitting %d events' % data.numEntries()
pdf.fitTo(data, NumCPU = 2, Timer = 1)


###########################################################################################################################################
## make some plots ##
#####################

from P2VVLoad import ROOTStyle
from ROOT import TCanvas
anglesCanv = TCanvas('anglesCanv', 'Angles')

for pad in anglesCanv.pads( 3, 1 ) :
  print pad

