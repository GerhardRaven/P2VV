from math import pi, sin, cos, sqrt

# job parameters
generateData = False
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

# constants
zero = ConstVar('zero', Value = 0.)

# variables
from parameterizations import JpsiphiHelicityAngles as HelAngles
angles = HelAngles(cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi')

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

