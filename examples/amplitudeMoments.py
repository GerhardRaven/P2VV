###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
generateData = True
nEvents = 50000

# data parameters
dataSetName = 'JpsiKstarData'
dataSetFile = 'amplitudeMoments.root'
NTuple = False

# values of transversity amplitudes
ampsToUse = [ 'A0', 'Apar', 'Aperp' ] #, 'AS' ]
A0Mag2Val    =  0.4
A0PhVal      =  0.
AparMag2Val  =  0.3
AparPhVal    = -2.4
AperpMag2Val =  0.3
AperpPhVal   = -0.79
ASMag2Val    =  0.3
ASPhVal      =  2.4

# plot options
angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{l})', '#phi' )
numBins    = ( 30, 30, 30 )
lineWidth  = 2
markStyle  = 8
markSize   = 0.4


###########################################################################################################################################
## build the angular PDF ##
###########################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles
angleFuncs = JpsiphiHelicityAngles( cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi' )

# variables in PDF
angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = list(angles)

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
pdfTerms = Amplitudes_AngularPdfTerms( AmpNames = ampsToUse, Amplitudes = transAmps, AngFunctions = angleFuncs.functions )
pdf = pdfTerms.buildSumPdf('AngularPDF')


###########################################################################################################################################
## generate/read data, fit, calculate moments, ... ##
#####################################################

# generate data
from P2VVLoad import RooFitOutput
if generateData :
  print 'amplitudeMoments: generating %d events' % nEvents
  data = pdf.generate( observables, nEvents )

  from P2VVGeneralUtils import writeData
  writeData( dataSetFile, dataSetName, data, NTuple )

else :
  from P2VVGeneralUtils import readData
  data = readData( dataSetFile, dataSetName, NTuple )

# fit data
print 'amplitudeMoments: fitting %d events' % data.numEntries()
pdf.fitTo( data, NumCPU = 2, Timer = 1 )

# calculate angular moments
createMoment = lambda i, l, m : RealMoment( P2VVAngleBasis( angleFuncs.angles, i, 0, l, m, 1. ), float( 2 * l + 1 ) / 2. )
moments  = [ createMoment( i, l, m ) for i in range(3)  for l in range(3) for m in range( -l, l + 1 ) ]
moments += [ createMoment( i, 2, m ) for i in range( 3, 10 ) for m in [ -2, 1 ] ]

computeRealMoments( data, moments )

for mom in moments :
    print mom.GetName(), mom.coefficient(), mom.variance(), mom.significance()


###########################################################################################################################################
## make some plots ##
#####################

# import ROOT plot style
from P2VVLoad import ROOTStyle

# create canvas
from ROOT import TCanvas
anglesCanv = TCanvas('anglesCanv', 'Angles')

# make plots
from P2VVGeneralUtils import plot
from ROOT import RooFit, RooCmdArg
for (pad, obs, nBins, plotTitle, xTitle) in zip(  anglesCanv.pads(2, 2)
                                                , angles
                                                , numBins
                                                , tuple( [ angle.GetTitle() for angle in angles ] )
                                                , angleNames
                                               ) :
    plot(  pad, obs, data, pdf, xTitle = xTitle
         , frameOpts = { 'Bins' : nBins, 'Title' : plotTitle }
         , dataOpts  = { 'MarkerStyle' : markStyle, 'MarkerSize' : markSize }
         , pdfOpts   = { 'LineWidth' : lineWidth }
        )

