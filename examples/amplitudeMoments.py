###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
nEvents = 50000

generateData   = True
fitData        = True
computeMoments = True
makePlots      = True

# data parameters
dataSetName = 'JpsiKstarData'
dataSetFile = 'amplitudeMoments.root'
NTuple = False

# angular moments
momentsFile = 'JpsiKstarMoments'

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

from P2VVLoad import RooFitOutput
if generateData :
    # generate data with PDF
    print 'amplitudeMoments: generating %d events' % nEvents
    data = pdf.generate( observables, nEvents )

    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, data, NTuple )

else :
    # read data from file
    from P2VVGeneralUtils import readData
    data = readData( dataSetFile, dataSetName, NTuple )

if fitData :
    # fit data
    print 'amplitudeMoments: fitting %d events' % data.numEntries()
    pdf.fitTo( data, NumCPU = 2, Timer = 1 )


# build angular moment basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
indices += [ ( PIndex, 2, YIndex1 ) for PIndex in range( 3, 10 ) for YIndex1 in [ -2, 1 ] ]

from P2VVGeneralUtils import RealMomentsBuilder
moments = RealMomentsBuilder()
moments.appendPYList( angleFuncs.angles, indices )

if computeMoments :
    # compute moments from data set
    moments.compute(data)
    moments.write(momentsFile)

else :
    # read moments from file
    moments.read(momentsFile)

# print moments to screen
moments.Print( MinSignificance = 3., Scale = ( 4. * sqrt(pi), 4. * sqrt(pi), 1 ) )


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import ROOT plot style
    from P2VVLoad import ROOTStyle

    # create canvas
    from ROOT import TCanvas
    anglesCanv = TCanvas( 'anglesCanv', 'Angles' )

    # make plots
    from P2VVGeneralUtils import plot
    from ROOT import RooFit, RooCmdArg
    for ( pad, obs, nBins, plotTitle, xTitle ) in zip(  anglesCanv.pads(2, 2)
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

