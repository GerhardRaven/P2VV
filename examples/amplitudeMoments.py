###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
readData = False
nEvents = 50000

fitDataOriginal = True
fitDataMoments  = False
fitDataCoefs    = True
makePlots       = True

angleNames = ( 'cthetaK', 'cthetal', 'phi' )

# data parameters
dataSetName = 'angleData'
dataSetFile = 'amplitudeMoments.root'
plotsFile   = 'amplitudeMoments.ps'
momentsFile = 'angularMoments'

# values of transversity amplitudes
ampsToUse = [ 'A0', 'Apar', 'Aperp', 'AS' ]
A0Mag2Val    =  0.4
A0PhVal      =  0.
AparMag2Val  =  0.3
AparPhVal    =  0.8
AperpMag2Val =  0.3
AperpPhVal   = -0.4
ASMag2Val    =  0.01
ASPhVal      =  2.4

# P_i( cos(theta_K) ) * Y_jk( cos(theta_l, phi) ) terms in coefficients PDF: ( ( i, j, k ), ( value, error, range [error] ) )
coefPDFParams = [  ( ( 0, 2,  0 ), -0.039, 0.004, 3. )
                 , ( ( 1, 0,  0 ), -0.158, 0.008, 3. )
                 , ( ( 1, 2,  0 ), +0.071, 0.007, 3. )
                 , ( ( 2, 0,  0 ), +0.165, 0.010, 3. )
                 , ( ( 2, 2,  0 ), -0.498, 0.009, 3. )
                ]

# moments overall scale
scale = 4. * sqrt(pi)

# fit options
fitOpts = dict( NumCPU = 2, Optimize = 2, Timer = True, Minimizer = 'Minuit2' )

# plot options
angleTitles = ( 'cos#kern[0.1]{#theta_{K}}', 'cos#kern[0.1]{#theta_{l}}', '#varphi [rad]' )
numBins     = ( 30, 30, 30 )
lineWidth   = 2
markStyle   = 8
markSize    = 0.4


###########################################################################################################################################
## build the angular PDF, generate/read data and do a fit ##
############################################################

# set RooFit output
from P2VV.Load import RooFitOutput

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'ampMomsWS' ).ws()

# angular functions
from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngularFunctions
angleFuncs = AngularFunctions( cpsi = angleNames[0], ctheta = angleNames[1], phi = angleNames[2] )

# variables in PDF
angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = list(angles)

# build terms for angular PDF
from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet
transAmps = JpsiVCarthesian_AmplitudeSet(  ReApar  = sqrt( AparMag2Val  / A0Mag2Val ) * cos( AparPhVal  - A0PhVal )
                                         , ImApar  = sqrt( AparMag2Val  / A0Mag2Val ) * sin( AparPhVal  - A0PhVal )
                                         , ReAperp = sqrt( AperpMag2Val / A0Mag2Val ) * cos( AperpPhVal - A0PhVal )
                                         , ImAperp = sqrt( AperpMag2Val / A0Mag2Val ) * sin( AperpPhVal - A0PhVal )
                                         , ReAS    = sqrt( ASMag2Val    / A0Mag2Val ) * cos( ASPhVal    - A0PhVal )
                                         , ImAS    = sqrt( ASMag2Val    / A0Mag2Val ) * sin( ASPhVal    - A0PhVal )
                                        )

from P2VV.Parameterizations.AngularPDFs import Amplitudes_AngularPdfTerms
pdfTerms = Amplitudes_AngularPdfTerms( AmpNames = ampsToUse, Amplitudes = transAmps, AngFunctions = angleFuncs.functions )

# build angular PDF
pdf = pdfTerms.buildSumPdf('AngularPDF')

if readData :
    # read data
    from P2VV.GeneralUtils import readData
    data = readData( dataSetFile, dataSetName = dataSetName )
    print 'amplitudeMoments: read %d events' % data.sumEntries()

else :
    # generate data with PDF
    print 'amplitudeMoments: generating %d events' % nEvents
    data = pdf.generate( observables, nEvents )

    from P2VV.GeneralUtils import writeData
    writeData( dataSetFile, dataSetName, data )

if fitDataOriginal :
    # fit amplitudes PDF to data
    print 'amplitudeMoments: fitting %d events' % data.numEntries()
    pdf.fitTo( data, **fitOpts )


###########################################################################################################################################
## compute angular moments and build moments PDFs ##
####################################################

# build angular moment basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(4) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]

# construct moment names strings
names0 = 'p2vvab_0000'
names1 = names0 + '|p2vvab_001.|p2vvab_100.|p2vvab_101.'

from P2VV.GeneralUtils import RealMomentsBuilder
moments = RealMomentsBuilder()
moments.appendPYList( angleFuncs.angles, indices )

# compute moments from data set
moments.compute(data)
moments.write( momentsFile, Scale = scale )

# print moments to screen
moments.Print( Scale = scale, MinSignificance = 0.                 )
moments.Print( Scale = scale, MinSignificance = 3., Names = names0 )
moments.Print( Scale = scale, MinSignificance = 3., Names = names1 )

# build new PDFs with angular moments
momPDFTerms  = moments.buildPDFTerms( MinSignificance = 0.                , Scale = scale                        , RangeNumStdDevs = 5. )
momPDFTerms0 = moments.buildPDFTerms( MinSignificance = 3., Names = names0, Scale = scale, CoefNamePrefix = 'C0_'                       )
momPDFTerms1 = moments.buildPDFTerms( MinSignificance = 3., Names = names1, Scale = scale, CoefNamePrefix = 'C1_'                       )

momPDF  = momPDFTerms.buildSumPdf('angMomentsPDF')
momPDF0 = momPDFTerms0.buildSumPdf('angMomentsPDF0')
momPDF1 = momPDFTerms1.buildSumPdf('angMomentsPDF1')

if fitDataMoments :
    # fit angular moments PDF to data
    momPDF.fitTo( data, **fitOpts )


###########################################################################################################################################
## build a PDF from angular basis functions and do a fit ##
###########################################################

# build new PDF with angular coefficients
from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
coefPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angleFuncs.angles
                                          , **dict( (  'C%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                     , {  'Name'    : 'Cab%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                        , 'Value'   : term[1]
                                                        , 'Error'   : term[2]
                                                        , 'MinMax'  : ( term[1] - term[3] * term[2], term[1] + term[3] * term[2] )
                                                        , 'Indices' : term[0]
                                                       }
                                                    ) for term in coefPDFParams
                                                  )
                                         )
coefPDF = coefPDFTerms.buildSumPdf('angCoefsPDF')

if fitDataCoefs :
    # fit data
    fitResult = coefPDF.fitTo( data, Save = True, **fitOpts )
    fitResult.PrintSpecial( text = True )
    fitResult.correlationMatrix().Print()


###########################################################################################################################################
## make some plots ##
#####################

if makePlots :
    # import ROOT plot style
    from P2VV.Load import LHCbStyle

    # create canvas
    from ROOT import TCanvas
    anglesCanv = TCanvas( 'anglesCanv', 'Angles' )

    # make plots
    from P2VV.GeneralUtils import plot
    from ROOT import RooFit, RooCmdArg
    for ( pad, obs, nBins, plotTitle, xTitle ) in zip(  anglesCanv.pads(2, 2)
                                                      , angles
                                                      , numBins
                                                      , tuple( [ angle.GetTitle() for angle in angles ] )
                                                      , angleTitles
                                                     ) :
        plot(  pad, obs, data, momPDF0, xTitle = xTitle, addPDFs = [ momPDF1, momPDF, coefPDF, pdf ]
             , frameOpts   = dict( Bins = nBins, Title = plotTitle )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts     = dict( LineWidth = lineWidth, LineColor = RooFit.kBlack )
             , addPDFsOpts = [  dict( LineWidth = lineWidth, LineColor = RooFit.kGreen + 2 )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kMagenta   )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kRed       )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kBlue      )
                             ]
            )

    anglesCanv.Print(plotsFile)
