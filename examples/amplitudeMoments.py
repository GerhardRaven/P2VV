###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
nEvents = 50000

generateData    = True
physicsPDF      = False

fitDataOriginal = True
fitDataMoments  = True
fitDataCoefs    = True
computeMoments  = True
makePlots       = True

# data parameters
dataSetName = 'JpsiKstarData'
dataSetFile = 'amplitudeMoments_hel.root'
#dataSetFile = 'amplitudeMoments_trans.root'
plotsFile = 'amplitudeMoments.ps'

# angular moments
momentsFile = 'JpsiKstarMoments'

if physicsPDF :
    # values of transversity amplitudes
    ampsToUse = [ 'A0', 'Apar', 'Aperp' ]#, 'AS' ]
    A0Mag2Val    =  0.45
    A0PhVal      =  0.
    AparMag2Val  =  0.35
    AparPhVal    =  pi
    AperpMag2Val =  0.2
    AperpPhVal   =  0.5 * pi
    ASMag2Val    =  0.4
    ASPhVal      =  0.

else :
    # P_i( cos(theta_K) ) * Y_jk( cos(theta_l, phi) ) terms in angular PDF: ( ( i, j, k ), ( value, min. value, max. value ) )
    angPDFParams = [  ( ( 0, 2,  0 ), ( 0.1, 0., 0.2 ) )
                    , ( ( 0, 2, -1 ), ( 0.1, 0., 0.2 ) )
                    , ( ( 0, 2,  2 ), ( 0.1, 0., 0.2 ) )
                    , ( ( 1, 0,  0 ), ( 0.1, 0., 0.2 ) )
                    , ( ( 2, 2,  1 ), ( 0.1, 0., 0.2 ) )
                   ]

# P_i( cos(theta_K) ) * Y_jk( cos(theta_l, phi) ) terms in coefficients PDF: ( ( i, j, k ), ( value, min. value, max. value ) )
coefPDFParams = [  ( ( 0, 2,  0 ), ( 0., -0.6, +0.6 ) )
                 , ( ( 0, 2, -1 ), ( 0., -0.6, +0.6 ) )
                 , ( ( 0, 2,  2 ), ( 0., -0.6, +0.6 ) )
                 , ( ( 1, 0,  0 ), ( 0., -0.6, +0.6 ) )
                 , ( ( 2, 2,  1 ), ( 0., -0.6, +0.6 ) )
                ]

# plot options
angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{l})', '#phi' )
#angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{tr})', '#phi_{tr}' )
numBins    = ( 30, 30, 30 )
lineWidth  = 2
markStyle  = 8
markSize   = 0.4

# moments overall scale
scale = 4. * sqrt(pi)


###########################################################################################################################################
## build the angular PDF, generate/read data and do a fit ##
############################################################

# import RooFit wrappers
from RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles
angleFuncs = JpsiphiHelicityAngles( cpsi = 'cthetaK', ctheta = 'cthetal', phi = 'phi' )
#from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles
#angleFuncs = JpsiphiTransversityAngles( cpsi = 'cpsi_tr', ctheta = 'ctheta_tr', phi = 'phi_tr' )

# variables in PDF
angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = list(angles)

# build terms for angular PDF
if physicsPDF :
    # terms with transversity amplitudes
    from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet
    transAmps = JpsiVCarthesian_AmplitudeSet(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                             , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                             , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                             , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                             , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                             , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                            )

    from P2VVParameterizations.AngularPDFs import Amplitudes_AngularPdfTerms
    pdfTerms = Amplitudes_AngularPdfTerms( AmpNames = ampsToUse, Amplitudes = transAmps, AngFunctions = angleFuncs.functions )

else :
    # terms with angular basis functions
    from P2VVParameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
    cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
    pdfTerms = AngleBasis_AngularPdfTerms(  Angles = angleFuncs.angles
                                          , **dict( (  'C%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                     , {  'Name'    : 'COab%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                        , 'Value'   : term[1][0]
                                                        , 'MinMax'  : term[1][ 1 : 3 ]
                                                        , 'Indices' : term[0]
                                                       }
                                                    ) for term in angPDFParams
                                                  )
                                         )

# build angular PDF
pdf = pdfTerms.buildSumPdf('AngularPDF')

from P2VVLoad import RooFitOutput
if generateData :
    # generate data with PDF
    print 'amplitudeMoments: generating %d events' % nEvents
    data = pdf.generate( observables, nEvents )

    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, data )

else :
    # read data from file
    from P2VVGeneralUtils import readData
    data = readData( dataSetFile, dataSetName = dataSetName )

if fitDataOriginal :
    # fit data
    print 'amplitudeMoments: fitting %d events' % data.numEntries()
    pdf.fitTo( data, NumCPU = 2, Timer = 1 )


###########################################################################################################################################
## compute angular moments and build moments PDFs ##
####################################################

# build angular moment basis functions
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
#indices += [ ( PIndex, 2, YIndex1 ) for PIndex in range( 3, 10 ) for YIndex1 in [ -2, 1 ] ]

# construct moment names strings
names0 = 'p2vvab_00..'
names1 = names0 + '|p2vvab_10..'
names2 = names1 + '|p2vvab_20..'

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
moments.Print( Scale = scale, MinSignificance = 3., Names = names0 )
moments.Print( Scale = scale, MinSignificance = 3., Names = names1 )
moments.Print( Scale = scale, MinSignificance = 3., Names = names2 )
moments.Print( Scale = scale, MinSignificance = 0.                 )

# build new PDFs with angular moments
momPDFTerms0 = moments.buildPDFTerms(MinSignificance = 3., Names = names0, Scale = scale, CoefNamePrefix = 'C0_')
momPDFTerms1 = moments.buildPDFTerms(MinSignificance = 3., Names = names1, Scale = scale, CoefNamePrefix = 'C1_')
momPDFTerms2 = moments.buildPDFTerms(MinSignificance = 3., Names = names2, Scale = scale, CoefNamePrefix = 'C2_')
momPDFTerms  = moments.buildPDFTerms(MinSignificance = 0.                , Scale = scale                        , RangeNumStdDevs = 5.)

momPDF0 = momPDFTerms0.buildSumPdf('angMomentsPDF0')
momPDF1 = momPDFTerms1.buildSumPdf('angMomentsPDF1')
momPDF2 = momPDFTerms2.buildSumPdf('angMomentsPDF2')
momPDF  = momPDFTerms.buildSumPdf('angMomentsPDF')

for event in range( data.numEntries() ) :
    varSet = data.get(event)
    angles[0].setVal( varSet.getRealValue('cthetaK') )
    angles[1].setVal( varSet.getRealValue('cthetal') )
    angles[2].setVal( varSet.getRealValue('phi') )
    if momPDF.getVal() < 0. : print angles[0].getVal(), angles[1].getVal(), angles[2].getVal(), momPDF.getVal()

if fitDataMoments :
    ## make some parameters constant
    #for var in momPDF.getVariables() :
    #    if var.GetName().startswith('C_p2vvab') : var.setConstant(True)

    # fit data
    momPDF.fitTo( data, NumCPU = 2, Timer = 1 )


###########################################################################################################################################
## build a PDF from angular basis functions and do a fit ##
###########################################################

# build new PDF with angular coefficients
from P2VVParameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
coefPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angleFuncs.angles
                                          , **dict( (  'C%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                     , {  'Name'    : 'Cab%d%d%s' % ( term[0][0], term[0][1], cnvrtInd(term[0][2]) )
                                                        , 'Value'   : term[1][0]
                                                        , 'MinMax'  : term[1][ 1 : 3 ]
                                                        , 'Indices' : term[0]
                                                       }
                                                    ) for term in coefPDFParams
                                                  )
                                         )
coefPDF = coefPDFTerms.buildSumPdf('angCoefsPDF')

if fitDataCoefs :
    # fit data
    coefPDF.fitTo( data, NumCPU = 2, Timer = 1 )


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
        plot(  pad, obs, data, pdf, xTitle = xTitle, addPDFs = [ coefPDF, momPDF0, momPDF1, momPDF2, momPDF ]
             , frameOpts   = dict( Bins = nBins, Title = plotTitle )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts     = dict( LineWidth = lineWidth, LineColor = RooFit.kBlack )
             , addPDFsOpts = [  dict( LineWidth = lineWidth, LineColor = RooFit.kGreen + 2 )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kCyan + 2  )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kBlue      )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kMagenta   )
                              , dict( LineWidth = lineWidth, LineColor = RooFit.kRed       )
                             ]
            )

    anglesCanv.Print(plotsFile)

