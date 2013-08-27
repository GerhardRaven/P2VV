###########################################################################################################################################
## set script parameters ##
###########################

from math import pi, sin, cos, sqrt

# job parameters
nEvents = 50000

generateData    = True
physicsPDF      = True

fitDataOriginal = False
fitDataMoments  = False
fitDataCoefs    = False
computeMoments  = True
makePlots       = True
transAngles     = False

# data parameters
dataSetName = 'JpsiKstarData'
#dataSetName = 'JpsiphiPHSP'
#dataSetName = 'Jpsiphi13144002'
dataSetFile = 'computeMoments.root'
#dataSetFile = '/data/bfys/jleerdam/Bs2Jpsiphi/Diego/jpsiphiTrans.root'
plotsFile = 'computeMoments.ps'

if transAngles : angleNames = ( 'cPsi',    'cTheta',    'Phi_b' ) # angleNames = ( 'cpsi_tr', 'ctheta_tr', 'phi_tr' )
else           : angleNames = ( 'cthetaK', 'cthetal',   'phi'   )

# angular moments
momentsFile = 'JpsiKstarMoments'
#momentsFile = 'JpsiphiPHSPMoments'
#momentsFile = 'Jpsiphi13144002Moments'

if physicsPDF :
    # values of transversity amplitudes
    ampsToUse = [ 'A0', 'Apar', 'Aperp' ]#, 'AS' ]
    A0Mag2Val    =  0.45 # 0.556
    A0PhVal      =  0.
    AparMag2Val  =  0.35 # 0.211
    AparPhVal    =  pi # -2.93
    AperpMag2Val =  0.2 # 0.233
    AperpPhVal   =  pi # 2.91
    ASMag2Val    =  0.
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
if transAngles : angleTitles = ( 'cos(#psi_{tr})',  'cos(#theta_{tr})', '#phi_{tr}' )
else           : angleTitles = ( 'cos(#theta_{K})', 'cos(#theta_{l})',  '#phi'      )
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
from P2VV.RooFitWrappers import *

# workspace
ws = RooObject(workspace = 'ws')

# angular functions
if transAngles : from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngularFunctions
else           : from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles     as AngularFunctions
angleFuncs = AngularFunctions( cpsi = angleNames[0], ctheta = angleNames[1], phi = angleNames[2] )

# variables in PDF
angles      = ( angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] )
observables = list(angles)

# build terms for angular PDF
if physicsPDF :
    # terms with transversity amplitudes
    from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet
    transAmps = JpsiVCarthesian_AmplitudeSet(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                             , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                             , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                             , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                             , ReAS    = sqrt(ASMag2Val    / A0Mag2Val) * cos(ASPhVal)
                                             , ImAS    = sqrt(ASMag2Val    / A0Mag2Val) * sin(ASPhVal)
                                            )

    from P2VV.Parameterizations.AngularPDFs import Amplitudes_AngularPdfTerms
    pdfTerms = Amplitudes_AngularPdfTerms( AmpNames = ampsToUse, Amplitudes = transAmps, AngFunctions = angleFuncs.functions )

else :
    # terms with angular basis functions
    from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
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

from P2VV.Load import RooFitOutput
if generateData :
    # generate data with PDF
    print 'computeMoments: generating %d events' % nEvents
    data = pdf.generate( observables, nEvents )

    from P2VV.Utilities.DataHandling import writeData
    writeData( dataSetFile, dataSetName, data )

else :
    # read data from file
    from P2VV.Utilities.DataHandling import readData
    data = readData( dataSetFile, dataSetName = dataSetName )

if fitDataOriginal :
    # fit data
    print 'computeMoments: fitting %d events' % data.numEntries()
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

from P2VV.Utilities.DataMoments import RealMomentsBuilder
moments = RealMomentsBuilder()
moments.appendPYList( angleFuncs.angles, indices )

if computeMoments :
    # compute moments from data set
    moments.compute(data)
    moments.write(momentsFile, Scale = scale)

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

#for event in range( data.numEntries() ) :
#    varSet = data.get(event)
#    angles[0].setVal( varSet.getRealValue('cthetaK') )
#    angles[1].setVal( varSet.getRealValue('cthetal') )
#    angles[2].setVal( varSet.getRealValue('phi') )
#    if momPDF.getVal() < 0. : print angles[0].getVal(), angles[1].getVal(), angles[2].getVal(), momPDF.getVal()

if fitDataMoments :
    ## make some parameters constant
    #for var in momPDF.getVariables() :
    #    if var.GetName().startswith('C_p2vvab') : var.setConstant(True)

    # fit data
    momPDF.fitTo( data, NumCPU = 2, Timer = 1 )


###########################################################################################################################################
## transform moments into angular observables ##
################################################

# measured values
coefNames = [ 'p2vvab_%d0%d%d' % inds for inds in [ ( 0, 2, 0), ( 2, 0, 0 ), ( 2, 2, 0 ) ] ]
coefNum = [ moments[cName].coefficient() for cName in coefNames ]
covNum  = [ [ 0. for cIt1 in range( len(coefNames) ) ] for cIt0 in range( len(coefNames) ) ]
for cIt, cName in enumerate(coefNames) : covNum[cIt][cIt] = moments[cName].variance()

# scale measured values so that C_000 = 4 * ( |A_0|^2 + |A_par|^2 + |A_perp|^2 ) = 4
for valIt in range( len(coefNum) ) :
    coefNum[valIt] *= 4. * scale
    for valItCov in range( len(coefNum) ) :
        covNum[valIt][valItCov] *= ( 4. * scale )**2

# transformation matrix values
#coef0Num = [ 0., 0., 0. ]
#RNum = [  [ 4. / sqrt(5.), 8. / sqrt(5.)  ]
#        , [ -8.,          -16.            ]
#        , [ 8. / sqrt(5.), -8. / sqrt(5.) ]
#       ]

coef0Num = [ 2. / sqrt(5), -4., -2. / sqrt(5.) ]
RNum = [  [ 6. / sqrt(5.) ]
        , [ -12.          ]
        , [ 6. / sqrt(5.) ]
       ]

# matrices
from ROOT import TMatrixD, TMatrixDSym, TVectorD
RMat    = TMatrixD( len(RNum), len(RNum[0]) )
RTrMat  = TMatrixD( len(RNum[0]), len(RNum) )
coefMat = TMatrixD( len(RNum), 1 )
covMat  = TMatrixDSym( len(RNum) )

for rowIt in range( len(RNum) ):
    for colIt in range( len(RNum[0]) ) : RMat[rowIt][colIt]    = RNum[rowIt][colIt]
    for colIt in range( len(RNum[0]) ) : RTrMat[colIt][rowIt]  = RNum[rowIt][colIt]
    for colIt in range(1)              : coefMat[rowIt][colIt] = coefNum[rowIt] - coef0Num[rowIt]
    for colIt in range( len(RNum) )    : covMat[rowIt][colIt]  = covNum[rowIt][colIt]
 
covInvMat = TMatrixD(covMat).Invert()

AMat = RTrMat * covInvMat * RMat
AInvMat = TMatrixD(AMat).Invert()
BMat = RTrMat * covInvMat
IMat = AInvMat * BMat * coefMat

delMat = coefMat - RMat * IMat
delTrMat = TMatrixD( 1, delMat.GetNrows() )
delTrMat.Transpose(delMat)
ChiSqMat = delTrMat * covInvMat * delMat

print 'computeMoments: angular observables from moments (Chi^2 / #dof = %.2f):' % ( ChiSqMat[0][0] / float( IMat.GetNrows() ) )
for rowIt in range( IMat.GetNrows() ) :
  print '  %d: %.4f +/- %.4f' % ( rowIt, IMat[rowIt][0], sqrt( AInvMat[rowIt][rowIt] ) )
AInvMat.Print()


###########################################################################################################################################
## build a PDF from angular basis functions and do a fit ##
###########################################################

# build new PDF with angular coefficients
from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
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
    from P2VV.Load import LHCbStyle

    # create canvas
    from ROOT import TCanvas
    anglesCanv = TCanvas( 'anglesCanv', 'Angles' )

    # make plots
    from P2VV.Utilities.Plotting import plot
    from ROOT import RooFit, RooCmdArg
    for ( pad, obs, nBins, plotTitle, xTitle ) in zip(  anglesCanv.pads(2, 2)
                                                      , angles
                                                      , numBins
                                                      , tuple( [ angle.GetTitle() for angle in angles ] )
                                                      , angleTitles
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

