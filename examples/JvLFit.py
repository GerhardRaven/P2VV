###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as pdfConfig

# job parameters
generateData   = False
doFit          = False

pdfConfig['makePlots']  = False
pdfConfig['SFit']       = False
pdfConfig['blind']      = True
pdfConfig['nominalPdf'] = True

plotsFile = 'JvLSFit.ps' if pdfConfig['SFit'] else 'JvLCFit.ps'
pdfConfig['angEffMomentsFile'] = angEffMomentsFile = 'effMomentsTransBasis' if pdfConfig['nominalPdf'] else 'effMomentsHelBasis'

pdfConfig['nTupleName'] = 'DecayTree'
pdfConfig['nTupleFile'] = '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20120120.root'
if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLFit.root'

pdfConfig['timeEffHistFile'] = '/data/bfys/dveijk/DataJpsiPhi/2012/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root'
pdfConfig['timeEffHistName'] = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_20bins'

# PDF options
pdfConfig['components']         = 'all'  # 'all' / 'signal' / 'background'
pdfConfig['transversityAngles'] = False
pdfConfig['bkgAnglePdf']        = ''     # '' / 'histPdf'
pdfConfig['perEventTimeRes']    = False
pdfConfig['multiplyByTimeEff']  = ''     # 'all' / 'signal'

pdfConfig['taggingConditionals'] = 'estWTag'

nEvents = 30000
pdfConfig['signalFraction'] = 0.67
pdfConfig['massRangeBackground'] = False

# transversity amplitudes
pdfConfig['amplitudeParam'] = 'phasesSWaveFrac'

pdfConfig['A0Mag2']    = 0.50
pdfConfig['AperpMag2'] = 0.25
pdfConfig['AparMag2']  = 1. - pdfConfig['A0Mag2'] - pdfConfig['AperpMag2']

pdfConfig['A0Phase']    = 0.
pdfConfig['AperpPhase'] = 3.
pdfConfig['AparPhase']  = 3.

pdfConfig['f_S']     = 0.02
pdfConfig['ASMag2']  = pdfConfig['f_S'] / ( 1. - pdfConfig['f_S'] )
pdfConfig['ASPhase'] = 3.

# CP violation parameters
pdfConfig['carthLambdaCP'] = False
pdfConfig['phiCP']         = 0.2 if pdfConfig['blind'] else 0.
pdfConfig['lambdaCPSq']    = 1.

# B lifetime parameters
pdfConfig['Gamma']      = 0.67
pdfConfig['deltaGamma'] = 0.1
pdfConfig['deltaM']     = 17.6

# asymmetries
pdfConfig['AProd'] = 0.

# fit options
fitOpts = dict(  NumCPU              = 1
               , Timer               = 1
#               , Minos               = False
#               , Hesse               = False
#               , Save                = True
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
if pdfConfig['nominalPdf'] :
    pdfConfig['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                               , ( 'trcostheta', 'cos(#theta_{tr})' )
                               , ( 'trphi',      '#phi_{tr}'        )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos(#theta_{K})' )
                               , ( 'helcosthetaL', 'cos(#theta_{l})' )
                               , ( 'helphi',       '#phi'            )
                              )
angleNames = pdfConfig['angleNames']

numBins      = ( 60, 30, 30, 30 )
numTimeBins  = ( 30, 30, 30 )
numAngleBins = ( 20, 20, 20 )
pdfConfig['numBkgAngleBins'] = ( 5, 7, 9 )

lineWidth = 2
markStyle = 8
markSize  = 0.4


###########################################################################################################################################
## build PDF ##
###############

from P2VVLoad import RooFitOutput

# workspace
from RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace')

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()


###########################################################################################################################################
## generate/read data and fit ##
################################

if generateData :
    # generate data
    print 'JvLFit: generating %d events' % nEvents
    fitData = pdf.generate(obsSetP2VV)

    from P2VVGeneralUtils import writeData
    writeData( dataSetFile, dataSetName, fitData )
elif pdfConfig['SFit'] :
    fitData = pdfBuild['sigSWeightData']
else :
    fitData = pdfBuild['data']

if doFit :
    # fix values of some parameters
    pdfBuild['lambdaCP'].setConstant('lambdaCPSq')
    for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] :
        CEvenOdd.setConstant('avgCEven.*')
        CEvenOdd.setConstant('avgCOdd.*')
    pdfBuild['taggingParams'].setConstant('tagCatCoef.*')
    for coef in pdfBuild['bkgAngCoefs'] : coef.setConstant()

    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if pdfConfig['SFit'] : fitResult = pdf.fitTo( fitData, SumW2Error = True, **fitOpts )
    else    : fitResult = pdf.fitTo( fitData,                    **fitOpts )

    print 120 * '=' + '\n'


###########################################################################################################################################
## make some plots ##
#####################

if pdfConfig['makePlots'] :
    # get variables
    obsSetP2VV = [ pdfBuild['observables'][obs] for obs in [  'time', angleNames[0][0], angleNames[1][0], angleNames[2][0]
                                                            , 'iTag', 'tagCatP2VV' ] ]
    time       = obsSetP2VV[0]
    angles     = obsSetP2VV[ 1 : 4 ]
    iTag       = obsSetP2VV[4]
    tagCatP2VV = obsSetP2VV[5]

    # tagging parameters
    numTagCats    = pdfBuild['tagCats']['numTagCats']
    tagCat5Min    = pdfBuild['tagCats'].traditionalCatRange(5)[0]
    taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
    tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

    # tagging category ranges
    tagCatP2VV.setRange( 'UntaggedRange', 'Untagged'    )
    tagCatP2VV.setRange( 'TaggedRange',   taggedCatsStr )
    tagCatP2VV.setRange( 'TagCat5Range',  tagCat5Str    )

    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

    if pdfConfig['SFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
                   , ( '', angleNames[0][1], angleNames[1][1], angleNames[1][1] )
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , ( True, ) + 3 * ( False, )
                  ) :
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth       )
             , components = comps
            )

    # plot lifetime
    timePlotTitles = tuple( [ time.GetTitle() + title for title in (  ' - linear'
                                                                    , ' - logarithmic'
                                                                    , ' - B/#bar{B} asymmetry'
                                                                   )
                            ] )
    timeCanv = TCanvas( 'timeCanv', 'Lifetime' )
    for ( pad, nBins, plotTitle, yTitle, yScale, dataCuts, pdfCuts, logY )\
            in zip(  timeCanv.pads( 2, 2 )
                   , numTimeBins
                   , timePlotTitles
                   , 2 * ( None, ) + ( 'B/#bar{B} asymmetry', )
                   , ( ( None, None ), ( 50., None ), ( None, None ) )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , ( False, True, False )
                  ) :
        plot(  pad, time, fitData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts        )
             , components = comps
            )

    # plot lifetime (tagged/untagged)
    timePlotTitles1 = tuple( [ time.GetTitle() + title for title in (  ' - untagged'
                                                                     , ' - tagging category 2'
                                                                     , ' - tagging category %d' % tagCat5Min
                                                                     , ' - B/#bar{B} asymmetry untagged'
                                                                     , ' - B/#bar{B} asymmetry tagging category 2'
                                                                     , ' - B/#bar{B} asymmetry tagging category %d' % tagCat5Min
                                                                    )
                            ] )
    timeCanv1 = TCanvas( 'timeCanv1', 'Lifetime' )
    for ( pad, nBins, plotTitle, yTitle, dataCuts, pdfCuts, logY )\
        in zip(  timeCanv1.pads( 3, 2 )
             , 2 * numTimeBins
             , timePlotTitles1
             , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
             ,   ( dict( Cut = '%s == 0'  % ( tagCatP2VV.GetName()             )                   ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VV.GetName()             )                   ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VV.GetName(), tagCat5Min )                   ), )
               + ( dict( Cut = '%s == 0'  % ( tagCatP2VV.GetName()             ), Asymmetry = iTag ), )
               + ( dict( Cut = '%s == 2'  % ( tagCatP2VV.GetName()             ), Asymmetry = iTag ), )
               + ( dict( Cut = '%s == %d' % ( tagCatP2VV.GetName(), tagCat5Min ), Asymmetry = iTag ), )
             ,   ( dict( Slice = ( tagCatP2VV, 'Untagged'              )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat2'               )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat%d' % tagCat5Min )                   ), )
               + ( dict( Slice = ( tagCatP2VV, 'Untagged'              ), Asymmetry = iTag ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat2'               ), Asymmetry = iTag ), )
               + ( dict( Slice = ( tagCatP2VV, 'TagCat%d' % tagCat5Min ), Asymmetry = iTag ), )
             , 3 * ( False, ) + 3 * ( False, )
            ) :
        plot(  pad, time, fitData, pdf, yTitle = yTitle, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'            )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts        )
             , components = comps
            )

    # plot angles
    anglePlotTitles =   tuple(  [ angle.GetTitle()                            for angle in angles ]\
                              + [ angle.GetTitle() + ' - B/#bar{B} asymmetry' for angle in angles ] )
    anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yTitle, dataCuts, pdfCuts )\
            in zip(  anglesCanv.pads( 3, 2 )
                   , 2 * angles
                   , 2 * numAngleBins
                   , anglePlotTitles
                   , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[1][1] )
                   , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
                   , 3 * ( dict(), ) + 3 * ( dict( Asymmetry = iTag ), )
                   , 3 * ( dict(), ) + 3 * ( dict( Asymmetry = iTag ), )
                  ) :
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yTitle = yTitle
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                             )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = lineWidth, **pdfCuts         )
             , components = comps
            )

    # print canvas to file
    timeAnglesCanv.Print(plotsFile + '(')
    timeCanv.Print(plotsFile)
    timeCanv1.Print(plotsFile)
    if pdfConfig['SFit'] :
      anglesCanv.Print(plotsFile + ')')
    else :
      anglesCanv.Print(plotsFile)
      pdfBuild['bkgAnglesCanv'].Print(plotsFile + ')')

