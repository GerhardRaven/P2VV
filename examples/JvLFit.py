###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
readData                = True
generateData            = False
doFit                   = True
fastFit                 = False
makeDLLPlots            = True
makeObservablePlots     = False
plotAnglesNoEff         = False
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['blind']      = False
pdfConfig['nominalPdf'] = False
sumW2Error              = False

pdfConfig['numEvents'] = 30000

plotsFile = 'plots/JvLSFitBank.ps' if pdfConfig['SFit'] else 'JvLCFit.ps'

if readData :
    pdfConfig['nTupleName'] = 'DecayTree'
    pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
else :
    pdfConfig['nTupleName'] = None
    pdfConfig['nTupleFile'] = None

if generateData :
    dataSetName = 'JpsiphiData'
    dataSetFile = 'JvLFit.root'

# fit options
fitOpts = dict(  NumCPU              = 1
               , Optimize            = 1
#               , Timer               = 1
#               , Minos               = False
#               , Hesse               = False
#               , Save                = True
#               , Minimizer           = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts

# plot options
lineWidth = 2
markStyle = 8
markSize  = 0.4

# PDF options
pdfConfig['transversityAngles'] = False
pdfConfig['bkgAnglePdf']        = ''
pdfConfig['sigTaggingPdf']      = 'tagCats'
pdfConfig['bkgTaggingPdf']      = 'tagCatsRelative'
pdfConfig['multiplyByTimeEff']  = ''

pdfConfig['conditionalTagging'] = False
pdfConfig['continuousEstWTag']  = False
pdfConfig['numEstWTagBins']     = 100

pdfConfig['eventTimeResolution'] = True
pdfConfig['numTimeResBins']      = 100

nEvents = 30000
pdfConfig['signalFraction'] = 0.67
pdfConfig['massRangeBackground'] = False

pdfConfig['amplitudeParam'] = 'bank' # 'phasesSWaveFrac'
pdfConfig['polarSWave']     = True
pdfConfig['AparParam']      = 'phase'

pdfConfig['carthLambdaCP'] = False

if not readData :
    pdfConfig['tagCats'] = [  ( 'Untagged',  0, 0.500001, 0.500, 0.505, 0., 0.648, 0. )
                            , ( 'TagCat1',   1, 0.499999, 0.484, 0.488, 0., 0.014, 0. )
                            , ( 'TagCat2',   2, 0.478,    0.467, 0.471, 0., 0.056, 0. )
                            , ( 'TagCat3',   3, 0.457,    0.447, 0.450, 0., 0.054, 0. )
                            , ( 'TagCat4',   4, 0.437,    0.427, 0.430, 0., 0.046, 0. )
                            , ( 'TagCat5',   5, 0.417,    0.408, 0.410, 0., 0.037, 0. )
                            , ( 'TagCat6',   6, 0.399,    0.390, 0.391, 0., 0.031, 0. )
                            , ( 'TagCat7',   7, 0.381,    0.372, 0.373, 0., 0.023, 0. )
                            , ( 'TagCat8',   8, 0.363,    0.354, 0.354, 0., 0.020, 0. )
                            , ( 'TagCat9',   9, 0.344,    0.334, 0.333, 0., 0.014, 0. )
                            , ( 'TagCat10', 10, 0.324,    0.313, 0.311, 0., 0.012, 0. )
                            , ( 'TagCat11', 11, 0.303,    0.291, 0.289, 0., 0.013, 0. )
                            , ( 'TagCat12', 12, 0.280,    0.269, 0.266, 0., 0.010, 0. )
                            , ( 'TagCat13', 13, 0.257,    0.246, 0.241, 0., 0.009, 0. )
                            , ( 'TagCat14', 14, 0.233,    0.222, 0.217, 0., 0.005, 0. )
                            , ( 'TagCat15', 15, 0.208,    0.194, 0.189, 0., 0.003, 0. )
                            , ( 'TagCat16', 16, 0.181,    0.168, 0.161, 0., 0.003, 0. )
                            , ( 'TagCat17', 17, 0.153,    0.141, 0.134, 0., 0.002, 0. )
                            , ( 'TagCat18', 18, 0.124,    0.113, 0.104, 0., 0.000, 0. )
                           ]

pdfConfig['timeEffHistFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root'
pdfConfig['timeEffHistName'] = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins'
#pdfConfig['timeEffHistName'] = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_20bins'

pdfConfig['angEffMomentsFile'] = 'effMomentsTransBasis' if pdfConfig['nominalPdf'] or pdfConfig['transversityAngles']\
                                 else 'effMomentsHelBasis'
#pdfConfig['angEffMomentsFile'] = 'effmoments_tcut_0.3_Feb.txt'
#pdfConfig['angEffMomentsFile'] = None

if pdfConfig['nominalPdf'] or pdfConfig['transversityAngles'] :
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

numBins = ( 60, 30, 30, 30 )
pdfConfig['numTimeBins'] = 30
numAngleBins = ( 20, 20, 20 )
pdfConfig['numAngleBins'] = ( 5, 7, 9 )


###########################################################################################################################################
## build PDF ##
###############

from P2VVLoad import RooFitOutput

# workspace
from RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

# get variables
obsSetP2VV = [ pdfBuild['observables'][obs] for obs in [ 'time', 'cpsi', 'ctheta', 'phi', 'iTag' ] ]
time       = obsSetP2VV[0]
angles     = obsSetP2VV[ 1 : 4 ]
iTag       = obsSetP2VV[4]
BMass      = pdfBuild['observables']['BMass']
estWTag    = pdfBuild['observables']['estWTag']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfBuild['iTagZeroTrick'] :
    tagCatP2VV = pdfBuild['observables']['tagCatP2VV']
    obsSetP2VV.append(tagCatP2VV)

    # tagging parameters
    numTagCats    = pdfBuild['tagCats']['numTagCats']
    tagCat5Min    = pdfBuild['tagCats'].traditionalCatRange(5)[0]
    taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
    tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

    # tagging category ranges
    tagCatP2VV.setRange( 'UntaggedRange', 'Untagged'    )
    tagCatP2VV.setRange( 'TaggedRange',   taggedCatsStr )
    tagCatP2VV.setRange( 'TagCat5Range',  tagCat5Str    )

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )


###########################################################################################################################################
## generate data ##
###################

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


###########################################################################################################################################
## fit data ##
##############

if doFit :
    # float/fix values of some parameters
    if pdfConfig['nominalPdf'] : pdfBuild['lambdaCP'].setConstant('lambdaCPSq')
    for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] :
        CEvenOdd.setConstant('avgCEven.*')
        if pdfConfig['nominalPdf'] : CEvenOdd.setConstant( 'avgCOdd.*', True )

    if pdfConfig['nominalPdf'] :
        pdfBuild['taggingParams'].setConstant( 'tagCatCoef.*', False )

    if fastFit :
        pdfBuild['lambdaCP'].setConstant('lambdaCPSq')
        for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
        pdfBuild['tagCats'].setConstant('.*')
        pdfBuild['lifetimeParams'].setConstant('dM|Gamma')
        pdfBuild['timeResModel'].setConstant('.*')
        pdfBuild['signalBMass'].setConstant('.*')
        if not pdfConfig['SFit'] :
            pdfBuild['backgroundBMass'].setConstant('.*')
            pdfBuild['backgroundTime'].setConstant('.*')

    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if pdfConfig['SFit'] : fitResult = pdf.fitTo( fitData, SumW2Error = sumW2Error, **fitOpts )
    else                 : fitResult = pdf.fitTo( fitData,                          **fitOpts )

    print 120 * '=' + '\n'


###########################################################################################################################################
## make some plots ##
#####################

if makeObservablePlots or pdfConfig['makePlots'] or makeDLLPlots :
    # import plotting tools
    from P2VVLoad import ROOTStyle
    from P2VVGeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

    # create projection data set for conditional observables
    if pdfConfig['SFit'] :
        comps = None
    else :
        comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                }

    projWDataSet = []
    if   pdfConfig['continuousEstWTag']   : projWDataSet += [ tagCatP2VV, estWTag, iTag ]
    elif pdfConfig['conditionalTagging']  : projWDataSet += [ tagCatP2VV, iTag ]
    if   pdfConfig['eventTimeResolution'] : projWDataSet += [ timeRes ]

    if projWDataSet :
        bulkData = fitData.reduce( CutRange = 'Bulk' )
        projWData     = dict( ProjWData = ( fitData.reduce(  ArgSet = projWDataSet ), True ) )
        projWDataBulk = dict( ProjWData = ( bulkData.reduce( ArgSet = projWDataSet ), True ) )
    else :
        projWData     = dict()
        projWDataBulk = dict()

if pdfConfig['makePlots'] :
    # plot background time
    bkgTimeCanv = TCanvas( 'bkgTimeCanv', 'Background Lifetime' )
    for ( pad, data, plotTitle, logY )\
          in zip(  bkgTimeCanv.pads( 2, 2 )
                 , 2 * [ pdfBuild['bkgRangeData'], pdfBuild['bkgSWeightData'] ]
                 , [  time.GetTitle() + ' - mass side bands - linear'
                    , time.GetTitle() + ' - mass S-weights - linear'
                    , time.GetTitle() + ' - mass side bands - logarithmic'
                    , time.GetTitle() + ' - mass S-weights - logarithmic'
                   ]
                 , 2 * [ False ] + 2 * [ True ]
                ) :
        plot(  pad, time, data, pdfBuild['bkg_t'], logy = logY
             , frameOpts  = dict( Title = plotTitle, Range = 'Bulk', Bins = pdfConfig['numTimeBins'] )
             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
             , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
            )

if makeObservablePlots and not pdfBuild['iTagZeroTrick'] :
    # plot lifetime and angles
    timeAnglesCanv = TCanvas( 'timeAnglesCanv', 'Lifetime and Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
            in zip(  timeAnglesCanv.pads( 2, 2 )
                   , obsSetP2VV[ : 5 ]
                   , numBins
                   , [ var.GetTitle() for var in obsSetP2VV[ : 5 ] ]
                   , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                   , ( True, ) + 3 * ( False, )
                  ) :
        plot(  pad, obs, fitData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle                                     )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize                      )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth )
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
                   , 3 * [ pdfConfig['numTimeBins'] ]
                   , timePlotTitles
                   , 2 * ( None, ) + ( 'B/#bar{B} asymmetry', )
                   , ( ( None, None ), ( 50., None ), ( None, None ) )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , 2 * ( dict(), ) + ( dict( Asymmetry = iTag ), )
                   , ( False, True, False )
                  ) :
        plot(  pad, time, fitData, pdf, yTitle = yTitle, yScale = yScale, logy = logY
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                    )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                         )
             , pdfOpts    = dict( list( projWDataBulk.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
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
             , 6 * [ pdfConfig['numTimeBins'] ]
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
             , frameOpts  = dict( Bins = nBins, Title = plotTitle, Range = 'Bulk'                                )
             , dataOpts   = dict( MarkerStyle = markStyle, MarkerSize = markSize, **dataCuts                     )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , components = comps
            )

    # plot angles
    if plotAnglesNoEff and pdfConfig['SFit'] and pdfConfig['multiplyByTimeEff'] not in [ 'all', 'signal' ]\
            and ( pdfConfig['nominalPdf'] or not pdfConfig['conditionalTagging'] ) :
        addPDFs = [ ws['sig_t_angles_tagCat_iTag'] ]
    else :
        addPDFs = [ ]

    anglePlotTitles =   tuple(  [ angle.GetTitle()                            for angle in angles ]\
                              + [ angle.GetTitle() + ' - B/#bar{B} asymmetry' for angle in angles ] )
    anglesCanv = TCanvas( 'anglesCanv', 'Decay Angles' )
    for ( pad, obs, nBins, plotTitle, xTitle, yTitle, dataCuts, pdfCuts )\
            in zip(  anglesCanv.pads( 3, 2 )
                   , 2 * angles
                   , 2 * numAngleBins
                   , anglePlotTitles
                   , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                   , 3 * ( None, ) + 3 * ( 'B/#bar{B} asymmetry', )
                   , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTag ), )
                   , 3 * ( dict( ), ) + 3 * ( dict( Asymmetry = iTag ), )
                  ) :
        plot(  pad, obs, fitData, pdf, addPDFs = addPDFs, xTitle = xTitle, yTitle = yTitle
             , frameOpts   = dict( Bins = nBins, Title = plotTitle                                                )
             , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize , **dataCuts                    )
             , pdfOpts     = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = lineWidth, **pdfCuts )
             , addPDFsOpts = [ dict( list( projWData.items() ), LineColor = kRed, LineWidth = lineWidth, **pdfCuts ) ]
             , components  = comps
            )

    if not pdfConfig['SFit'] :
        # plot signal mass
        pad = pdfBuild['massCanv'].cd(2)
        plot(  pad, BMass, fitData, pdf
             , frameOpts  = dict( Range = 'Signal', Bins = pdfConfig['numBMassBins'][0], Title = BMass.GetTitle() + ' full fit - signal' )
             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                                      )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 2                                            )
             , components = comps
            )

    # print canvas to file
    timeAnglesCanv.Print(plotsFile + '(')
    timeCanv.Print(plotsFile)
    timeCanv1.Print(plotsFile)
    if pdfConfig['makePlots'] :
        anglesCanv.Print(plotsFile)
        pdfBuild['massCanv'].Print(plotsFile)
        bkgTimeCanv.Print(plotsFile)
        pdfBuild['bkgAnglesCanv'].Print(plotsFile)
        pdfBuild['estWTagCanv'].Print(plotsFile + ')')

    else :
        anglesCanv.Print(plotsFile + ')')

elif pdfConfig['makePlots'] :
    pdfBuild['massCanv'].Print(plotsFile + '(')
    bkgTimeCanv.Print(plotsFile)
    pdfBuild['bkgAnglesCanv'].Print(plotsFile)
    pdfBuild['estWTagCanv'].Print(plotsFile + ')')


if makeDLLPlots :
    wsPars =\
        {  'ReApar'       : ( '#DeltaNLL Re(A_{#parallel})',           'Re(A_{#parallel})',           -0.480, -0.455, 10, 0.001, 0.01 )
         , 'ImApar'       : ( '#DeltaNLL Im(A_{#parallel})',           'Im(A_{#parallel})',           -0.2,    0.2,   10, 0.001, 0.01 )
         , 'cosAparPhase' : ( '#DeltaNLL cos(#delta_{#parallel})',     'cos(#delta_{#parallel})',     -1.,    -0.94,  10, 0.001, 0.01 )
         , 'AparPhase'    : ( '#DeltaNLL #delta_{#parallel}',          '#delta_{#parallel}',           2.8,    3.5,   10, 0.001, 0.01 )
         , 'sqrtfS_Re'    : ( '#DeltaNLL #sqrt{f_{S}}^{R}',            '#sqrt{f_{S}}^{R}',            -0.22,  -0.10,  10, 0.001, 0.01 )
         , 'sqrtfS_Im'    : ( '#DeltaNLL #sqrt{f_{S}}^{I}',            '#sqrt{f_{S}}^{I}',            -0.060,  0.085, 10, 0.001, 0.01 )
         , 'ReASOdd'      : ( '#DeltaNLL Re(A_{S} / A_{#perp})',       'Re(A_{S} / A_{#perp})',        0.20,   0.44,  10, 0.001, 0.01 )
         , 'ImASOdd'      : ( '#DeltaNLL Im(A_{S} / A_{#perp})',       'Im(A_{S} / A_{#perp})',       -0.06,   0.04,  10, 0.001, 0.01 )
         , 'f_S'          : ( '#DeltaNLL f_{S}',                       'f_{S}',                        0.00,   0.05,  10, 0.001, 0.01 )
         , 'ASPhase'      : ( '#DeltaNLL #delta_{S}',                  '#delta_{S}',                   2.75,   3.30,  10, 0.001, 0.01 )
         , 'ASOddMag2'    : ( '#DeltaNLL A_{S}^{2} / A_{#perp}^{2}',   'A_{S}^{2} / A_{#perp}^{2}',    0.0,    0.2,   10, 0.001, 0.01 )
         , 'ASOddPhase'   : ( '#DeltaNLL #delta_{S} - #delta_{#perp}', '#delta_{S} - #delta_{#perp}', -0.2,    0.2,   10, 0.001, 0.01 )
        }

    dllPars = [ ]
    if   pdfConfig['AparParam'] == 'real' : dllPars.append('ReApar')
    elif pdfConfig['AparParam'] == 'ReIm' : dllPars.append('ImApar')
    elif pdfConfig['AparParam'] == 'cos'  : dllPars.append('cosAparPhase')
    else                                  : dllPars.append('AparPhase')

    if pdfConfig['polarSWave'] :
        if pdfConfig['amplitudeParam'] == 'bank' : dllPars += [ 'ASOddMag2', 'ASOddPhase' ]
        else                                     : dllPars += [ 'f_S', 'ASPhase' ]
    else :
        if pdfConfig['amplitudeParam'] == 'bank' : dllPars += [ 'ReASOdd', 'ImASOdd' ]
        else                                     : dllPars += [ 'sqrtfS_Im', 'sqrtfS_Im' ]

    # float/fix values of some parameters
    pdfBuild['lambdaCP'].setConstant('lambdaCPSq')
    for CEvenOdd in pdfBuild['taggingParams']['CEvenOdds'] : CEvenOdd.setConstant('avgCEven.*|avgCOdd.*')
    pdfBuild['tagCats'].setConstant('.*')
    pdfBuild['lifetimeParams'].setConstant('dM|Gamma')
    pdfBuild['timeResModel'].setConstant('.*')
    pdfBuild['signalBMass'].setConstant('.*')
    if not pdfConfig['SFit'] :
        pdfBuild['backgroundBMass'].setConstant('.*')
        pdfBuild['backgroundTime'].setConstant('.*')

    # create DNLL/PLL plots
    from ROOT import RooFit, RooArgSet, TCanvas
    nll = pdf.createNLL( fitData, **fitOpts )
    dllCanvs = [ ]
    canvFileName = plotsFile[ : -3 ] + 'DLLs.ps'
    for parIter, par in enumerate(dllPars) :
        pll = nll.createProfile( RooArgSet( ws[par] ) )
        parFrame = ws[par].frame(  RooFit.Range( wsPars[par][2], wsPars[par][3] )
                                 , RooFit.Bins( wsPars[par][4] )
                                 , RooFit.Title( wsPars[par][0] )
                                )

        print 'JvLFit: plotting Delta -log(L) for %s' % par
        nll.plotOn( parFrame, RooFit.ShiftToZero(), RooFit.LineColor(kBlue), RooFit.Precision( wsPars[par][5] ) )

        print 'JvLFit: plotting profile Delta -log(L) for %s' % par
        pll.plotOn( parFrame, RooFit.LineColor(kRed), RooFit.Precision( wsPars[par][6] ) )

        parFrame.GetXaxis().SetTitle( wsPars[par][1] )
        parFrame.GetYaxis().SetTitle('#DeltaNLL')

        dllCanvs.append( TCanvas( 'dllCanv%d' % parIter , 'DLL canvas' ) )
        parFrame.Draw()

        for canvIter, canv in enumerate(dllCanvs) :
            if len(dllCanvs) == 1 or canvIter not in [ 0, len(dllCanvs) - 1 ] : namePF = ''
            elif canvIter == 0 : namePF = '('
            else : namePF = ')'
            canv.Print( canvFileName + namePF )
