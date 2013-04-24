###########################################################################################################################################
## set script parameters ##
###########################

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
doFit                   = True
randomParVals           = ( ) # ( 1., 12346 ) # ( 2., 12345 )

#OutputPath for the plots in line 

pdfConfig['nTupleName'] = 'DecayTree'
pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'

# fit options
fitOpts = dict(  NumCPU    = 8
               , Optimize  = 2
               , Timer     = True
               , Minimizer = 'Minuit2'
              )
pdfConfig['fitOptions'] = fitOpts

# PDF options
pdfConfig['transversityAngles'] = False
pdfConfig['bkgAnglePdf']          = 'hybrid'            # 'hybrid'
pdfConfig['sigTaggingPdf']        = 'tagUntag'          # 'tagUntag' # nominal: 'tagCats'
pdfConfig['bkgTaggingPdf']        = 'tagUntagRelative'  # 'tagUntagRelative' # 'tagCatsRelative'
pdfConfig['multiplyByTagPdf']     = False
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012'         # 'paper2012' # 'HLT1Unbiased'
pdfConfig['multiplyByAngEff']     = 'basis012'          # 'basis012'
pdfConfig['parameterizeKKMass']   = 'simultaneous'      # 'simultaneous'
pdfConfig['ambiguityParameters']  = False
pdfConfig['lifetimeRange']        = ( 0.3, 14. )
pdfConfig['SWeightsType']         = 'simultaneousFreeBkg'  # 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ] # [ 1008., 1032. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.33, 0.09), (0.073, 0.030), (0.009, 0.012), (0.012, 0.010), (0.061, 0.027), (0.18, 0.04) ]
                                    , [ (1.1,  0.5 ), (0.7,   0.2  ), (0.4,   0.4  ), (-0.6,  0.3  ), (-0.4, 0.2   ), (-0.7, 0.2 ) ] )
pdfConfig['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ] # [ 0.498 ] # [ 0.326 ] # [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

pdfConfig['sameSideTagging']    = True
pdfConfig['conditionalTagging'] = True
pdfConfig['continuousEstWTag']  = True
pdfConfig['numEstWTagBins']     = 20
pdfConfig['constrainTagging']   = 'constrain'  # 'constrain'

pdfConfig['timeResType']           = 'eventNoMean' # 'event' # 'eventNoMean'
pdfConfig['numTimeResBins']        = 50
pdfConfig['constrainTimeResScale'] = 'constrain'  # nominal: 'constrain'

pdfConfig['numEvents'] = 10000
pdfConfig['signalFraction'] = 0.45
pdfConfig['massRangeBackground'] = True

pdfConfig['amplitudeParam'] = 'phasesSWaveFrac' # 'bank' # 'phasesSWaveFrac'
pdfConfig['ASParam']        = 'deltaPerp'  # 'deltaPerp'
pdfConfig['AparParam']      = 'phase' # 'Mag2ReIm' # 'phase'

pdfConfig['constrainDeltaM'] = 'constrain'  # nominal: 'constrain'

pdfConfig['lambdaCPParam'] = 'lambPhi'  # 'lambPhi'

constWTagAsyms = 'P1'

pdfConfig['timeEffHistFile']      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/timeAcceptanceStartValues.root'\
                                    if pdfConfig['timeEffType'] == 'fit' else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomentsFile']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'\
                                    if and pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'

if pdfConfig['transversityAngles'] :
    pdfConfig['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                               , ( 'trcostheta', 'cos(#theta_{tr})' )
                               , ( 'trphi',      '#phi_{tr}'        )
                              )
else :
    pdfConfig['angleNames'] = (  ( 'helcosthetaK', 'cos #theta_{K}' )
                               , ( 'helcosthetaL', 'cos #theta_{#mu}' )
                               , ( 'helphi',       '#phi_{h} [rad]'            )
                              )

numBins = ( 60, 30, 30, 30 )
pdfConfig['numTimeBins'] = 30
numAngleBins = ( 20, 20, 20 )
pdfConfig['numAngleBins'] = ( 5, 7, 9 )


###########################################################################################################################################
## build PDF ##
###############

from P2VV.Load import RooFitOutput

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

# get variables
obsSetP2VV = [ pdfBuild['observables'][obs] for obs in [ 'time', 'cpsi', 'ctheta', 'phi', 'iTagOS' ] ]
time       = obsSetP2VV[0]
angles     = obsSetP2VV[ 1 : 4 ]
iTagOS     = obsSetP2VV[4]
iTagSS     = pdfBuild['observables']['iTagSS']
BMass      = pdfBuild['observables']['BMass']
mumuMass   = pdfBuild['observables']['mumuMass']
KKMass     = pdfBuild['observables']['KKMass']
estWTagOS  = pdfBuild['observables']['estWTagOS']
timeRes    = pdfBuild['observables']['timeRes']

if not pdfConfig['SFit'] : obsSetP2VV.append(BMass)

if not pdfBuild['iTagZeroTrick'] :
    tagCatP2VVOS = pdfBuild['observables']['tagCatP2VVOS']
    tagCatP2VVSS = pdfBuild['observables']['tagCatP2VVSS']
    obsSetP2VV.append(tagCatP2VVOS)

    # tagging parameters
    numTagCats    = pdfBuild['tagCatsOS']['numTagCats']
    tagCat5Min    = pdfBuild['tagCatsOS'].traditionalCatRange(5)[0]
    taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
    tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

    # tagging category ranges
    tagCatP2VVOS.setRange( 'UntaggedRange', 'Untagged'    )
    tagCatP2VVOS.setRange( 'TaggedRange',   taggedCatsStr )
    tagCatP2VVOS.setRange( 'TagCat5Range',  tagCat5Str    )


###########################################################################################################################################
## get data ##
##############

if pdfConfig['SFit'] :
    defData = pdfBuild['sigSWeightData']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']
    from P2VV.GeneralUtils import correctSWeights
    fitData = correctSWeights( pdfBuild['sigSWeightData'], 'N_bkgMass_sw'
                              , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else '' )

else :
    defData = pdfBuild['data']
    sigData = pdfBuild['sigSWeightData']
    bkgData = pdfBuild['bkgSWeightData']
    fitData = pdfBuild['data']

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(defData)
pdfPars = pdf.getParameters(defData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
from math import sqrt
for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
    if not pdfConfig['sameSideTagging'] :
        CEvenOdds.setConstant('avgCEven.*')
        CEvenOdds.setConstant( 'avgCOdd.*', True )
    else :
        for CEvenOdd in CEvenOdds :
            CEvenOdd.setConstant('avgCEven.*')
            CEvenOdd.setConstant( 'avgCOdd.*', True )

pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

pdfBuild['amplitudes'].setConstant('C_SP')

if randomParVals :
    import random
    # give parameters random offsets
    print 'JvLFit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'JvLFit: fit data:'
fitData.Print()
print 'JvLFit: observables in PDF:'
pdfObs.Print('v')
print 'JvLFit: parameters in PDF:'
pdfPars.Print('v')

if doFit :
    # fit data
    print 120 * '='
    print 'JvLFit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    if pdfConfig['SFit'] : fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, **fitOpts )
    else                 : fitResult = pdf.fitTo( fitData,                     Save = True, **fitOpts )

    fitResult.Print()


###########################################################################################################################################
## make plots ##
#####################

# import plotting tools
from P2VV.Load import LHCbStyle
from P2VV.GeneralUtils import plot, CPcomponentsPlotingToolkit
from ROOT import TCanvas, kRed, kGreen, kMagenta, kBlue, kSolid

#Initialaze the CP components ploting toolkit
CpPlotsKit = CPcomponentsPlotingToolkit(pdf,defData)

#Get dictionary with all the pdfs in the KKmass bins
# {'bin_i', {'total'=... , 'even'=... , 'odd'=... , 'swave'=...}}
pdfsDict = CpPlotsKit.getCPcompPdfKKbins()

#Get some useful stuff ncessesary for looping
KKbins      = CpPlotsKit.getNumKKbins()  #Get nummber of KKmass bins 
binNames    = CpPlotsKit.getKKbinNames() #Get list of KKmass bin names
CPcomps     = CpPlotsKit.getCpCompNames()#Get list of names of the CP components
angleNames  = pdfConfig['angleNames']
observables = [time] + angles

#Set plot options      
markStyle = 8
markSize  = 0.5
lineWidth = 4
CpPlotsKit.setLineColors( dict(total = kBlue, even=kRed, odd=kGreen+3, swave=kMagenta+3) )
CpPlotsKit.setLineStyles( dict(total = kSolid, even=9   , odd=7       , swave=5         ) )
CpPlotsKit.setLineWidth(lineWidth)


#LHCbLabel
from ROOT import TPaveText, gStyle
lhcbName = TPaveText(gStyle.GetPadLeftMargin() + 0.14,
                         0.87 - gStyle.GetPadTopMargin(),
                         gStyle.GetPadLeftMargin() + 0.20,
                         0.95 - gStyle.GetPadTopMargin(),
                         "BRNDC")
lhcbName.AddText("LHCb")
lhcbName.SetFillColor(0)
lhcbName.SetTextAlign(12)
lhcbName.SetBorderSize(0)


#Get pdf and AdPdf options in case of bining and No bining the projData 
PDFoptsBin   = CpPlotsKit.getPdfOptsSixKKbins(BinData=True)
PDFoptsNOBin =  CpPlotsKit.getPdfOptsSixKKbins(BinData=False)
addPDFotpsBin = CpPlotsKit.getAddPdfsOptsSixKKbins(BinData=True)
addPDFotpsNOBin = CpPlotsKit.getAddPdfsOptsSixKKbins(BinData=False)

#Plot and save  lifetime and angles in all the KK mass bins
for bin in binNames:
         print '\n\nP2VV - INFO: Plotting decay time and angular distributions of ' + bin
         dataSlice = defData.reduce('KKMassCat==KKMassCat::' + bin)
         pdfSlice  = pdfsDict[bin]['total']
         binIdx = binNames.index(bin)
         for ( pad, obs, nBins, xTitle, yScaleRel, logY )\
                 in zip(  [ TCanvas(o.GetName()+bin) for o in observables ]
                        , observables
                        , numBins
                        , ( time.GetTitle()+' [ps]', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                        , 2 *(( 1., 1.2 ),) + 2 *(( 1. , 1.2 ),) 
                        , ( True, ) + 3 * ( False, )
                       ) :           
             plot(  pad, obs, dataSlice, pdfSlice, xTitle=xTitle, yScaleRel=yScaleRel, logy=logY
                  , frameOpts   = dict( Bins = nBins, Name = bin + obs.GetName() + 'Histo'      )
                  , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize          )
                  , pdfOpts     = PDFoptsNOBin[bin] if obs==time else PDFoptsBin[bin]
                  , addPDFs     = [ pdfsDict[bin][c] for c in CPcomps ]
                  , addPDFsOpts = [addPDFotpsNOBin[binIdx][c] for c in CPcomps ] if obs==time\
                            else  [addPDFotpsBin[binIdx][c] for c in CPcomps ]
                   )
             lhcbName.Draw()
        ## print canvas to file
             fName =  bin + '_' + obs.GetName() + '_sFit.ps' if pdfConfig['SFit'] \
                else  bin + '_' + obs.GetName() + '_cFit.ps'
             pad.Print(fName)
                      
# Save all the plots in a root file as RooPlot objects.
from P2VV.GeneralUtils import _P2VVPlotStash as rooplots
from ROOT import TFile
plotsFile = TFile('RooPlots6x4.root','recreate')
for plot in rooplots: plot.Write()
plotsFile.Close()

