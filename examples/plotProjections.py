###########################################################################################################################################
## set script parameters ##
###########################

from P2VVParameterizations.FullPDFs import Bs2Jpsiphi_Winter2012 as PdfConfig
pdfConfig = PdfConfig()

# job parameters
pdfConfig['selection']  = 'paper2012' # 'paper2012' # 'HLT1Unbiased'
pdfConfig['makePlots']  = False
pdfConfig['SFit']       = True
pdfConfig['nominalPdf'] = False  # nominal PDF option does not work at the moment
doFit                   = False
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

# plot options
lineWidth = 2
markStyle = 8
markSize  = 0.4

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
                                    if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] else\
                                    '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis'

if not pdfConfig['nominalPdf'] and pdfConfig['transversityAngles'] :
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
    from P2VVGeneralUtils import correctSWeights
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
## make some plots ##
#####################

# import plotting tools
from P2VVLoad import LHCbStyle
from P2VVGeneralUtils import plot, CPcomponentsPlotingToolkit
from ROOT import TCanvas, kBlack, kBlue, kRed, kGreen, kDashed, kFullCircle, kFullSquare
from RooFitWrappers import SimultaneousPdf

#Initialaze the CP components ploting toolkit
CpPlotsKit = CPcomponentsPlotingToolkit(pdf,defData)

KKbins = CpPlotsKit.getNumKKbins()
binNames = CpPlotsKit.getKKbinNames()
CPcomps = CpPlotsKit.getCpCompNames()

#Choose line style and Color:
lineStyles = dict(even=9, odd=3, swave=5)
lineColors = dict(even='',odd='',swave='kRed')
CpPlotsKit.setLineColors(lineColors)
CpPlotsKit.setLineStyles(lineStyles)

#Construct a SuperDictionary of the pdf of the f0orm:
  # {'bin_i' : {'total':pdfTotal_i, 'even':pdfEven_i, 'odd':pdfOdd_i, 'swave:pdfSwave_i'}   }
#pdfsDict = CpPlotsKit.getCPcompPdfKKbins()
  # also include an entry for the total complete pdf: {'complete', .}
#pdfs = CpPlotsKit.getCPcompPdf()

#Get the relative normalitaion fractions of pdf components in a SuperDictionary
#CPfracs = CpPlotsKit.getCPnormFracs()

#Get the relative normalitaion fractions of Simultaneous Pdf components
#SliceFracs = CpPlotsKit.getKKslicesNormFracs()

#Get pdf dwaing options (ProjectCondObs and LineStyle)
CpPlotsKit.setLineWidth(lineWidth)
PDFopts = CpPlotsKit.getPdfOpts()

#Get the additional pdfs ( Simulpdf.getPdf(KKbin) )
ADDpdfs =  CpPlotsKit.getAddPdfs()

#Get pdf dwaing options of the additional pdfs
# (ProjectCondObs, Norms, AddTo, Invisible and LineStyle)
ADDpdfsOpts = CpPlotsKit.getAddPdfOpts()


#Report: 2 assumptions, Negative normalization integral.



        
assert(False)

projDSet = defData.reduce('KKMassCat==KKMassCat::bin0')
PDFopts.pop('ProjWData')
PDFopts.update(dict(ProjWData = (projDSet,False)))
add = [ADDpdfs[0],ADDpdfs[1]]
addOpts = [ADDpdfsOpts[0],ADDpdfsOpts[1]]

addOpts[0].pop('ProjWData')
addOpts[1].pop('ProjWData')
addOpts[1].pop('Invisible')

c =TCanvas()
#Plot and Save
#Make a single canvas for all KK mass bins (Using the add six curves feature)
timeAnglesCanv = TCanvas('timeAnglesCanv')
OutputFilename = 'plotProjections_sFit.ps' if pdfConfig['SFit'] else 'plotProjections_cFit.ps'
         
#Plot background if cFit.
if pdfConfig['SFit'] : comps = None
else :
    comps = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
             , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )               
            }


#plot aobservables
for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
                 in zip(  timeAnglesCanv.pads( 2, 2 )
                        , obsSetP2VV[ : 4]
                        , numBins
                        , [ var.GetTitle() for var in obsSetP2VV[ : 4 ] ]
                        , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                        , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                        , ( True, ) + 3 * ( False, )
                       ) :
             plot(  pad, obs, defData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
                  , frameOpts   = dict( Bins = nBins, Title = plotTitle                    )
                  , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize     )
                  , pdfOpts     = PDFopts
                  , addPDFs     = [(pdfsDict['complete']['even'].getPdf(bin))for bin in binNames] +\
                                  [(pdfsDict['complete']['odd'].getPdf(bin))for bin in binNames]  +\
                                  [(pdfsDict['complete']['swave'].getPdf(bin))for bin in binNames]
                  , addPDFsOpts = ADDpdfsOpts
                  , components  = comps
                 )
             

             timeAnglesCanv.Print(OutputFilename)






# Start from here tomorroww

plot(  c, angles[0], defData, pdf = bin0, pdfOpts = dict(LineColor = 2, Invisible = ()))











assert(False)




plot(c3, angles[0], defData, addPDFs = ADDpdfs, addPDFsOpts = ADDpdfsOpts )      
#Plot lifetime and angles in all the KK mass bins

#Draw options for individual KK mass bins ploting 
## CpCompDrawOptsDict = dict((bin, [dict(list(projWDataDict[bin].items()), Normalization = fracs[bin]['even'] , LineStyle = 9 ) , 
##                                  dict(list(projWDataDict[bin].items()), Normalization = fracs[bin]['odd']  , LineStyle = 3 ) ,
##                                  dict(list(projWDataDict[bin].items()), Normalization = fracs[bin]['swave'], LineStyle = 5, LineColor = 2 )]) for bin in binNames)


for bin in xrange(binNames):
         print '\n\nplotProjections: plotting lifetime and angular distributions of bin{0}'.format(bin)
         timeAnglesCanv = TCanvas( 'timeAnglesCanv_bin{0}'.format(bin), 'Lifetime and Decay Angles bin{0}'.format(bin) )

         #dataslice = defData.reduce('KKMassCat==KKMassCat::bin{0}'.format(bin))
         #pdf = pdfsDict['bin{0}'.format(bin)]['total']
         #OutputFilename = 'plotProjections_sFit_bin{0}.ps'.format(bin) if pdfConfig['SFit'] else 'plotProjections_cFit_bin{0}.ps'.format(bin)
         
         for ( pad, obs, nBins, plotTitle, xTitle, yScale, logY )\
                 in zip(  timeAnglesCanv.pads( 2, 2 )
                        , obsSetP2VV[ : 4 ]
                        , numBins
                        , [ var.GetTitle() for var in obsSetP2VV[ : 4 ] ]
                        , ( '', angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                        , ( ( 0.1, None ), ) + 3 * ( ( None, None ), )
                        , ( True, ) + 3 * ( False, )
                       ) :
             plot(  pad, obs, defData, pdf, xTitle = xTitle, yScale = yScale, logy = logY
                  , frameOpts   = dict( Bins = nBins, Title = plotTitle                    )
                  , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize     )
                 #, pdfOpts     = projWDataDict['bin{0}'.format(bin)]                          
                  , pdfOpts     = dict(LineWidth = lineWidth                               )
                 #, addPDFs     = [pdfsDict['bin{0}'.format(bin)]['even'] ,
                 #                 pdfsDict['bin{0}'.format(bin)]['odd']  ,
                 #                 pdfsDict['bin{0}'.format(bin)]['swave']] 
                 #, addPDFsOpts = CpCompDrawOptsDict['bin{0}'.format(bin)]
                  , components  = comps
                 )
             # print canvas to file
             timeAnglesCanv.Print(OutputFilename)


assert(False)










#How to plot slices of data.
slices = { 'Slices':[(trgCategory, 'exclB' ),(trgCategory , 'notExclB' )]  }






##Conditional Observables plots in KK mass bins
## from P2VVGeneralUtils import getCondObsPlotsInKKbins
## CondObsCanv = TCanvas('CondObsCanvInKKbins')
## condObsCanvs = getCondObsPlotsInKKbins(pdf, defData, CondObsCanv )
## assert(False)

