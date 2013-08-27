
###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
doFit                   = False
pdfConfig['selection']  = 'paper2012'

parFileIn  = ''
parFileOut = ''

pdfConfig['nTupleName'] = 'DecayTree'
pdfConfig['nTupleFile'] = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'


pdfConfig['timeEffHistFile']      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomentsFile']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/trans_UB_UT_trueTime_BkgCat050_KK30_Basis'




# fit options
fitOpts = dict(  NumCPU    = 2
               , Optimize  = 2
               , Minimizer = 'Minuit2'
#               , Offset    = True
#               , Hesse     = False
               , Timer     = True
#               , Verbose   = True
              )
pdfConfig['fitOptions'] = fitOpts

fitRange      = ''
corrSFitErr   = 'sumWeight' # '' / 'sumWeight' / ( 0.887, [ 0.566, 0.863, 0.956, 0.948, 0.855, 0.662 ] ) / 'matrix'
randomParVals = ( ) # ( 1., 12345 )
MinosPars     = [#  'AparPhase'
                 #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
                 #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
                ]

# PDF options
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012'
pdfConfig['multiplyByAngEff']     = 'basis012Plus'
pdfConfig['parameterizeKKMass']   = 'simultaneous'
pdfConfig['SWeightsType']         = 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.067, 0.029), (0.008, 0.011), (0.016, 0.011), (0.055, 0.026), (0.17,  0.04) ]
                                     , [ (1.3,  0.7 ), (0.77,  0.28 ), (0.50,  0.47 ), (-0.51, 0.25 ), (-0.46, 0.21 ), (-0.65, 0.20) ] )
pdfConfig['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

pdfConfig['sameSideTagging']    = True
pdfConfig['conditionalTagging'] = True
pdfConfig['continuousEstWTag']  = True
pdfConfig['constrainTagging']   = 'constrain'

pdfConfig['timeResType']           = 'eventNoMean'
pdfConfig['numTimeResBins']        = 40
pdfConfig['constrainTimeResScale'] = 'fixed'

pdfConfig['constrainDeltaM'] = 'constrain'

pdfConfig['lambdaCPParam'] = 'lambPhi'


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# signal and background data sets
sigData = pdfBuild['sigSWeightData']
bkgData = pdfBuild['bkgSWeightData']

# data set with weights corrected for background dilution: for phi_s fit only!
if corrSFitErr == 'sumWeight'\
        or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
    from P2VV.Utilities.DataHandling import correctSWeights
    fitData = correctSWeights( pdfBuild['sigSWeightData'], 'N_bkgMass_sw'
                              , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else ''
                              , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

else :
    fitData = pdfBuild['sigSWeightData']

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


###########################################################################################################################################
## fit data ##
##############

# float/fix values of some parameters
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
    # give parameters random offsets
    import random
    print 'Bs2JpsiKK2011Fit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
          % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
    random.seed( randomParVals[1] if randomParVals[1] else None )
    for par in pdfPars :
        if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() )

# print parameters
print 120 * '='
print 'Bs2JpsiKK2011Fit: fit data:'
fitData.Print()
print 'Bs2JpsiKK2011Fit: observables in PDF:'
pdfObs.Print('v')
print 'Bs2JpsiKK2011Fit: parameters in PDF:'
pdfPars.Print('v')

if doFit :
    # fit data
    print 120 * '='
    print 'Bs2JpsiKK2011Fit: fitting %d events (%s)' % ( fitData.numEntries(), 'weighted' if fitData.isWeighted() else 'not weighted' )

    RooMinPars = [ ]
    if MinosPars :
        print 'Bs2JpsiKK2011Fit: running Minos for parameters',
        for parName in MinosPars :
            RooMinPars.append( pdfPars.find(parName) )
            print '"%s"' % RooMinPars[-1],
        print

    fitResult = pdf.fitTo( fitData, SumW2Error = True if corrSFitErr == 'matrix' else False
                          , Minos = RooMinPars, Save = True, Range = fitRange
                          , **fitOpts
                         )

    # print parameter values
    from P2VV.Imports import parNames, parValues
    print 'Bs2JpsiKK2011Fit: parameters:'
    fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    fitResult.covarianceMatrix().Print()
    fitResult.correlationMatrix().Print()

    print 120 * '=' + '\n'

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )



###########################################################################################################################################
## make plots ##
#####################

# import plotting tools
from P2VV.Load import LHCbStyle
from P2VV.Utilities.Plotting import plot, CPcomponentsPlotingToolkit
from ROOT import TCanvas, kRed, kGreen, kMagenta, kBlue, kSolid

#LHCbLabel
from ROOT import TPaveText, gStyle
lhcbName = TPaveText(0.28, 0.77, 0.38, 0.90, "BRNDC")
lhcbName.AddText("LHCb")
lhcbName.SetFillColor(0)
lhcbName.SetTextAlign(12)
lhcbName.SetBorderSize(0)




assert(False)


data = fitData   # sigData fitData









#Initialaze the CP components ploting toolkit
CpPlotsKit = CPcomponentsPlotingToolkit(pdf,data)

#Get some useful stuff ncessesary for looping
angleNames  = pdfConfig['angleNames']
observables = list(pdf.Observables()-pdf.ConditionalObservables()) # DecayTime + angles
numBins = ( 30, 30, 30, 60 )

#Make the y-axis titles look nicer
yTitles = []
for obs in observables:
    range = obs.getRange()[1] - obs.getRange()[0]
    binWidth = round( range / numBins[observables.index(obs)], 2) 
    if obs.getUnit(): yTitles.append( 'Candidates / ' + str(binWidth) + ' ' + obs.getUnit() )
    else:             yTitles.append( 'Candidates / ' + str(binWidth) )

#Set plot options      
markStyle = 8
markSize  = 0.5
CpPlotsKit.setLineColors( dict(total = kBlue , even=kRed, odd=kGreen+3, swave=kMagenta+3) )
CpPlotsKit.setLineStyles( dict(total = kSolid, even=9   , odd=7       , swave=5         ) )
CpPlotsKit.setLineWidth(4)

##Plot and Save
for ( pad, obs, nBins, xTitle, yTitle, yScale, logY )\
        in zip(  [ TCanvas(o.GetName()) for o in observables ]
               , observables
               , numBins
               , ( angleNames[0][1], angleNames[1][1], angleNames[2][1], 'B_{s}^{0} decay time [ps]' )
               , yTitles
               , 3 * ( ( None, 1400 ), ) + ( ( 0.1, 10e4 ), )
               , 3 * ( False, ) + ( True, )
                ) :
    print '\n\n\n Ploting Observable ' + obs.GetName() +  ' {0}/{1} '.format(observables.index(obs)+1, len(observables)) + '\n\n\n'
    plot(  pad, obs, data, pdf, xTitle=xTitle, yTitle=yTitle, yScale=yScale, logy=logY
           , frameOpts   = dict( Bins = nBins, Name = obs.GetName() + 'Histo'   )
           , dataOpts    = dict( MarkerStyle = markStyle, MarkerSize = markSize )
           , pdfOpts     = CpPlotsKit.getPdfOpts(BinData=False) if 'time' in obs.GetName()\
                      else CpPlotsKit.getPdfOpts(BinData=True )
           , addPDFs     = CpPlotsKit.getAddPdfs()
           , addPDFsOpts = CpPlotsKit.getAddPdfsOpts(BinData=False) if 'time' in obs.GetName()\
                      else CpPlotsKit.getAddPdfsOpts(BinData=True )
           )
    lhcbName.Draw()
    filename = obs.GetName() + '_sFit.ps' if pdfConfig['SFit'] else obs.GetName() + '_cFit.ps'
    pad.Print(filename)

#Save all the plots in a root file as RooPlot objects.
from P2VV.Utilities.Plotting import _P2VVPlotStash as rooplots
from ROOT import TFile
plotsFile = TFile('RooPlotsFinal.root','recreate')
for plot in rooplots: plot.Write()
plotsFile.Close()



