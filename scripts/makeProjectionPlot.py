###########################################################################################################################################
## set script parameters ##
###########################

plotVar = 'timeLin'  # 'timeLin' / 'timeLog' / 'ctk' / 'ctl' / 'phi'
plotComp = 'total'   # 'total' / 'phi' / 'even' / 'odd' / 'S'

import sys
if len(sys.argv) > 1 :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('plotVar')
    parser.add_argument('plotComp')

    args = parser.parse_args()
    plotVar = args.plotVar
    plotComp = args.plotComp

assert plotVar  in [ 'timeLin', 'timeLog', 'ctk', 'ctl', 'phi' ]
assert plotComp in [ 'total', 'phi', 'even', 'odd', 'S' ]

print 'job parameters:'
print '  plot variable:  %s' % plotVar
print '  plot component: %s' % plotComp

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
#parFileIn   = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/parVals.par'
parFileIn   = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/LHCb/parVals.par'
dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetName = 'JpsiKK_sigSWeight'
#dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140822_fromNTuple_v2_corrWeights.root'
#plotsFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/projPlot_%s_%s.ps' % ( plotVar, plotComp )
#ROOTFilePath  = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/projPlot_%s_%s.root' % ( plotVar, plotComp )
plotsFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/LHCb/projPlot_%s_%s.ps' % ( plotVar, plotComp )
ROOTFilePath  = '/project/bfys/jleerdam/softDevel/P2VV2/test/plots/Reco14/LHCb/projPlot_%s_%s.root' % ( plotVar, plotComp )

# PDF options
#pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = dataPath + 'timeAcceptanceFit_2011.root'
#pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = dataPath + 'timeAcceptanceFit_2012.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = '/project/bfys/jleerdam/softDevel/P2VV2/test/fitResults/Reco14/effWeights/timeAcceptanceFit_2011.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = '/project/bfys/jleerdam/softDevel/P2VV2/test/fitResults/Reco14/effWeights/timeAcceptanceFit_2012.root'
pdfConfig['anglesEffType'] = 'basisSig6'
pdfConfig['angEffMomsFiles'] = dataPath + 'angEffNominalRew_5thOrder.par'

pdfConfig['lambdaCPParam'] = 'lambPhi'

# variables
varNames = [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ]
varDict = dict( timeLin = 'time', timeLog = 'time', ctk = 'helcosthetaK', ctl = 'helcosthetaL', phi = 'helphi' )
numBins = dict( timeLin = 47,     timeLog = 69,     ctk = 40,             ctl = 40,             phi = 40       )
parDict = dict( A0Mag2 = 'A0Mag2', AparMag2 = 'AparMag2', AperpMag2 = 'AperpMag2', f_S = 'f_S' )
rangeName = 'Bulk' if plotVar == 'timeLin' else ''

# plot options
from ROOT import gStyle, kBlack, kBlue, kRed, kGreen, kMagenta, kSolid, kDashed, kFullCircle, kFullSquare, kFullDotLarge
from P2VV.Load import LHCbStyle
#gStyle.SetEndErrorSize(3)
lineWidth     = 3
lineColor     = kBlue if plotComp == 'total' else kRed if plotComp in [ 'phi', 'even' ] else kGreen+3 if plotComp == 'odd' else kMagenta+3
lineStyle     = kSolid if plotComp == 'total' else 7 if plotComp in [ 'phi', 'odd' ] else 9 if plotComp == 'even' else 5
markStyle     = kFullDotLarge
markSize      = 0.7
markColor     = kBlack
markLineWidth = 3
gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
gStyle.SetLineStyleString( 7, ' 40 20'        )
gStyle.SetLineStyleString( 9, ' 100 20'       )

obsXTitles = dict(  timeLin = 'Decay time [ps]'
                  , timeLog = 'Decay time [ps]'
                  , ctk     = 'cos#kern[0.1]{#theta_{K}}'
                  , ctl     = 'cos#kern[0.1]{#theta_{#mu}}'
                  , phi     = '#varphi_{h} [rad]'
                 )
obsYTitles = dict(  timeLin = 'Candidates / (%.2g ps)'
                  , timeLog = 'Candidates / (%.2g ps)'
                  , ctk     = 'Candidates / %.2g'
                  , ctl     = 'Candidates / %.2g'
                  , phi     = 'Candidates / (%.2g#kern[0.3]{#pi} rad)'
                 )
xTitleOffset = 1.10
yTitleOffsets = dict(  timeLin = 1.15
                     , timeLog = 1.00
                     , ctk     = 1.15
                     , ctl     = 1.15
                     , phi     = 1.15
                    )
LHCbLabel = 'LHCb'
LHCbCoords = dict(  timeLin = ( 0.88, 0.85 )
                  , timeLog = ( 0.88, 0.85 )
                  , ctk     = ( 0.80, 0.85 )
                  , ctl     = ( 0.88, 0.85 )
                  , phi     = ( 0.88, 0.85 )
                 )


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

# build the PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)

# print parameters
print 120 * '='
print 'makeProjectionPlot: data:'
dataSet.Print()
print 'makeProjectionPlot: observables in PDF:'
pdfObs.Print('v')
print 'makeProjectionPlot: parameters in PDF:'
pdfPars.Print('v')

## fit data
#from P2VV.Utilities.DataHandling import correctWeights
#fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )
#fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = 8, Optimize = 2, Timer = True, Minimizer = 'Minuit2'
#                      , Strategy = 1, Offset = True )
#from P2VV.Imports import parNames, parValues
#fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )


###########################################################################################################################################
## build PDF for plotting ##
############################

if plotComp in [ 'phi', 'even', 'odd', 'S' ] :
    from ROOT import RooRealVar, RooConstVar, RooCustomizer
    from P2VV.RooFitWrappers import __dref__
    pdfCust = RooCustomizer( __dref__(pdf), plotComp )
    zeroCust = RooConstVar( 'zeroCust', 'zeroCust', 1.e-6 )

    if plotComp == 'even' :
        AparMag2Cust = RooRealVar( 'AparMag2Cust', 'AparMag2Cust', ws[ parDict['AparMag2'] ].getVal(), 0., 1. )
        pdfCust.replaceArg( ws[ parDict['AparMag2'] ], AparMag2Cust )
        pdfCust.replaceArg( ws[ parDict['AperpMag2'] ], zeroCust )

    elif plotComp == 'odd' :
        pdfCust.replaceArg( ws[ parDict['AparMag2'] ], zeroCust )
        pdfCust.replaceArg( ws[ parDict['A0Mag2'] ], zeroCust )

    if plotComp == 'S' :
        pdfCust.replaceArg( ws[ parDict['AparMag2'] ], zeroCust )
        pdfCust.replaceArg( ws[ parDict['AperpMag2'] ], zeroCust )
        pdfCust.replaceArg( ws[ parDict['A0Mag2'] ], zeroCust )

    else :
        if pdfConfig['paramKKMass'] :
            fSNames = [ '%s_bin%d' % ( parDict['f_S'], bin ) for bin in range( pdfBuild['KKMassBinning'].numBins() ) ]
        else :
            fSNames = [ parDict['f_S'] ]
        for name in fSNames :
            pdfCust.replaceArg( ws[name], zeroCust )

    func = pdfCust.build()

else :
    func = pdf

from P2VV.RooFitWrappers import SimultaneousPdf
projDataObs = [ obs for obs in pdf.ConditionalObservables() ]
if isinstance( pdf, SimultaneousPdf ) :
    if pdf.indexCat().isFundamental() :
        projDataObs += [ pdf.indexCat() ]
    else :
        projDataObs += [ cat for cat in pdf.indexCat().inputCatList() ]

if projDataObs :
    projData = dataSet.reduce( ArgSet = projDataObs )
    print 'makeProjectionPlot: plot projection data set:'
    projData.Print()
else :
    projData = None

pdfObs = dataSet.get()
from P2VV.RooFitWrappers import ExplicitNormPdf
obs = pdfObs.find( varDict[plotVar] )
xMin = obs.getMin(rangeName)
xMax = obs.getMax(rangeName)
normFac = ( xMax - xMin ) * dataSet.sumEntries() / float( numBins[plotVar] )
plotPdf = ExplicitNormPdf(  Name           = '%sPdf_%s' % ( plotVar, plotComp )
                          , Observables    = [ obs ]
                          , IntegrationObs = [ pdfObs.find(name) for name in varNames if name != varDict[plotVar] ]
                          , Function       = func
                          , NormFunction   = pdf
                          , NormFactor     = normFac
                          , ProjectionData = projData
                         )


###########################################################################################################################################
## make plot ##
###############

from ROOT import TFile, TCanvas
plotsFile = TFile.Open( ROOTFilePath, 'RECREATE' )
canv = TCanvas( '%sCanv' % plotVar )
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)

print 'makeProjectionPlot: plotting variable "%s" for component "%s"' % ( plotVar, plotComp )
binWidth = ( xMax - xMin ) / float( numBins[plotVar] ) / ( pi if plotVar == 'phi' else 1. )
from P2VV.Utilities.Plotting import plot
#redDataSet = dataSet.reduce( ArgSet = [ ws['time'], ws['helcosthetaK'], ws['helcosthetaL'], ws['helphi'] ] )
#plots = plot(  canv, obs, redDataSet, plotPdf
plots = plot(  canv, obs, dataSet, plotPdf
             , xTitle       = obsXTitles[plotVar]
             , yTitle       = obsYTitles[plotVar] % binWidth, logy = True if plotVar == 'timeLog' else False
             , xTitleOffset = xTitleOffset
             , yTitleOffset = yTitleOffsets[plotVar]
             , frameOpts    = dict( Range = rangeName, Bins = numBins[plotVar], Title = '' )
             , dataOpts     = dict( MarkerStyle = markStyle, MarkerSize = markSize, LineWidth = markLineWidth )
             , pdfOpts      = dict( LineColor = lineColor, LineWidth = lineWidth, LineStyle = lineStyle )
            )

from ROOT import TLatex
label = TLatex()
label.SetTextAlign(32)
label.SetTextSize(0.072)
drawLabel = label.DrawLatexNDC( LHCbCoords[plotVar][0], LHCbCoords[plotVar][1], LHCbLabel )
plotsFile.Add(drawLabel)

canv.Print(plotsFilePath)
from ROOT import TObject
plotsFile.Write( ROOTFilePath, TObject.kOverwrite )
plotsFile.Close()
