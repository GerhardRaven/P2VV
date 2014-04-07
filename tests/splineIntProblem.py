getPdfVals = False
createNLL = False
doFit = False
useTimeResMean = False

dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'
dataSetName = 'JpsiKK_sigSWeight'
parFileIn   = '/project/bfys/jleerdam/softDevel/P2VV2/test/20112012Reco14DataFitValues_6KKMassBins.par'

fitOpts = dict(  NumCPU    = 8
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Strategy  = 1
               , Offset    = True
              )

# PDF configuration
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig( RunPeriods = '3fb' )

pdfConfig['splitParams']['hlt1_excl_biased'] = [ 'tagCatCoef0_1' ]
pdfConfig['timeEffParameters'] = dict( Parameterization = 'Spline', Fit = False )

timeEffFile2011 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins_TOS.root'
timeEffFile2012 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins_TOS.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = dataPath + 'angEffNominalRew_moms.par'

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

# data set with weights corrected for background dilution: for phi_s fit only!
from P2VV.Utilities.DataHandling import correctWeights
fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )

# build PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

if not useTimeResMean :
    # remove time resolution offset
    ws['timeResMu_p2011'].setVal(0.)
    ws['timeResMu_p2012'].setVal(0.)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)

# print parameters
print 120 * '='
print 'splineIntProblem: data:'
dataSet.Print()
print 'splineIntProblem: observables in PDF:'
pdfObs.Print('v')
print 'splineIntProblem: parameters in PDF:'
pdfPars.Print('v')
print 'splineIntProblem: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()

# get PDF value at large decay time
pdfObs.setRealValue( 'time', 13. )
pdf.getVal()

# plot time-acceptance functions
from P2VV.Load import LHCbStyle
from ROOT import TCanvas
canv = TCanvas('canv')
canv.SetLeftMargin(0.18)
canv.SetRightMargin(0.05)
canv.SetBottomMargin(0.18)
canv.SetTopMargin(0.05)
canv.Divide( 2, 2 )
for yIt, year in enumerate( [ 2011, 2012 ] ) :
    for tIt, trig in enumerate( [ 'notExclB', 'exclB' ] ) :
        accFun = pdfBuild['timeResModels'][ '{bin0;p%d;%s}' % ( year, trig ) ]._shape
        plot = pdfObs.find('time').frame()
        accFun.plotOn(plot)
        plot.SetXTitle('Decay time [ps]')
        plot.SetYTitle( 'HLT1 %s acc. %d' % ( 'unbiased' if trig == 'notExclB' else 'excl. biased', year ) )
        plot.SetTitleOffset( 1.1, 'x' )
        plot.SetTitleOffset( 1.0, 'y' )
        plot.SetMinimum(-1.)
        plot.SetMaximum(2.)
        canv.cd( 2 * tIt + yIt + 1 )
        plot.Draw()
canv.Print('acceptance.ps')

if getPdfVals :
    # get PDF values for events in data
    for it, obsSet in enumerate(fitData) :
        pdfObs.assignValueOnly(obsSet)
        pdfVal = pdf.getVal()
        if pdfVal <= 0. : break

if createNLL :
    # create NLL
    nll = pdf.createNLL(fitData)
    nll.getVal()

if doFit :
    # fit data
    fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, **fitOpts )
    from P2VV.Imports import parNames, parValues
    fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )
