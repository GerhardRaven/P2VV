fit = True
useTimeResMean = True

dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetFile = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
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

pdfConfig['splitParams']['hlt1_excl_biased_dec'] = [ 'tagCatCoef0_1' ]
pdfConfig['timeEffParameters'] = dict( Parameterization = 'Spline', Fit = False )

timeEffFile2011 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
timeEffFile2012 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

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

# print PDF values
from ROOT import RooArgSet
angSet = RooArgSet( ws[var] for var in [ 'helcosthetaK', 'helcosthetaL', 'helphi' ] )
timeAngSet = RooArgSet( ws[var] for var in [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ] )
print '\n\n' + '-' * 80
print 'PDF values:'
print 'unnormalized:', pdf.getVal()
print 'angle-normalized:', pdf.getVal(angSet)
print 'time-angle-normalized:', pdf.getVal(timeAngSet)
print '-' * 80 + '\n\n'

if fit :
    # fit data
    fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, **fitOpts )
    from P2VV.Imports import parNames, parValues
    fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )
