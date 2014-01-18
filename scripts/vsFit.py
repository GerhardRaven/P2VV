from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d', '--FitData',    dest='FitData',    default='',   help='' )
parser.add_option('-a', '--AngAccFile', dest='AngAccFile', default='',   help='' )
parser.add_option('-i', '--ParFileIn',  dest='ParFileIn',  default='',   help='' )
parser.add_option('-o', '--ParFileOut', dest='ParFileOut', default='',   help='' )
(options, args) = parser.parse_args()

# data
dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
dataSetFile = options.FitData if options.FitData else dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
dataSetName = 'JpsiKK_sigSWeight'

myPath = '/project/bfys/vsyropou/data/'

# read / write fited parameters from file
parFileIn = options.ParFileIn if options.ParFileIn \
    else myPath + 'nominalFitResults/20112012Reco14DataFitValues_6KKMassBins.par'
parFileOut = options.ParFileOut if options.ParFileOut else '20112012Reco14DataFitValues_6KKMassBins.par'

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


pdfConfig['timeEffParameters'] = dict()

timeEffFile2011 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
timeEffFile2012 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = options.AngAccFile if options.AngAccFile \
    else myPath + 'uncorrecteEffMoments/MC20112012_Sim08/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_weights'

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

# # fix values of some parameters
# for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
#         if not pdfConfig['SSTagging'] :
# 		CEvenOdds.setConstant( 'avgCEven.*')
#                 CEvenOdds.setConstant( 'avgCOdd.*', True )
#         else :
#             for CEvenOdd in CEvenOdds :
# 		    CEvenOdd.setConstant('avgCEven.*')
#                     CEvenOdd.setConstant( 'avgCOdd.*', True )
# pdfBuild['amplitudes'].setConstant('C_SP')

if parFileIn:
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)
print '-' * 80 + '\n\n'

# fit data
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, **fitOpts )
from P2VV.Imports import parNames, parValues
fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )
fitResult.SetName( parFileOut.replace('.par','') )
from ROOT import TFile
resultFile = TFile.Open( parFileOut.replace('.par','.root'), 'recreate')
resultFile.cd()
fitResult.Write()
resultFile.Close()

if parFileOut :
    pdfConfig.getParametersFromPdf( pdf,  fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )

print 120*'='
