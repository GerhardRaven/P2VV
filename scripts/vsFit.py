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

# read / write fited parameters from file
parFileIn = options.ParFileIn if options.ParFileIn \
    else '/project/bfys/jleerdam/softDevel/P2VV2/test/20112012Reco14DataFitValues_6KKMassBins.par'
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

pdfConfig['splitParams']['hlt1_excl_biased_dec'] = [ 'tagCatCoef0_1' ]
pdfConfig['timeEffParameters'] = dict()

timeEffFile2011 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
timeEffFile2012 = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeEffFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeEffFile2012
pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = options.AngAccFile if options.AngAccFile \
    else dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

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
# for CEvenOdds in self._pdfBuild['taggingParams']['CEvenOdds'] :
#         if not self._pdfConfig['SSTagging'] :
#             if namePF:
#                 CEvenOdds.setConstant( 'avgCEven.*')
#                 CEvenOdds.setConstant( 'avgCOdd.*', True )
#             else: 
#                 CEvenOdds.setConstant( 'avgCEven.*')
#                 CEvenOdds.setConstant( 'avgCOdd.*', True )
#         else :
#             for CEvenOdd in CEvenOdds :
#                 if namePF:
#                     CEvenOdd.setConstant( 'avgCEven.*')
#                     CEvenOdd.setConstant( 'avgCOdd.*', True )
#                 else: 
#                     CEvenOdd.setConstant('avgCEven.*')
#                     CEvenOdd.setConstant( 'avgCOdd.*', True )
# self._pdfBuild['amplitudes'].setConstant('C_SP')

if parFileIn:
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# print PDF values
from ROOT import RooArgSet
angSet = RooArgSet( ws[var] for var in [ 'helcosthetaK', 'helcosthetaL', 'helphi' ] )
timeAngSet = RooArgSet( ws[var] for var in [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ] )
print '\n\n' + '-' * 80
for period in [ 2011, 2012 ] :
    ws['runPeriod'].setIndex(period)
    for cat in [ 0, 1 ] :
        ws['hlt1_excl_biased_dec'].setIndex(cat)
        print 'PDF values "%d"/"%s":' % ( ws['runPeriod'].getIndex(), ws['hlt1_excl_biased_dec'].getLabel() )
        print 'unnormalized:', pdf.getVal()
        print 'angle-normalized:', pdf.getVal(angSet)
        print 'time-angle-normalized:', pdf.getVal(timeAngSet)
        print
print '-' * 80 + '\n\n'

if fit :
    # fit data
    fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, **fitOpts )
    from P2VV.Imports import parNames, parValues
    fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )

if parFileOut :
    pdfConfig.getParametersFromPdf( pdf,  fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut )
