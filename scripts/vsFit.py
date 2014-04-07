#####################################################################################################################
## configuration  ##
############################

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-d', '--FitData',      dest='FitData',      default='',        help='' )
parser.add_option('-a', '--AngAccFile',   dest='AngAccFile',   default='',        help='' )
parser.add_option('-i', '--ParFileIn',    dest='ParFileIn',    default='',        help='' )
parser.add_option('-o', '--ParFileOut',   dest='ParFileOut',   default='',        help='' )
parser.add_option('-b', '--writeUnBlPar', dest='writeUnBlPar', default='True',    help='' )
parser.add_option('-u', '--unblind',      dest='unblind',      default='False',   help='' )
parser.add_option('-c', '--NumCpu',       dest='NumCpu',       default=8,         help='', type=int )
(options, args) = parser.parse_args()

# data
dataPath    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/'
myPath = '/project/bfys/vsyropou/data/Bs2JpsiPhi/'
dataSetFile = options.FitData if options.FitData else \
              dataPath + 'angEff/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
dataSetName = 'JpsiKK_sigSWeight'

# read / write fited parameters from file
parFileIn = options.ParFileIn if options.ParFileIn \
    else myPath + 'nominalFitResults/corrAngAccDEC/20112012Reco14DataFitValues_6KKMassBins.root'
parFileOut = options.ParFileOut if options.ParFileOut else '20112012Reco14DataFitValues_6KKMassBins.par'

# PDF configuration
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# swich blidn on/off
if 'True' in options.unblind: pdfConfig['blind'] = {}

# cpv parametrization
pdfConfig['lambdaCPParam'] = 'lambPhi'

# acceptances
 # decay time
pdfConfig['timeEffType']   = 'paper2012'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
        = dataPath + 'Reco14/timeAcceptanceFit_2011.root' # Bs_HltPropertimeAcceptance_Data_2011_40bins_TOS.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
        = dataPath + 'Reco14/timeAcceptanceFit_2012.root' # Bs_HltPropertimeAcceptance_Data_2012_40bins_TOS.root'
 # decay angles
pdfConfig['anglesEffType']   = 'weights'
pdfConfig['angEffMomsFiles'] = options.AngAccFile if options.AngAccFile \
    else myPath   + 'angEffMoments/correctedEffMoms/DEC/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm_fromAna' # Iter. proc. corrected acc. of DEC dataset  
#   else myPath   + 'angEffMoments/uncorrecteEffMoments/TOS/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'     # Uncorrected acc of TOS dataset

#####################################################################################################################
## read data and build pdf ##
############################

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

if parFileIn:
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)
print '-' * 80 + '\n\n'

# fix-float lambda
ws['lambdaCP'].setConstant(False)


###########################################################################################################################################
## fit data ##
##############

# fit
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = options.NumCpu, Optimize = 2, Timer = True, 
                       Minimizer = 'Minuit2', Strategy = 1, Offset = True 
                       )

# save
from P2VV.Imports import parNames, parValues
fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )
fitResult.SetName( parFileOut.replace('.par','') )
from ROOT import TFile
resultFile = TFile.Open( parFileOut.replace('.par','.root'), 'recreate')
resultFile.cd()
fitResult.Write()
resultFile.Close()

if parFileOut:
    pdfConfig.getParametersFromPdf( pdf,  fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut, Floating = True )
    if 'True' in options.writeUnBlPar:
        for parName in ['phiCP', 'dGamma']:
            par = pdfConfig.parameters().pop('__%s__'%parName)
            pdfConfig.parameters()[parName] = (ws[parName].getVal(), ) + par[1:]
        pdfConfig.writeParametersToFile( filePath = parFileOut.replace('.par','_unbl.par'), Floating = True )
