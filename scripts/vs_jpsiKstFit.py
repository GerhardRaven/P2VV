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
parser.add_option('-c', '--NumCpu',       dest='NumCpu',       default=1,         help='', type=int )
(options, args) = parser.parse_args()

# data
path        = '/project/bfys/vsyropou/data/Bs2JpsiPhi/'
dataSetFile = options.FitData if options.FitData else path + 'MC/Bd_MCT_p_correctAngles.root'
dataSetName = 'T'

# read / write fited parameters from file
parFileIn  = options.ParFileIn if options.ParFileIn else path + 'nominalFitResults/<writedetails>.par'
parFileOut = options.ParFileOut if options.ParFileOut else '<detailsoutput>.par'

# PDF configuration
from P2VV.Parameterizations.FullPDFs import Bs2JpsiKst_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# swich blidn on/off
if 'True' in options.unblind: pdfConfig['blind'] = {}

# angular acceptance input
pdfConfig['anglesEffType']   = 'weights'
pdfConfig['angEffMomsFiles'] = options.AngAccFile if options.AngAccFile \
    else path   + 'angEffMoments/correctedEffMoms/DEC/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm_fromAna' # Iter. proc. corrected acc. of DEC dataset  

#####################################################################################################################
## read data and build pdf ##
############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiKstWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = True

## TODO:: make correct sWeights for jpsiKst.
# data set with weights corrected for background dilution
# from P2VV.Utilities.DataHandling import correctWeights
# fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )

# build PDF
from P2VV.Parameterizations.FullPDFs import Bs2JpsiKst_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

assert False
if parFileIn:
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)
print '-' * 80 + '\n\n'

###########################################################################################################################################
## fit data ##
##############

# fit
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = options.NumCpu, Optimize = 2, Timer = True, 
                       Minimizer = 'Minuit2', Strategy = 1, Offset = True 
                       )

## TODO: Make the imports for jpsiKst parameters
# save
# from P2VV.Imports import parNames, parValues
# fitResult.PrintSpecial( text = True, ParNames = parNames, ParValues = parValues )
# fitResult.SetName( parFileOut.replace('.par','') )
# from ROOT import TFile
# resultFile = TFile.Open( parFileOut.replace('.par','.root'), 'recreate')
# resultFile.cd()
# fitResult.Write()
# resultFile.Close()

if parFileOut:
    pdfConfig.getParametersFromPdf( pdf,  fitData )
    pdfConfig.writeParametersToFile( filePath = parFileOut, Floating = True )
    if 'True' in options.writeUnBlPar:
        for parName in ['phiCP', 'dGamma']:
            par = pdfConfig.parameters().pop('__%s__'%parName)
            pdfConfig.parameters()[parName] = (ws[parName].getVal(), ) + par[1:]
        pdfConfig.writeParametersToFile( filePath = parFileOut.replace('.par','_unbl.par'), Floating = True )
