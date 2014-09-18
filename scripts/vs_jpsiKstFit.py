import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-B', '--BsFit',         default = 's',   type = str )
parser.add_argument( '-c', '--NumCpu',        default = 1,     type = int )   
parser.add_argument( '-s', '--sWave',         default = False, action = 'store_true'  )               
parser.add_argument( '-d', '--FitData',       default = False, action = 'store_true'  )
parser.add_argument( '-a', '--AngAccFile',    default = False, action = 'store_true'  )
parser.add_argument( '-i', '--ParFileIn',     default = False, action = 'store_true'  )
parser.add_argument( '-o', '--ParFileOut',    default = False, action = 'store_true'  )
parser.add_argument( '-b', '--writeUnBlPar',  default = True,  action = 'store_false' )
parser.add_argument( '-u', '--unblind',       default = False, action = 'store_true'  )
options = parser.parse_args()

#####################################################################################################################
## configuration  ##
############################

bsFit = options.BsFit

# dataset
path = '/project/bfys/vsyropou/data/Bs2JpsiKst'
dataSetFile = path + '/RealData/P2VVDataSet_20112012Reco14_Bs2JpsiKst_allKaons_fitNtuple_120614_%s_weighted.root'%('Bs'if bsFit else 'Bd')
dataSetName = 'jpsiKst'

# read / write fited parameters from file
parFileIn  = options.ParFileIn if options.ParFileIn else path + 'nominalFitResults/<writedetails>.par'
parFileOut = options.ParFileOut if options.ParFileOut else '<detailsoutput>.par'

# PDF configuration
from P2VV.Parameterizations.FullPDFs import Bs2JpsiKst_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# Bd/Bs weights
pdfConfig['sWeights'] = 'Bs' if bsFit else 'Bd'

# KK mass binning
# pdfConfig[''] = 

# swich blidn on/off
if options.unblind: pdfConfig['blind'] = {}

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
dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False, ImportIntoWS = False )
pdfConfig['signalData'] = dataSet
pdfConfig['readFromWS'] = False

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

## TODO:update csp
