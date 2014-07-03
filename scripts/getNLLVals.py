###########################################################################################################################################
## set script parameters ##
###########################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobID', '-I', type = int, default = 1000000 )
parser.add_argument( '--scanPar', '-p', required = True )
parser.add_argument( '--scanParVals', '-v', required = True, nargs = '+', type = float )
parser.add_argument( '--model', '-m', default = 'lamb_phi' )  # 'phi' / 'lamb_phi' / 'polarDep'
parser.add_argument( '--blind', '-b', default = True )
parser.add_argument( '--numCPU', '-c', type = int, default = 2 )
parser.add_argument( '--dataPath', '-d', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14' )
parser.add_argument( '--workPath', '-w', default = '/project/bfys/jleerdam/softDevel/P2VV2/test' )
parser.add_argument( '--dataSetFile', '-f', default = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root' )
parser.add_argument( '--dataSetName', '-n', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--parFileIn', '-i' )
parser.add_argument( '--timeAccFile2011', '-x', default = 'timeAcceptanceFit_2011.root' )
parser.add_argument( '--timeAccFile2012', '-y', default = 'timeAcceptanceFit_2012.root' )
parser.add_argument( '--angAccType', '-t', default = 'weights' )  # 'weights' / 'basis012' / 'basis01234' / 'basisSig6'
parser.add_argument( '--angAccFile', '-z', default = 'angEffNominalRew_moms.par' )

args = parser.parse_args()
assert args.model in [ 'phi', 'lamb_phi', 'polarDep' ]
blind = False if not args.blind or str( args.blind ).lower() in [ 'false', '0' ] else True
assert type(args.numCPU) == int and args.numCPU > 0 and args.numCPU < 20
dataPath = args.dataPath
workPath = args.workPath
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
dataSetFile = dataPath + args.dataSetFile
parFileIn = args.parFileIn
if not parFileIn :
    parFileIn = workPath + '20112012Reco14DataFitValues_6KKMassBins%s.par'\
                           % ( '_CPVDecay' if args.model == 'polarDep' else '_fixedLamb' if args.model == 'phi' else '' )
else :
    parFileIn = workPath + parFileIn
timeAccFile2011 = dataPath + args.timeAccFile2011
timeAccFile2012 = dataPath + args.timeAccFile2012
assert args.angAccType in [ 'weights', 'basis012', 'basis01234', 'basisSig6' ]
angAccFile = dataPath + args.angAccFile

print 'job parameters:'
print '  job ID: %s' % args.jobID
print '  blind analysis: %s' % ( 'true' if blind else 'false' )
print '  scan parameter: %s' % args.scanPar
print '  scan parameter values: %s' % args.scanParVals
print '  model: %s' % args.model
print '  number of cores: %d' % args.numCPU
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  dataset file: %s' % dataSetFile
print '  dataset name: %s' % args.dataSetName
print '  input parameter file: %s' % parFileIn
print '  time acceptance file 2011: %s' % timeAccFile2011
print '  time acceptance file 2012: %s' % timeAccFile2012
print '  angular acceptance type: %s' % args.angAccType
print '  angular acceptance file: %s' % angAccFile

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()
if not blind :
    pdfConfig['blind'] = { }

# PDF options
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = dataPath + args.timeAccFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = dataPath + args.timeAccFile2012
pdfConfig['anglesEffType'] = args.angAccType
pdfConfig['angEffMomsFiles'] = dataPath + args.angAccFile

pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if args.model == 'polarDep' else 'lambPhi'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataSetFile, dataSetName = args.dataSetName,  NTuple = False )
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

# fix of float |lambda|
if args.model == 'phi' :
    ws['lambdaCP'].setVal(1.)# 0.96340)
    ws['lambdaCP'].setConstant(True)
elif args.model == 'lamb_phi' :
    ws['lambdaCP'].setConstant(False)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)
assert args.scanPar in ws and ws[args.scanPar] in pdfPars

# prevent tagging dilution = 0 in initialization
if abs( 1. - pdfPars.getRealValue('wTagP1OS') ) < 1.e-5 and not pdfPars.find('wTagP1OS').isConstant() :
    pdfPars.setRealValue( 'wTagP1OS', 0.999 )
if abs( 1. - pdfPars.getRealValue('wTagP1SS') ) < 1.e-5 and not pdfPars.find('wTagP1SS').isConstant() :
    pdfPars.setRealValue( 'wTagP1SS', 0.999 )

# print parameters
print 120 * '='
print 'getNLLVals: data:'
dataSet.Print()
print 'getNLLVals: observables in PDF:'
pdfObs.Print('v')
print 'getNLLVals: parameters in PDF:'
pdfPars.Print('v')

# correct weights in data
from P2VV.Utilities.DataHandling import correctWeights
fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )


###########################################################################################################################################
## get NLL values ##
####################

print 'getNLLVals: create NLL variable'
nll = pdf.createNLL( fitData, NumCPU = args.numCPU, Optimize = 2 )
nll.getVal()
print '\nNLL value = %.3f' % nll.getVal()

# open output file for NLL values
nllFilePath = workPath + 'NLLVals_%s_%d.par' % ( args.scanPar, args.jobID )
try :
    nllFile = open( nllFilePath, 'w' )
except :
    raise RuntimeError( 'getNLLVals: ERROR: unable to open file "%s"' % nllFilePath )

# loop over values of scan parameter
for val in args.scanParVals :
    # set parameter values
    pdfConfig.setParametersInPdf(pdf)
    if ws[args.scanPar].getMin() > val :
        print 'getNLLVals: WARNING: setting minimum of parameter "%s" to %.5g' % ( args.scanPar, val )
        ws[args.scanPar].setMin(val)
    if ws[args.scanPar].getMax() < val :
        print 'getNLLVals: WARNING: setting maximum of parameter "%s" to %.5g' % ( args.scanPar, val )
        ws[args.scanPar].setMax(val)
    ws[args.scanPar].setVal(val)
    ws[args.scanPar].setConstant(True)
    print 'getNLLVals: %s = %.6g\n' % ( args.scanPar, val )
    nllVal = nll.getVal()

    # do profile likelihood fit
    fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = args.numCPU, Optimize = 2, Timer = True
                          , Minimizer = 'Minuit2', Strategy = 1, Offset = True, Hesse = False, Minos = False )
    profVal = nll.getVal()

    # write NLL values to file
    print '\ngetNLLVals: parameters in PDF:'
    pdfPars.Print('v')
    nllStr = '%s = %.6g : NLL = %.6f : profiled NLL = %.6f' % ( args.scanPar, val, nllVal, profVal )
    print '\ngetNLLVals: %s' % nllStr
    print 'writing NLL values to "%s"\n' % nllFilePath
    nllFile.write( nllStr + '\n' )
nllFile.close()
