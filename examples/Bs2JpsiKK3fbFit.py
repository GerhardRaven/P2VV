###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobName', '-N', default = 'Bs2JpsiKKFit' )
parser.add_argument( '--model', '-m', default = 'lamb_phi' )  # 'phi' / 'lamb_phi' / 'polarDep'
parser.add_argument( '--fixLowAcc', '-l', default = True )
parser.add_argument( '--fixUpAcc', '-u', default = False )
parser.add_argument( '--numCPU', '-c', type = int, default = 2 )
parser.add_argument( '--dataPath', '-d', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14' )
parser.add_argument( '--workPath', '-w', default = '/project/bfys/jleerdam/softDevel/P2VV2/test' )
parser.add_argument( '--dataSetFile', '-f', default = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root' )
parser.add_argument( '--accDataSetFile', '-af', default =  'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root' )
parser.add_argument( '--dataSetName', '-n', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--parFileIn', '-i' )
parser.add_argument( '--parFileOut', '-o' )
parser.add_argument( '--timeAccFile2011', '-t11', default = 'timeAcceptanceFit_2011.root' )
parser.add_argument( '--timeAccFile2012', '-t12', default = 'timeAcceptanceFit_2012.root' )
parser.add_argument( '--angAccFile', '-a', default = 'angEffNominalRew_moms.par' )

args = parser.parse_args()
assert args.model in [ 'phi', 'lamb_phi', 'polarDep' ]
fixLowAcc = False if not args.fixLowAcc or str( args.fixLowAcc ).lower() in [ 'false', '0' ] else True
fixUpAcc = False if not args.fixUpAcc or str( args.fixUpAcc ).lower() in [ 'false', '0' ] else True
assert type(args.numCPU) == int and args.numCPU > 0 and args.numCPU < 20
dataPath = args.dataPath
workPath = args.workPath
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
dataSetFile = dataPath + args.dataSetFile
accDataSetFile = dataPath + args.accDataSetFile
parFileIn = args.parFileIn
if parFileIn == None :
    parFileIn = workPath + '20112012Reco14DataFitValues_6KKMassBins%s.par' % ( '_CPVDecay' if args.model == 'polarDep' else '' )
elif parFileIn :
    parFileIn = workPath + parFileIn
parFileOut = args.parFileOut
if parFileOut == None :
    parFileOut = workPath + '%s%s_%sLow_%sUp.par'\
                 % ( args.jobName + ( '_' if args.jobName else '' ), args.model, 'fix' if fixLowAcc else 'float'
                    , 'fix' if fixUpAcc else 'float' )
elif parFileOut :
    parFileOut = workPath + parFileOut
timeAccFile2011 = dataPath + args.timeAccFile2011
timeAccFile2012 = dataPath + args.timeAccFile2012
angAccFile = dataPath + args.angAccFile

# print script settings
print 'job parameters:'
if args.jobName :
    print '  job name: %s' % args.jobName
print '  model: %s' % args.model
print '  fix lower decay-time acceptance: %s' % ( 'true' if fixLowAcc else 'false' )
print '  fix upper decay-time acceptance: %s' % ( 'true' if fixUpAcc  else 'false' )
print '  number of cores: %d' % args.numCPU
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  dataset file: %s' % dataSetFile
print '  acceptance dataset file: %s' % accDataSetFile
print '  dataset name: %s' % args.dataSetName
print '  input parameter file: %s' % parFileIn
print '  output parameter file: %s' % parFileOut
print '  time acceptance file 2011: %s' % timeAccFile2011
print '  time acceptance file 2012: %s' % timeAccFile2012
print '  angular acceptance file: %s' % angAccFile

# PDF options
from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()
pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if args.model == 'polarDep' else 'lambPhi'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]

pdfConfig['timeEffType'] = 'paper2012' if fixLowAcc else 'fit_uniformUB'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeAccFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeAccFile2012
if pdfConfig['timeEffType'].startswith('fit') :
    from P2VV.Parameterizations.FullPDFs import SimulCatSettings
    pdfConfig['timeEffData']['file'] = accDataSetFile
    pdfConfig['externalConstr']['acceptance'] = SimulCatSettings('acceptanceConstr')
    pdfConfig['externalConstr']['acceptance'].addSettings(  [ 'runPeriod' ], [ [ 'p2011', 'p2012' ] ]
                                                          , {  ( 'hlt1_excl_biased_dec', 'exclB' ) : ( 0.65, 0.01 )
                                                             , ( 'hlt2_biased', 'B' )              : ( 0.65, 0.01 )
                                                            }
                                                         )

if fixUpAcc :
    for it, sett in enumerate( pdfConfig['externalConstr']['betaTimeEff'] ) :
        per = sett[0]['runPeriod']
        assert per in [ [ 'p2011' ], [ 'p2012' ] ]
        pdfConfig['externalConstr']['betaTimeEff'][it] = ( sett[0], ( -0.008639, 0. ) if per == [ 'p2011' ] else ( -0.012669, 0. ) )

pdfConfig['anglesEffType'] = 'weights'
pdfConfig['angEffMomsFiles'] = angAccFile


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
    ws['lambdaCP'].setVal(1.)
    ws['lambdaCP'].setConstant(True)
elif args.model == 'lamb_phi' :
    ws['lambdaCP'].setConstant(False)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)

# print parameters
print 120 * '='
print 'Bs2JpsiKK3fbFit: data:'
dataSet.Print()
print 'Bs2JpsiKK3fbFit: observables in PDF:'
pdfObs.Print('v')
print 'Bs2JpsiKK3fbFit: parameters in PDF:'
pdfPars.Print('v')
print 'Bs2JpsiKK3fbFit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()


###########################################################################################################################################
## fit data ##
##############

# fit data
from P2VV.Utilities.DataHandling import correctWeights
fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = args.numCPU, Optimize = 2, Timer = True, Minimizer = 'Minuit2'
                      , Strategy = 1, Offset = True )

if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    from P2VV.Imports import parNames, parValuesCPVDecay as parValues
else :
    from P2VV.Imports import parNames, parValues
print 120 * '-'
print 'parameter values:'
fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
print 120 * '-'
print 'correlation matrix:'
fitResult.correlationMatrix().Print()
print 120 * '-'
print 'covariance matrix:'
fitResult.covarianceMatrix().Print()
print 120 * '-' + '\n'

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( pdf, fitData )
    pdfConfig.writeParametersToFile(  filePath = parFileOut
                                    , FitStatus = ( fitResult.status(), fitResult.minNll(), fitResult.edm() )
                                   )
