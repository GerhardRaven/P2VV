###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobName', '-N', default = 'Bs2JpsiphiMCFit' )
parser.add_argument( '--runPeriod', '-R', default = '2011' )
parser.add_argument( '--model', '-m', default = 'lamb_phi' )  # 'phi' / 'lamb_phi' / 'polarDep'
parser.add_argument( '--fixUpAcc', '-u', default = True )
parser.add_argument( '--trueTags', '-b', default = True )
parser.add_argument( '--fixTagging', '-a', default = False )
parser.add_argument( '--fixTagAsym', '-k', default = False )
parser.add_argument( '--numCPU', '-c', type = int, default = 2 )
parser.add_argument( '--runHesse', '-e', default = True )
parser.add_argument( '--runMinos', '-s', default = '' )
parser.add_argument( '--minosPars', '-p', default = None )
parser.add_argument( '--dataPath', '-d', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14' )
parser.add_argument( '--workPath', '-w', default = '/project/bfys/jleerdam/softDevel/P2VV2/test' )
parser.add_argument( '--dataSetFile', '-f', default = 'DS2011/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_000.root' )
parser.add_argument( '--dataSetName', '-n', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--parFileIn', '-i' )
parser.add_argument( '--parFileOut', '-o' )
parser.add_argument( '--resultFileOut', '-r' )
parser.add_argument( '--timeAccFile2011', '-x', default = 'DS2011/timeAcc_2011.root' )
parser.add_argument( '--timeAccFile2012', '-y', default = 'DS2012/timeAcc_2012.root' )
parser.add_argument( '--angAccType', '-t', default = 'weights' )  # 'weights' / 'basis012' / 'basis01234' / 'basisSig6'
parser.add_argument( '--angAccFile', '-z', default = 'DS2011/Sim08_2011_hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights' )
parser.add_argument( '--constAngAcc', '-q', default = True )

args = parser.parse_args()
assert args.model in [ 'phi', 'lamb_phi', 'polarDep' ]
assert args.runPeriod in [ '2011', '2012' ]
fixUpAcc = False if not args.fixUpAcc or str( args.fixUpAcc ).lower() in [ 'false', '0' ] else True
trueTags = False if not args.trueTags or str( args.trueTags ).lower() in [ 'false', '0' ] else True
fixTagging = False if ( not args.fixTagging or str( args.fixTagging ).lower() in [ 'false', '0' ] ) and not trueTags else True
fixTagAsym = False if ( not args.fixTagAsym or str( args.fixTagAsym ).lower() in [ 'false', '0' ] ) and not trueTags else True
assert type(args.numCPU) == int and args.numCPU > 0 and args.numCPU < 20
runHesse = False if not args.runHesse or str( args.runHesse ).lower() in [ 'false', '0' ] else True
minosPars = args.minosPars.split(',') if args.minosPars and str(args.minosPars) != 'None' else [ ]
runMinos = False if ( args.runMinos == '' and not minosPars ) or str( args.runMinos ).lower() in [ 'false', '0' ] else True
dataPath = args.dataPath
workPath = args.workPath
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
dataSetFile = dataPath + args.dataSetFile
parFileIn = args.parFileIn
if parFileIn == None or parFileIn == 'None' :
    parFileIn = workPath + 'MCParValues.par'
elif parFileIn :
    parFileIn = workPath + parFileIn
parFileOut = args.parFileOut
if parFileOut == None or parFileOut == 'None' :
    parFileOut = workPath + '%s%s.par' % ( args.jobName + ( '_' if args.jobName else '' ), args.model )
elif parFileOut :
    parFileOut = workPath + parFileOut
resultFileOut = args.resultFileOut
if resultFileOut == None or resultFileOut == 'None' :
    resultFileOut = workPath + '%s%s.root' % ( args.jobName + ( '_' if args.jobName else '' ), args.model )
elif resultFileOut :
    resultFileOut = workPath + resultFileOut
timeAccFile2011 = dataPath + args.timeAccFile2011
timeAccFile2012 = dataPath + args.timeAccFile2012
assert args.angAccType in [ 'weights', 'basis012', 'basis01234', 'basisSig6' ]
angAccFile = dataPath + args.angAccFile
constAngAcc = False if not args.constAngAcc or str( args.constAngAcc ).lower() in [ 'false', '0' ] else True

# print script settings
print 'job parameters:'
if args.jobName :
    print '  job name: %s' % args.jobName
print '  model: %s' % args.model
print '  runPeriod: %s' % args.runPeriod
print '  fix upper decay-time acceptance: %s' % ( 'true' if fixUpAcc  else 'false' )
print '  true flavour tags: %s' % ( 'true' if trueTags else 'false' )
print '  fix tagging calibration: %s' % ( 'true' if fixTagging else 'false' )
print '  fix tagging calibration asymmetries: %s' % ( 'true' if fixTagAsym else 'false' )
print '  number of cores: %d' % args.numCPU
print '  run Hesse: %s' % ( 'true' if runHesse else 'false' )
print '  run Minos: %s' % ( 'true' if runMinos else 'false' )
print '  Minos parameters:', ', '.join( '"%s"' % par for par in minosPars )
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  dataset file: %s' % dataSetFile
print '  dataset name: %s' % args.dataSetName
print '  input parameter file: %s' % parFileIn
print '  output parameter file: %s' % parFileOut
print '  output fit result file: %s' % resultFileOut
print '  time acceptance file 2011: %s' % timeAccFile2011
print '  time acceptance file 2012: %s' % timeAccFile2012
print '  angular acceptance type: %s' % args.angAccType
print '  angular acceptance file: %s' % angAccFile
print '  constant angular acceptance parameters: %s' % ( 'true' if constAngAcc else 'false' )

# clear command-line options
import sys
sys.argv = sys.argv[ 0 : 1 ]

# PDF options
from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig( RunPeriods = args.runPeriod )
pdfConfig['blind'] = { }
pdfConfig['splitParams'] = { }
pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if args.model == 'polarDep' else 'lambPhi'

pdfConfig['paramKKMass']     = ''
pdfConfig['KKMassBinBounds'] = [ 990., 1050. ]
pdfConfig['CSPValues']       = [ 1. ]

pdfConfig['timeEffParameters'] = dict( Cache = False )
pdfConfig['timeEffHistFiles']['file'] = timeAccFile2011 if args.runPeriod == '2011' else timeAccFile2012
if args.runPeriod == '2011' :
    #timeAccHistHLT1UB2011 = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
    timeAccHistHLT1UB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
    pdfConfig['timeEffHistFiles']['hlt1UB'] = timeAccHistHLT1UB

#pdfConfig['timeResType'] = ''
#pdfConfig['timeEffParameters'] = dict( Cache = False )
#for par in [ 'sf_placeholder', 'timeResMu', 'timeResFrac2', 'sf_mean_slope', 'sf_mean_quad', 'sf_sigma_slope', 'sf_sigma_quad' ] :
#    pdfConfig['externalConstr'].pop(par)
pdfConfig['timeResType'] += '_mean_param'
pdfConfig['externalConstr']['mu_placeholder']   = (  0.0350  if args.runPeriod == '2011' else  0.0349,  0. )
pdfConfig['externalConstr']['timeResMu_offset'] = ( -0.00174 if args.runPeriod == '2011' else -0.00169, 0. )
pdfConfig['externalConstr']['timeResMu_slope']  = ( -0.1350  if args.runPeriod == '2011' else -0.127,   0. )
pdfConfig['externalConstr']['timeResMu_quad']   = ( -3.50    if args.runPeriod == '2011' else -2.29,    0. )
pdfConfig['externalConstr']['sf_placeholder']   = (  0.0350  if args.runPeriod == '2011' else  0.0349,  0. )
pdfConfig['externalConstr']['timeResFrac2']     = (  0.191   if args.runPeriod == '2011' else  0.229,   0. )
pdfConfig['externalConstr']['sf_mean_slope']    = (  1.239   if args.runPeriod == '2011' else  1.2766,  0. )
pdfConfig['externalConstr']['sf_mean_quad']     = (  1.98    if args.runPeriod == '2011' else  1.34,    0. )
pdfConfig['externalConstr']['sf_sigma_slope']   = (  0.247   if args.runPeriod == '2011' else  0.232,   0. )
pdfConfig['externalConstr']['sf_sigma_quad']    = (  3.32    if args.runPeriod == '2011' else  2.10,    0. )

if fixUpAcc :
    pdfConfig['externalConstr']['betaTimeEff'] = ( 0., 0. )

if trueTags :
    pdfConfig['externalConstr']['wTagP0OS']    = ( 0.3791, 0. )
    pdfConfig['externalConstr']['wTagP1OS']    = ( 1.,     0. )
    pdfConfig['externalConstr']['wTagP0SS']    = ( 0.445,  0. )
    pdfConfig['externalConstr']['wTagP1SS']    = ( 1.,     0. )
    pdfConfig['externalConstr']['wTagDelP0OS'] = ( 0.,     0. )
    pdfConfig['externalConstr']['wTagDelP1OS'] = ( 0.,     0. )
    pdfConfig['externalConstr']['wTagDelP0SS'] = ( 0.,     0. )
    pdfConfig['externalConstr']['wTagDelP1SS'] = ( 0.,     0. )
else :
    if fixTagging :
        pdfConfig['externalConstr']['wTagP0OS'] = ( 0.38152946, 0. )#( 0.3791, 0.00001 )
        pdfConfig['externalConstr']['wTagP1OS'] = ( 1.0118512,  0. )#( 1.,     0.00001 )
        pdfConfig['externalConstr']['wTagP0SS'] = ( 0.44585594, 0. )#( 0.445,  0.00001 )
        pdfConfig['externalConstr']['wTagP1SS'] = ( 0.95813206, 0. )#( 1.,     0.00001 )
    if fixTagAsym :
        pdfConfig['externalConstr']['wTagDelP0OS'] = (  0.014023729,  0. )#(  0.0140, 0.00001 )
        pdfConfig['externalConstr']['wTagDelP1OS'] = (  0.065743477,  0. )#(  0.066,  0.00001 )
        pdfConfig['externalConstr']['wTagDelP0SS'] = ( -0.015786075,  0. )#( -0.0158, 0.00001 )
        pdfConfig['externalConstr']['wTagDelP1SS'] = (  0.0081990069, 0. )#(  0.008,  0.00001 )

pdfConfig['anglesEffType'] = args.angAccType
pdfConfig['constAngEffCoefs'] = constAngAcc
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
    ws['lambdaCP'].setVal(1.)# 0.96340)
    ws['lambdaCP'].setConstant(True)
elif args.model == 'lamb_phi' :
    ws['lambdaCP'].setConstant(False)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(dataSet)
pdfPars = pdf.getParameters(dataSet)

# make S-wave parameters constant
for name in [ 'f_S', 'ASOddPhase' ] :
    par = pdfPars.find(name)
    if par.getMin() > 0. : par.setMin(0.)
    par.setVal(0.)
    par.setConstant(True)
pdfPars.find('C_SP').setConstant(True)

# prevent tagging dilution = 0 in initialization
if abs( 1. - pdfPars.getRealValue('wTagP1OS') ) < 1.e-5 and not pdfPars.find('wTagP1OS').isConstant() :
    pdfPars.setRealValue( 'wTagP1OS', 0.999 )
if abs( 1. - pdfPars.getRealValue('wTagP1SS') ) < 1.e-5 and not pdfPars.find('wTagP1SS').isConstant() :
    pdfPars.setRealValue( 'wTagP1SS', 0.999 )

# print parameters
print 120 * '='
print 'Bs2JpsiphiMCFit: data:'
dataSet.Print()
print 'Bs2JpsiphiMCFit: observables in PDF:'
pdfObs.Print('v')
print 'Bs2JpsiphiMCFit: parameters in PDF:'
pdfPars.Print('v')
print 'Bs2JpsiphiMCFit: constraints in PDF:'
for constr in pdf.ExternalConstraints() : constr.Print()


###########################################################################################################################################
## fit data ##
##############

# get parameters for Minos
if runMinos and not minosPars :
    minosPars += [ 'Gamma', 'dGamma', 'dM', 'AparPhase', 'AperpPhase' ]
    if args.model in [ 'phi', 'lamb_phi' ] :
        minosPars += [ 'phiCP', 'A0Mag2', 'AperpMag2' ]
    if args.model == 'lamb_phi' :
        minosPars += [ 'lambdaCP' ]
    elif args.model == 'polarDep' :
        minosPars += [ 'phiCPAv', 'phiCPRel_Apar', 'phiCPRel_AperpApar', 'phiCPRel_AS'
                      , 'CCPAv', 'CCPRel_Apar', 'CCPRel_Aperp', 'CCPAv_AS', 'avA02', 'avAperp2' ]
minosPars = [ pdfPars.find(par) for par in minosPars ]
if minosPars :
    print 'Bs2JpsiphiMCFit: running Minos for parameters', ', '.join( '"%s"' % par.GetName() for par in minosPars )

# fit data
if dataSet.isWeighted() :
    from P2VV.Utilities.DataHandling import correctWeights
    fitData = correctWeights( dataSet, [ ] )
else :
    fitData = dataSet
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = args.numCPU, Optimize = 2, Timer = True, Minimizer = 'Minuit2'
                      , Strategy = 1, Offset = True, Hesse = runHesse, Minos = minosPars )

if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    from P2VV.Imports import parNames, parValuesCPVDecay as parValues
elif args.model == 'phi' :
    from P2VV.Imports import parNames, parValuesFixLamb as parValues
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

if resultFileOut :
    # write fit result to file
    from ROOT import TFile, TObject
    resultFile = TFile.Open( resultFileOut, 'RECREATE' )
    resultFile.Append(fitResult)
    resultFile.Write( resultFileOut, TObject.kOverwrite )
    resultFile.Close()
