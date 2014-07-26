###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobName', '-N', default = 'Bs2JpsiKKFit' )
parser.add_argument( '--model', '-m', default = 'lamb_phi' )  # 'phi' / 'lamb_phi' / 'polarDep'
parser.add_argument( '--blind', '-b', default = True )
parser.add_argument( '--tagCatsType', '-T', default = 'linearCats' )
parser.add_argument( '--contEstWTag', '-C', default = True )
parser.add_argument( '--fixLowAcc', '-l', default = True )
parser.add_argument( '--fixUpAcc', '-u', default = False )
parser.add_argument( '--fixTagging', '-a', default = False )
parser.add_argument( '--fixTagAsym', '-k', default = False )
parser.add_argument( '--constrDelM', '-M', default = False )
parser.add_argument( '--numCPU', '-c', type = int, default = 2 )
parser.add_argument( '--runHesse', '-e', default = True )
parser.add_argument( '--runMinos', '-s', default = '' )
parser.add_argument( '--minosPars', '-p', default = None )
parser.add_argument( '--dataPath', '-d', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14' )
parser.add_argument( '--workPath', '-w', default = '/project/bfys/jleerdam/softDevel/P2VV2/test' )
parser.add_argument( '--dataSetFile', '-f', default = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root' )
parser.add_argument( '--accDataSetFile', '-g', default =  'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root' )
parser.add_argument( '--timeAccConstr', '-j', default = 'poisson' )
parser.add_argument( '--dataSetName', '-n', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--parFileIn', '-i' )
parser.add_argument( '--parFileOut', '-o' )
parser.add_argument( '--resultFileOut', '-r' )
parser.add_argument( '--timeAccFile2011', '-x', default = 'timeAcceptanceFit_2011.root' )
parser.add_argument( '--timeAccFile2012', '-y', default = 'timeAcceptanceFit_2012.root' )
parser.add_argument( '--angAccType', '-t', default = 'weights' )  # 'weights' / 'basis012' / 'basis01234' / 'basisSig6'
parser.add_argument( '--angAccFile', '-z', default = 'angEffNominalRew_moms.par' )
parser.add_argument( '--constAngAcc', '-q', default = True )
parser.add_argument( '--timeResSyst', '-R', default = '' )

args = parser.parse_args()
assert args.model in [ 'phi', 'lamb_phi', 'polarDep' ]
blind = False if not args.blind or str( args.blind ).lower() in [ 'false', '0' ] else True
assert args.tagCatsType in [ '', 'linearCats' ]
contEstWTag = False if not args.tagCatsType or str( args.contEstWTag ).lower() in [ 'false', '0' ] else True
fixLowAcc = False if not args.fixLowAcc or str( args.fixLowAcc ).lower() in [ 'false', '0' ] else True
fixUpAcc = False if not args.fixUpAcc or str( args.fixUpAcc ).lower() in [ 'false', '0' ] else True
fixTagging = False if not args.fixTagging or str( args.fixTagging ).lower() in [ 'false', '0' ] else True
fixTagAsym = False if not args.fixTagAsym or str( args.fixTagAsym ).lower() in [ 'false', '0' ] else True
constrDelM = False if not args.constrDelM or str( args.constrDelM ).lower() in [ 'false', '0' ] else True
assert type(args.numCPU) == int and args.numCPU > 0 and args.numCPU < 20
runHesse = False if not args.runHesse or str( args.runHesse ).lower() in [ 'false', '0' ] else True
minosPars = args.minosPars.split(',') if args.minosPars and str(args.minosPars) != 'None' else [ ]
runMinos = False if ( args.runMinos == '' and not minosPars ) or str( args.runMinos ).lower() in [ 'false', '0' ] else True
dataPath = args.dataPath
workPath = args.workPath
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
dataSetFile = dataPath + args.dataSetFile
accDataSetFile = dataPath + args.accDataSetFile
assert args.timeAccConstr in [ 'poisson', 'poisson_minimal', 'multinomial', 'average' ]
parFileIn = args.parFileIn
if parFileIn == None or parFileIn == 'None' :
    parFileIn = workPath + '20112012Reco14DataFitValues_6KKMassBins%s.par'\
                           % ( '_CPVDecay' if args.model == 'polarDep' else '_fixedLamb' if args.model == 'phi' else '' )
elif parFileIn :
    parFileIn = workPath + parFileIn
parFileOut = args.parFileOut
if parFileOut == None or parFileOut == 'None' :
    parFileOut = workPath + '%s%s_%sLow_%sUp_%sTag.par'\
                 % ( args.jobName + ( '_' if args.jobName else '' ), args.model, 'fix' if fixLowAcc else 'float'
                    , 'fix' if fixUpAcc else 'float', 'fix' if fixTagging else 'float' )
elif parFileOut :
    parFileOut = workPath + parFileOut
resultFileOut = args.resultFileOut
if resultFileOut == None or resultFileOut == 'None' :
    resultFileOut = workPath + '%s%s_%sLow_%sUp_%sTag.root'\
                 % ( args.jobName + ( '_' if args.jobName else '' ), args.model, 'fix' if fixLowAcc else 'float'
                    , 'fix' if fixUpAcc else 'float', 'fix' if fixTagging else 'float' )
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
print '  blind analysis: %s' % ( 'true' if blind else 'false' )
print '  tagging categories type: %s' % args.tagCatsType
print '  continuous estimated wrong-tag variable: %s' % contEstWTag
print '  fix lower decay-time acceptance: %s' % ( 'true' if fixLowAcc else 'false' )
print '  fix upper decay-time acceptance: %s' % ( 'true' if fixUpAcc  else 'false' )
print '  fix tagging calibration: %s' % ( 'true' if fixTagging else 'false' )
print '  fix tagging calibration asymmetries: %s' % ( 'true' if fixTagAsym else 'false' )
print '  constrain DeltaM: %s' % ( 'true' if constrDelM else 'false' )
print '  number of cores: %d' % args.numCPU
print '  run Hesse: %s' % ( 'true' if runHesse else 'false' )
print '  run Minos: %s' % ( 'true' if runMinos else 'false' )
print '  Minos parameters:', ', '.join( '"%s"' % par for par in minosPars )
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  dataset file: %s' % dataSetFile
print '  acceptance dataset file: %s' % accDataSetFile
print '  time acceptance constraint: %s' % args.timeAccConstr
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
pdfConfig = PdfConfig()

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )

if not blind :
    pdfConfig['blind'] = { }

#pdfConfig['KKMassBinBounds'] = [ 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12. ]
#pdfConfig['CSPValues'] = [ 0.9022, 0.8619, 0.8875, 0.9360 ]
#pdfConfig['KKMassBinBounds'] = [ 1020. - 12., 1020. + 12. ]
#pdfConfig['CSPValues'] = [ 0.498 ]

pdfConfig['tagCatsType'] = args.tagCatsType
pdfConfig['contEstWTag'] = contEstWTag
if pdfConfig['tagCatsType'] != 'linearCats' :
    for name in [ 'wTagP0OS', 'wTagP1OS', 'wTagDelP0OS', 'wTagDelP1OS', 'wTagP0SS', 'wTagP1SS', 'wTagDelP0SS', 'wTagDelP1SS' ] :
        pdfConfig['externalConstr'].pop(name)

pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if args.model == 'polarDep' else 'lambPhi'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]

pdfConfig['timeEffType'] = 'paper2012' if fixLowAcc else 'fit_uniformUB'
#timeAccHistHLT1UB2011 = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
timeAccHistHLT1UB2011 = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = timeAccFile2011
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = timeAccFile2012
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['hlt1UB'] = timeAccHistHLT1UB2011
if pdfConfig['timeEffType'].startswith('fit') :
    pdfConfig['timeEffConstraintType'] = args.timeAccConstr
    from P2VV.Parameterizations.FullPDFs import SimulCatSettings
    pdfConfig['timeEffData']['file'] = accDataSetFile
    pdfConfig['externalConstr']['acceptance'] = SimulCatSettings('acceptanceConstr')
    pdfConfig['externalConstr']['acceptance'].addSettings(  [ 'runPeriod' ], [ [ 'p2011', 'p2012' ] ]
                                                          , {  ( 'hlt1_excl_biased_dec', 'exclB' ) : ( 0.65, 0.01 )
                                                             , ( 'hlt2_biased', 'B' )              : ( 0.65, 0.01 )
                                                            }
                                                         )

timeResSystType = args.timeResSyst
from collections import defaultdict
splitConstr = defaultdict(dict)
if timeResSystType == 'mean_param':
    pdfConfig['timeResType'] += '_mean_param'
    splitConstr['mu_placeholder']['2011']   = (  0.0350, 0. )
    splitConstr['timeResMu_offset']['2011'] = ( -0.0020, 0. )
    splitConstr['timeResMu_slope']['2011']  = ( -0.1475, 0. )
    splitConstr['timeResMu_quad']['2011']   = ( -6.2781, 0. )

    splitConstr['mu_placeholder']['2012']  =  (  0.0349, 0. )
    splitConstr['timeResMu_offset']['2012'] = ( -0.0031, 0. )
    splitConstr['timeResMu_slope']['2012']  = ( -0.1634, 0. )
    splitConstr['timeResMu_quad']['2012']   = ( -6.9108, 0. )
    
    split_runPeriod = set(pdfConfig['splitParams']['runPeriod'])
    for pn in ('timeResMu',):
        pdfConfig['externalConstr'].pop(pn)
        pdfConfig['splitParams']['runPeriod'].remove(pn)
    for par, vals in splitConstr.iteritems() :
        if par not in split_runPeriod:
            pdfConfig['splitParams']['runPeriod'].append(par)
elif timeResSystType == 'quadratic_no_offset':
    from collections import defaultdict
    splitConstr = defaultdict(dict)
    splitConstr['sf_placeholder']['2011']  = (  0.0350,  0. )
    splitConstr['sf_mean_slope']['2011']   = (  1.436, 0. )
    splitConstr['sf_mean_quad']['2011']    = ( -8.145, 0. )
    splitConstr['timeResFrac2']['2011']    = (  0.151, 0. )
    splitConstr['sf_sigma_slope']['2011']  = (  0.348, 0. )
    splitConstr['sf_sigma_quad']['2011']   = (  2.518, 0. )

    splitConstr['sf_placeholder']['2012']  =  ( 0.0349,  0. )
    splitConstr['sf_mean_slope']['2012']   =  ( 1.425, 0. )
    splitConstr['sf_mean_quad']['2012']    =  ( 3.167, 0. )
    splitConstr['timeResFrac2']['2012']    =  ( 0.265, 0. )
    splitConstr['sf_sigma_slope']['2012']  =  ( 0.271, 0. )
    splitConstr['sf_sigma_quad']['2012']   =  ( 3.136, 0. )

for par, vals in splitConstr.iteritems() :
    from P2VV.Parameterizations.FullPDFs import SimulCatSettings
    constr = SimulCatSettings( '%sConstr' % par )
    constr.addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ], vals['2011'] )
    constr.addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ], vals['2012'] )
    pdfConfig['externalConstr'][par] = constr

    
if fixUpAcc :
    for it, sett in enumerate( pdfConfig['externalConstr']['betaTimeEff'] ) :
        per = sett[0]['runPeriod']
        assert per in [ [ 'p2011' ], [ 'p2012' ] ]
        betaVal = pdfConfig.parameters().pop( 'betaTimeEff_%s' % per[0] )[0]
        pdfConfig['externalConstr']['betaTimeEff'][it] = ( sett[0], ( betaVal, 0. ) )

if fixTagging :
    pdfConfig['externalConstr']['wTagP0OS'] = ( pdfConfig.parameters().pop('wTagP0OS')[0], 0. )
    pdfConfig['externalConstr']['wTagP1OS'] = ( pdfConfig.parameters().pop('wTagP1OS')[0], 0. )
    pdfConfig['externalConstr']['wTagP0SS'] = ( pdfConfig.parameters().pop('wTagP0SS')[0], 0. )
    pdfConfig['externalConstr']['wTagP1SS'] = ( pdfConfig.parameters().pop('wTagP1SS')[0], 0. )
if fixTagAsym :
    pdfConfig['externalConstr']['wTagDelP0OS'] = ( pdfConfig.parameters().pop('wTagDelP0OS')[0], 0. )
    pdfConfig['externalConstr']['wTagDelP1OS'] = ( pdfConfig.parameters().pop('wTagDelP1OS')[0], 0. )
    pdfConfig['externalConstr']['wTagDelP0SS'] = ( pdfConfig.parameters().pop('wTagDelP0SS')[0], 0. )
    pdfConfig['externalConstr']['wTagDelP1SS'] = ( pdfConfig.parameters().pop('wTagDelP1SS')[0], 0. )

if constrDelM :
    pdfConfig['externalConstr']['dM'] = ( 17.768, 0.024 )

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
    # set parameters in PDF
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

# prevent tagging dilution = 0 in initialization
if abs( 1. - pdfPars.getRealValue('wTagP1OS') ) < 1.e-5 and not pdfPars.find('wTagP1OS').isConstant() :
    pdfPars.setRealValue( 'wTagP1OS', 0.999 )
if abs( 1. - pdfPars.getRealValue('wTagP1SS') ) < 1.e-5 and not pdfPars.find('wTagP1SS').isConstant() :
    pdfPars.setRealValue( 'wTagP1SS', 0.999 )

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

# get parameters for Minos
if runMinos and not minosPars :
    blindAdd = '__' if pdfConfig['blind'] else ''
    minosPars += [ 'Gamma', '{0}dGamma{0}'.format(blindAdd), 'dM', 'AparPhase', 'AperpPhase' ]
    if args.model in [ 'phi', 'lamb_phi' ] :
        minosPars += [ '{0}phiCP{0}'.format(blindAdd), 'A0Mag2', 'AperpMag2' ]
    if args.model == 'lamb_phi' :
        minosPars += [ 'lambdaCP' ]
    elif args.model == 'polarDep' :
        minosPars += [ '{0}phiCPAv{0}'.format(blindAdd), '{0}phiCPRel_Apar{0}'.format(blindAdd)
                      , '{0}phiCPRel_AperpApar{0}'.format(blindAdd), '{0}phiCPRel_AS{0}'.format(blindAdd)
                      , 'CCPAv', 'CCPRel_Apar', 'CCPRel_Aperp', 'CCPAv_AS', 'avA02', 'avAperp2' ]
minosPars = [ pdfPars.find(par) for par in minosPars ]
if minosPars :
    print 'Bs2JpsiKK3fbFit: running Minos for parameters', ', '.join( '"%s"' % par.GetName() for par in minosPars )

# fit data
from P2VV.Utilities.DataHandling import correctWeights
fitData = correctWeights( dataSet, [ 'runPeriod', 'KKMassCat' ] )
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
