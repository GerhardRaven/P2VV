###########################################################################################################################################
## set script parameters ##
###########################

model = 'lamb_phi'  # 'phi' / 'lamb_phi' / 'polarDep'
fixLowAcc = True
fixUpAcc  = False
numCPU = 7
dataPath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
workPath = '/project/bfys/jleerdam/softDevel/P2VV2/test/'
dataSetFile = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root'

import sys
if len(sys.argv) > 1 :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('model')
    parser.add_argument('fixLowAcc')
    parser.add_argument('fixUpAcc')
    parser.add_argument( 'numCPU', type = int )
    if len(sys.argv) > 5 :
        parser.add_argument('dataPath')
        parser.add_argument('workPath')
        parser.add_argument('dataSetFile')

    args = parser.parse_args()
    model = args.model
    fixLowAcc = True if args.fixLowAcc.lower() == 'true' else False
    fixUpAcc  = True if args.fixUpAcc.lower()  == 'true' else False
    numCPU = args.numCPU
    if len(sys.argv) > 5 :
        if args.dataPath : dataPath = args.dataPath
        if args.workPath : workPath = args.workPath
        if args.dataSetFile : dataSetFile = args.dataSetFile

assert model in [ 'phi', 'lamb_phi', 'polarDep' ]
assert type(fixLowAcc) == bool
assert type(fixUpAcc)  == bool
assert type(numCPU) == int and numCPU > 0 and numCPU < 20
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
assert dataSetFile

print 'job parameters:'
print '  model: %s' % model
print '  fix lower decay-time acceptance: %s' % ( 'true' if fixLowAcc else 'false' )
print '  fix upper decay-time acceptance: %s' % ( 'true' if fixUpAcc  else 'false' )
print '  number of cores: %d' % numCPU
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  dataset file: %s' % dataSetFile

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
parFileIn  = workPath + '20112012Reco14DataFitValues_6KKMassBins%s.par' % ( '_CPVDecay' if model == 'polarDep' else '' )
parFileOut = workPath + '%s_%sLow_%sUp.par' % ( model, 'fix' if fixLowAcc else 'float', 'fix' if fixUpAcc else 'float' )
dataSetName = 'JpsiKK_sigSWeight'

# PDF options
pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if model == 'polarDep' else 'lambPhi'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]

pdfConfig['timeEffType'] = 'paper2012' if fixLowAcc else 'fit_uniformUB'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
        = dataPath + 'timeAcceptanceFit_2011.root' # 'Bs_HltPropertimeAcceptance_Data_2011_40bins_TOS.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
        = dataPath + 'timeAcceptanceFit_2012.root' # 'Bs_HltPropertimeAcceptance_Data_2012_40bins_TOS.root'
if pdfConfig['timeEffType'].startswith('fit') :
    from P2VV.Parameterizations.FullPDFs import SimulCatSettings
    pdfConfig['timeEffData']['file'] = dataPath + 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root'
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
pdfConfig['angEffMomsFiles'] = dataPath + 'angEffNominalRew_moms.par'


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = dataPath + dataSetFile, dataSetName = dataSetName,  NTuple = False )
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
if model == 'phi' :
    ws['lambdaCP'].setVal(1.)
    ws['lambdaCP'].setConstant(True)
elif model == 'lamb_phi' :
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
fitResult = pdf.fitTo( fitData, SumW2Error = False, Save = True, NumCPU = numCPU, Optimize = 2, Timer = True, Minimizer = 'Minuit2'
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
