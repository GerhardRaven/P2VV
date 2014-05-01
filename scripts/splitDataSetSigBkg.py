###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobID', '-I', type = int, default = 1000000 )
parser.add_argument( '--generateMass', '-m', default = False )
parser.add_argument( '--weightedData', '-w', default = False )
parser.add_argument( '--dataPathIn', '-f', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140421.root' )
parser.add_argument( '--dataPathOut', '-o', default = 'genProtoData.root' )
parser.add_argument( '--datasetName', '-n', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--parFilePathIn', '-i', default = 'JpsiKKMassPars.par' )
args = parser.parse_args()
genMass = True if args.generateMass else False
weightData = True if args.weightedData else False

# print script settings
print 'job parameters:'
print '  job ID: %d' % args.jobID
print '  generate mass: %s' % ( 'true' if genMass else 'false' )
print '  weighted data: %s' % ( 'true' if weightData else 'false' )
print '  data file in: %s' % args.dataPathIn
print '  data path out: %s' % args.dataPathOut
print '  dataset name: %s' % args.datasetName
print '  parameter file in: %s' % args.parFilePathIn

constSplitPars = [ 'm_sig_alpha_1', 'm_sig_alpha_2', 'm_sig_n_1', 'm_sig_n_2', 'm_sig_lambda', 'm_sig_zeta' ]
condObs = dict(  runPeriod = [ 'runPeriod'            ]
               , HLT1Cat   = [ 'HLT1Cat', 'runPeriod' ]
               , KKMass    = [ 'KKMass'               ]
               , sigmat    = [ 'sigmat'               ]
               , etaTag    = [ 'etaTag'               ]
               , OSSSTagEq = [ 'OSSSTagEq', 'etaTag'  ]
              )
condObsVars = dict(  runPeriod = ( 'runPeriod', )
                   , HLT1Cat   = ( 'hlt1_excl_biased', )
                   , KKMass    = ( 'mdau2', )
                   , sigmat    = ( 'sigmat', )
                   , etaTag    = ( 'tagomega_os_cb', 'tagomega_ss_nn' )
                   , OSSSTagEq = ( 'iTagOS', 'iTagSS' )
                  )


###########################################################################################################################################
## build J/psiKK mass PDF ##
############################

# create workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiKKWS' ).ws()

# read data set
from P2VV.Utilities.DataHandling import readData
dataSet = readData( filePath = args.dataPathIn, dataSetName = args.datasetName,  NTuple = False )

# get observables
from P2VV.RooFitWrappers import RealVar, Category
observables = dict( mass = RealVar('mass'), runPeriod = Category('runPeriod'), KKMassCat = Category('KKMassCat') )

# initialize PDF components
from P2VV.RooFitWrappers import Component
nSig = dataSet.sumEntries()
nTot = float( dataSet.numEntries() )
sigMassComps  = Component( 'sigMass',  [ ], Yield = ( nSig,        0., nTot ) )
cbkgMassComps = Component( 'cbkgMass', [ ], Yield = ( nTot - nSig, 0., nTot ) )
massComps = [ sigMassComps, cbkgMassComps ]

# build the signal mass PDF
sigMassArgs = dict( Name = 'sig_m', mass = observables['mass'] )
from P2VV.Parameterizations.MassPDFs import Ipatia2_Signal_Mass as SignalBMass
sigMassArgs['m_sig_lambda']  = -2.5
sigMassArgs['m_sig_zeta']    =  0.01
sigMassArgs['m_sig_alpha_1'] =  3.0
sigMassArgs['m_sig_alpha_2'] =  2.5
sigMassArgs['m_sig_n_1']     =  1.5
sigMassArgs['m_sig_n_2']     =  2.0
signalBMass = SignalBMass( **sigMassArgs )
sigMassComps += signalBMass.pdf()

# build the combinatorial background mass PDF
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass
backgroundBMass = BackgroundBMass( Name = 'cbkg_m', mass = observables['mass'] )
cbkgMassComps += backgroundBMass.pdf()

# build total mass PDF
from P2VV.RooFitWrappers import buildPdf
massPdf = buildPdf( massComps, Observables = [ observables['mass'] ], Name = 'JpsiKKMass' )

# get splitting categories and parameters
yieldNames  = [ par.GetName() for par in massPdf.Parameters() if par.getAttribute('Yield') ]
allParNames = [ par.GetName() for par in massPdf.Parameters() if not par.getAttribute('Yield') if not par.isConstant() ]
sigParNames = [ par.GetName() for par in signalBMass.parameters() if not par.isConstant() ]
bkgParNames = [ par.GetName() for par in backgroundBMass.parameters() if not par.isConstant() ]
splitCats = {  observables['runPeriod'].GetName() : set( yieldNames + bkgParNames )
             , observables['KKMassCat'].GetName() : set( yieldNames               )
            }

# split Ipatia 2 signal parameters
KKCatName = observables['KKMassCat'].GetName()
for parName in sigParNames :
    splitCats[KKCatName].add(parName)
for parName in constSplitPars :
    splitCats[KKCatName].add(parName)

splitParsDict = { }
pdfVars = massPdf.getVariables()
for cat, params in splitCats.iteritems() :
    assert ws.cat(cat), 'P2VV - ERROR: splitDataSetSigBkg: category "%s" not in workspace' % cat
    for par in params :
        assert ws.var(par), 'P2VV - ERROR: splitDataSetSigBkg: no variable "%s" in workspace' % par
        assert pdfVars.find(par), 'P2VV - ERROR: splitDataSetSigBkg: variable "%s" not in PDF' % par
        if ws[par] not in splitParsDict :
            splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
        else :
            splitParsDict[ ws[par] ].add( ws[cat] )

# create list of split parameters and categories
from P2VV.Utilities.General import createSplitParsList
splitPars = createSplitParsList(splitParsDict)

# build simultaneous mass PDF
print 'P2VV - INFO: splitDataSetSigBkg: building simultaneous PDF "%s":' % ( massPdf.GetName() + '_simul' )
print 13 * ' ' + 'split parameters:'
for it, pars in enumerate(splitPars) :
    print 13 * ' ' + '%2d: pars: [ %s ]' % ( it, ', '.join( par.GetName() for par in pars[0] ) )
    print 13 * ' ' + '    cats: [ %s ]' % ', '.join( cat.GetName() for cat in pars[1] )
from P2VV.RooFitWrappers import SimultaneousPdf
simMassPdf = SimultaneousPdf(  massPdf.GetName() + '_simul'
                             , MasterPdf       = massPdf
                             , SplitCategories = [ pars[1] for pars in splitPars ]
                             , SplitParameters = [ pars[0] for pars in splitPars ]
                            )
massPdfObs  = simMassPdf.getObservables(dataSet)
massPdfPars = simMassPdf.getParameters(dataSet)

# create PDF configuration object
from P2VV.Parameterizations.FullPDFs import PdfConfiguration as PdfConfig
pdfConfig = PdfConfig()

# read parameters from file
print 120 * '='
pdfConfig.readParametersFromFile( filePath = args.parFilePathIn )
pdfConfig.setParametersInPdf(simMassPdf)
print 'P2VV - INFO: splitDataSetSigBkg: observables in mass PDF:'
massPdfObs.Print('v')
print 'P2VV - INFO: splitDataSetSigBkg: parameters in mass PDF:'
massPdfPars.Print('v')


###########################################################################################################################################
## split data into signal and background ##
###########################################

# create integration sets
from ROOT import RooArgSet
intSet = RooArgSet()
normSet = RooArgSet( massPdfObs.find('mass') )

# get signal and background mass PDFs
comps = [ 'sig', 'cbkg' ]
yields = dict( sig = { }, cbkg = { } )
pdfs   = dict( sig = { }, cbkg = { } )
from P2VV.Utilities.General import getSplitPar
splitCat = simMassPdf.indexCat()
for state in splitCat :
    stateName = state.GetName()
    yields['sig'][stateName]  = getSplitPar( 'N_sigMass',  stateName, massPdfPars ).getVal()
    yields['cbkg'][stateName] = getSplitPar( 'N_cbkgMass', stateName, massPdfPars ).getVal()
    if genMass :
        pdfs['sig'][stateName]  = simMassPdf.getPdf(stateName).pdfList().at(0)
        pdfs['cbkg'][stateName] = simMassPdf.getPdf(stateName).pdfList().at(1)
    else :
        pdfs['sig'][stateName]  = simMassPdf.getPdf(stateName).pdfList().at(0).createIntegral( intSet, normSet )
        pdfs['cbkg'][stateName] = simMassPdf.getPdf(stateName).pdfList().at(1).createIntegral( intSet, normSet )
for comp in comps :
    yields[comp]['total'] = sum( y for y in yields[comp].itervalues() )

print 'P2VV - INFO: signal and background yields:'
for comp in comps :
    print '  %4s: %9.2f' % ( comp, yields[comp]['total'] )

# create datasets
from ROOT import RooFit, RooDataSet, RooRealVar
dataObs = RooArgSet( dataSet.get() )
if genMass :
    var = dataObs.find('mass')
    if var : dataObs.remove(var)
    splitData = dict( [ ( comp, dict( [ ( state.GetName(), RooDataSet('%s_%s' % (dataSet.GetName(), comp), dataSet.GetTitle(), dataObs) )\
                                        for state in splitCat ] ) ) for comp in comps ] )
else :
    sigWeight = RooRealVar( 'sigWeight', 'sigWeight', 1. )
    dataObs.add(sigWeight)
    if weightData :
        splitData = dict( [ ( comp, dict( [ ( 'total', RooDataSet( '%s_%s' % ( dataSet.GetName(), comp ), dataSet.GetTitle()
                                                                  , dataObs, RooFit.WeightVar(sigWeight) ) ) ] ) ) for comp in comps ] )
    else :
        splitData = dict( [ ( comp, dict( [ ( 'total', RooDataSet( '%s_%s' % ( dataSet.GetName(), comp ), dataSet.GetTitle()
                                                                  , dataObs ) ) ] ) ) for comp in comps ] )

# function to get observable value from dataset
def getObsVal( names, obsSet ) :
    vals = [ ]
    for name in names :
        catInd = obsSet.getCatIndex( name, -999 )
        vals.append( catInd if catInd != -999 else obsSet.getRealValue(name) )
    return tuple(vals) if len(vals) > 1 else vals[0]

# fill datasets
from P2VV.Utilities.Studies import condObsDist
from ROOT import TRandom3
randGen = TRandom3(args.jobID)
for evObs in dataSet :
    # get values of observables
    massPdfObs.assignValueOnly(evObs)
    dataObs.assignValueOnly(evObs)
    stateName = splitCat.getLabel()
    if not genMass : sigWeight.setVal( dataSet.weight() )

    # get PDF values
    pdfVals = dict( [ ( comp, yields[comp]['total'] * ( 1. if genMass else pdfs[comp][stateName].getVal() ) ) for comp in comps ] )
    for varList in condObs.itervalues() :
        valList = [ getObsVal( condObsVars[var], evObs ) for var in varList ]
        for comp in comps :
            pdfVals[comp] *= condObsDist( varList, valList, comp )

    # decide if event is signal or background
    totVal = sum( pdfVals.itervalues() )
    addVal = 0.
    genVal = randGen.Rndm()
    for comp in comps :
        addVal += pdfVals[comp]
        if addVal < genVal * totVal : continue
        if genMass :
            splitData[comp][stateName].add(dataObs)
        elif weightData :
            splitData[comp]['total'].add( dataObs, sigWeight.getVal() )
        else :
            splitData[comp]['total'].add(dataObs)
        break

if genMass :
    # generate J/psiKK mass
    genObs = RooArgSet( massPdfObs.find('mass') )
    dataObs.add(genObs)
    genData = RooDataSet( dataSet.GetName(), dataSet.GetTitle(), dataObs, WeightVar = 'wMC' )
    for comp in comps :
        splitData[comp]['total'] = RooDataSet( '%s_%s' % ( dataSet.GetName(), comp ), dataSet.GetTitle(), dataObs )
        for state in splitCat :
            stateName = state.GetName()
            splitData[comp][stateName] = pdfs[comp][stateName].generate( genObs, ProtoData = ( splitData[comp][stateName], False, False ) )
            splitData[comp]['total'].append( splitData[comp][stateName] )
        splitData[comp]['total'] = RooDataSet( '%s_%s' % ( dataSet.GetName(), comp ), dataSet.GetTitle(), dataObs
                                              , Import = splitData[comp]['total'], WeightVar = ( 'wMC', True ) )
        genData.append( splitData[comp]['total'] )

# print datasets and write them to file
from ROOT import TFile, TObject
dataFile = TFile.Open( args.dataPathOut, 'RECREATE' )
print 'P2VV - INFO: produced datasets (%d entries, sum of weights = %.1f):'\
      % (  sum( ds['total'].numEntries() for ds in splitData.itervalues() )
         , sum( ds['total'].sumEntries() for ds in splitData.itervalues() )
        )
if genMass :
    genData.Print()
    dataFile.Add(genData)
for comp in comps :
    splitData[comp]['total'].Print()
    dataFile.Add( splitData[comp]['total'] )
    if genMass :
        for state in splitCat :
            ds = splitData[comp][ state.GetName() ]
            print '%s %s: %d entries (%.1f weighted)' % ( comp, state.GetName(), ds.numEntries(), ds.sumEntries() )

dataFile.Write( args.dataPathOut, TObject.kOverwrite )
dataFile.Close()
