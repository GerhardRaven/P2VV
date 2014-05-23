###########################################################################################################################################
## set script parameters ##
###########################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobID', '-I', type = int, default = 1000000 )
parser.add_argument( '--model', '-m', default = 'lamb_phi' )  # 'phi' / 'lamb_phi' / 'polarDep'
parser.add_argument( '--inputData', '-f', default = 'genProtoData.root' )
parser.add_argument( '--outputData', '-o' )
parser.add_argument( '--dataName', '-n', default = 'JpsiKK_sigSWeight_sig' )
parser.add_argument( '--parFileIn', '-i' )

args = parser.parse_args()
assert args.model in [ 'phi', 'lamb_phi', 'polarDep' ]
parFileIn = args.parFileIn
if parFileIn == None or parFileIn == 'None' :
    parFileIn = '20112012Reco14DataFitValues_6KKMassBins%s.par'\
                % ( '_CPVDecay' if args.model == 'polarDep' else '_fixedLamb' if args.model == 'phi' else '' )
outputDataFile = args.outputData
if not outputDataFile : outputDataFile = 'genData_sig_%s.root' % args.model

print 'job parameters:'
print '  job ID:', args.jobID
print '  model: %s' % args.model
print '  input data file: %s' % args.inputData
print '  output data file: %s' % outputDataFile
print '  dataset name: %s' % args.dataName
print '  input parameter file: %s' % parFileIn

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
plotsFilePath = '' #'genDists_sig_%s.pdf' % args.model

# PDF options
pdfConfig['lambdaCPParam'] = 'observables_CPVDecay' if args.model == 'polarDep' else 'lambPhi'
if pdfConfig['lambdaCPParam'] == 'observables_CPVDecay' :
    pdfConfig['splitParams']['KKMassCat'] = [ 'av' + par if par == 'f_S' else par for par in pdfConfig['splitParams']['KKMassCat'] ]

pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file'] = 'timeAcceptanceFit_2011.root'
pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file'] = 'timeAcceptanceFit_2012.root'
pdfConfig['anglesEffType'] = 'basisSig6'
pdfConfig['angEffMomsFiles'] = 'angEffNominalRew_5thOrder.par'


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
inputData = readData( filePath = args.inputData, dataSetName = args.dataName,  NTuple = False )
pdfConfig['signalData'] = inputData
pdfConfig['readFromWS'] = True

# build the PDF
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

if args.model == 'phi' :
    # set |lambda| = 1
    ws['lambdaCP'].setVal(1.)
    ws['lambdaCP'].setConstant(True)

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(inputData)
pdfPars = pdf.getParameters(inputData)

# print parameters
print 120 * '='
print 'genBs2JpsiKKSignal: observables in PDF:'
pdfObs.Print('v')
print 'genBs2JpsiKKSignal: parameters in PDF:'
pdfPars.Print('v')


###########################################################################################################################################
## generate data ##
###################

# set seed of RooFit random-number generator
from ROOT import RooRandom
RooRandom.randomGenerator().SetSeed(args.jobID)

# get observables to generate and create prototype dataset
genObsSet = [ ws[var] for var in [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ] ]
condObsSet = [ var for var in inputData.get() if not any( var.GetName() == genVar.GetName() for genVar in genObsSet ) ]
if condObsSet :
    protoData = inputData.reduce( ArgSet = condObsSet )
    #protoData = protoData.reduce( Cut = 'hlt1_excl_biased==1' )
    print 'genBs2JpsiKKSignal: proto-data set:'
    protoData.Print()
else :
    protoData = inputData

# set maximum value of PDF function
from P2VV.RooFitWrappers import SimultaneousPdf
if isinstance( pdf, SimultaneousPdf ) :
    for cat in pdf.indexCat() : pdf.getPdf( cat.GetName() ).setMaxVal(0.1)
else :
    pdf.setMaxVal(0.1)

# generate data
outputData = pdf.generate( genObsSet, ProtoData = ( protoData, False, False ), Name = args.dataName )
print 'genBs2JpsiKKSignal: saving generated data set to %s:' % outputDataFile
outputData.Print()
from ROOT import TFile, TObject
outFile = TFile.Open( outputDataFile, 'RECREATE' )
outFile.Add(outputData)
outFile.Write( outputDataFile, TObject.kOverwrite )

if plotsFilePath :
    # make some plots
    frames = [ ]
    canvs  = [ ]
    from P2VV.Load import LHCbStyle
    from ROOT import TCanvas
    for var in genObsSet :
        if var.GetName() == 'time' :
            frames.append( var.frame( 0.3, 5., 40 ) )
        else :
            frames.append( var.frame(30) )
        frames[-1].SetMarkerSize(0.1)
        outputData.plotOn( frames[-1], MarkerSize = 0.7 )
        frames[-1].SetMinimum(0.)
    
        canvs.append( TCanvas( var.GetName() + '_canv' ) )
        frames[-1].Draw()

    for it, canv in enumerate(canvs) :
        canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(canvs) - 1 else '' ) )
