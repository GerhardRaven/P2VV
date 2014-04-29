# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--inputFiles', '-f', nargs = '+', default = [ 'genData_sig.root', 'genData_cbkg.root' ] )
parser.add_argument( '--inputNames', '-n', nargs = '+', default = [ 'JpsiKK_sigSWeight_sig', 'JpsiKK_sigSWeight_cbkg' ] )
parser.add_argument( '--outputFile', '-o', default = 'genData.root' )
parser.add_argument( '--outputName', '-m', default = 'JpsiKK_sigSWeight' )
parser.add_argument( '--weightName', '-w', default = 'sigWeight' )
args = parser.parse_args()
assert len(args.inputFiles) > 1 and len(args.inputNames) == len(args.inputFiles)

# print script settings
print 'job parameters:'
print '  input files: %s' % args.inputFiles
print '  input names: %s' % args.inputNames
print '  output file: %s' % args.outputFile
print '  output name: %s' % args.outputName
print '  event weight name: %s' % args.weightName

# get input datasets
import P2VV.RooFitWrappers
from ROOT import TFile
inData = [ ]
for fileName, dataName in zip( args.inputFiles, args.inputNames ) :
    inFile = TFile.Open(fileName)
    inData.append( inFile.Get(dataName) )
    inFile.Close()
    assert inData[-1]
    print 'P2VV - INFO: mergeDataSets: read dataset "%s" from file "%s"' % ( dataName, fileName )
    inData[-1].Print()

# merge files
from ROOT import RooDataSet, RooArgSet
obsSet = RooArgSet( inData[0].get() )
outData = RooDataSet( args.outputName, args.outputName, obsSet )
for data in inData : outData.append(data)
outData = RooDataSet( args.outputName, args.outputName, obsSet, Import = outData, WeightVar = ( args.weightName, True ) )

# write dataset to file
print 'P2VV - INFO: mergeDataSets: writing dataset to file "%s"' % args.outputFile
outData.Print()
from ROOT import TObject
outFile = TFile.Open( args.outputFile, 'RECREATE' )
outFile.Add(outData)
outFile.Write( args.outputFile, TObject.kOverwrite )
outFile.Close()
