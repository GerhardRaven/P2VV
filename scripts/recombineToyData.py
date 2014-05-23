# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--inputID', '-i', type = int, default = 0 )
parser.add_argument( '--numFiles', '-n', type = int, default = 1 )
args = parser.parse_args()

print 'job parameters:'
print '  input ID: %d' % args.inputID
print '  number of files to process: %d' % args.numFiles

# other options
import P2VV.RooFitWrappers
inputFilePaths = [ '/data/bfys/jleerdam/Bs2Jpsiphi/combBkgGen/Bkg_sample_2014_100K_%d.root' % id\
                   for id in range( args.inputID * args.numFiles, ( args.inputID + 1 ) * args.numFiles ) ]
outputFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/Bs2JpsiKKToys/cbkgDataRecomb/cbkg_%03d.root' % args.inputID
nTupleNameIn = 'T'
nTupleNameOut = 'DecayTree'

obsNames = [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ]
splitCats = [ 'etaOS', 'etaSS', 'hlt1Cat', 'tResBin' ]
splitPars = dict(  etaOS   = ( 'cat',   'tagdecision_os_cb', [ 0, 1 ], { 0 : 0, 1 : 1, -1 : 1 } )
                 , etaSS   = ( 'cat',   'tagdecision_ss_nn', [ 0, 1 ], { 0 : 0, 1 : 1, -1 : 1 } )
                 , hlt1Cat = ( 'cat',   'hlt1_excl_biased',  [ 0, 1 ], { 0 : 0, 1 : 1 }         )
                 , tResBin = ( 'float', 'sigmat',            [ 0, 1, 2 ], [ ( 0.028, 0 ), ( 0.036, 1 ) ] )
                )
from itertools import product
states = [ stList for stList in product( *tuple( [ st for st in splitPars[cat][2] ] for cat in splitCats ) ) ]

from array import array
obsAdds = dict( [ ( name, array( 'd', [ 0. ] ) ) for name in obsNames ] )

# open output file and create output trees
from ROOT import TFile, TTree, TObject
outputFile = TFile.Open( outputFilePath, 'RECREATE' )
outputNTuples = { }
for state in states :
    nTName = '%s_%s' % ( nTupleNameOut, '_'.join( '%s%d' % ( name, st ) for name, st in zip( splitCats, state ) ) )
    outputNTuples[state] = TTree( nTName, nTName )
    for name in obsNames : outputNTuples[state].Branch( name, obsAdds[name], name + '/D' )

# loop over subset of input files
nFilesRead = 0
for path in inputFilePaths :
    inputFile = TFile.Open(path)
    if not inputFile : continue
    inputNTuple = inputFile.Get(nTupleNameIn)
    assert inputNTuple

    print 'P2VV - INFO: recombineToyData(): reading file %s' % inputFile.GetName()
    #for name in obsNames : inputNTuple.SetBranchAddress( name, obsAdds[name] )
    for ev in inputNTuple :
        state = [ ]
        for obs in obsNames :
            obsAdds[obs][0] = getattr( ev, obs )
        for cat in splitCats :
            val = getattr( ev, splitPars[cat][1] )
            if splitPars[cat][0] == 'cat' :
                bin = splitPars[cat][3][val]
            else :
                bin = splitPars[cat][2][-1]
                for thresh in splitPars[cat][3] :
                    if val < thresh[0] :
                        bin = thresh[1]
                        break
            state.append(bin)
        outputNTuples[ tuple(state) ].Fill()
    inputFile.Close()
    nFilesRead += 1
print 'P2VV - INFO: recombineToyData(): processed %d input files' % nFilesRead

# write output file
print 'P2VV - INFO: recombineToyData(): writing %d events to %s'\
      % ( sum( nt.GetEntries() for nt in outputNTuples.itervalues() ), outputFilePath )
outputFile.Write( outputFilePath, TObject.kOverwrite )
outputFile.Close()
