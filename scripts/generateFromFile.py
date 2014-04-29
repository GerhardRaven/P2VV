###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobID', '-I', type = int, default = 1000000 )
parser.add_argument( '--inputFile', '-f', default = 'cbkgInput.root' )
parser.add_argument( '--inputNTupleName', '-n', default = 'DecayTree' )
parser.add_argument( '--protoFile', '-p', default = 'genProtoData.root' )
parser.add_argument( '--protoDataName', '-m', default = 'JpsiKK_sigSWeight_cbkg' )
parser.add_argument( '--outputFile', '-o', default = 'genData_cbkg.root' )
args = parser.parse_args()

# print script settings
print 'job parameters:'
print '  job ID: %d' % args.jobID
print '  input file: %s' % args.inputFile
print '  input n-tuple name: %s' % args.inputNTupleName
print '  proto-data file: %s' % args.protoFile
print '  proto-dataset name: %s' % args.protoDataName
print '  output file: %s' % args.outputFile

obsNames = [ 'time', 'helcosthetaK', 'helcosthetaL', 'helphi' ]
splitCats = [ 'etaOS', 'etaSS', 'hlt1Cat', 'tResBin' ]
splitPars = dict(  etaOS   = ( 'cat',   'tagCatP2VVOS',      [ 0, 1 ], { 0 : 0, 1 : 1 } )
                 , etaSS   = ( 'cat',   'tagCatP2VVSS',      [ 0, 1 ], { 0 : 0, 1 : 1 } )
                 , hlt1Cat = ( 'cat',   'hlt1_excl_biased',  [ 0, 1 ], { 0 : 0, 1 : 1 } )
                 , tResBin = ( 'float', 'sigmat',            [ 0, 1, 2 ], [ ( 0.028, 0 ), ( 0.036, 1 ) ] )
                )
from itertools import product
states = [ stList for stList in product( *tuple( [ st for st in splitPars[cat][2] ] for cat in splitCats ) ) ]


###########################################################################################################################################
## select events from input file ##
###################################

# read prototype data
from ROOT import TFile, RooArgSet
protoFile = TFile.Open(args.protoFile)
protoData = protoFile.Get(args.protoDataName)
protoFile.Close()
assert protoData and all( protoData.get().find(obs) for obs in obsNames )
protoSet = RooArgSet( protoData.get() )
print 'P2VV - INFO: generateFromFile(): prototype data:'
protoData.Print()

# read input data
inputFile = TFile.Open(args.inputFile)
assert inputFile
inNTuples = dict( [ ( state, ( inputFile.Get( '%s_%s' % ( args.inputNTupleName
                              , '_'.join( '%s%d' % ( cat, st ) for cat, st in zip( splitCats, state ) ) ) ), [ ] ) )\
                    for state in states ] )
print 'P2VV - INFO: generateFromFile(): numbers of input events:'
for state in states :
    print '  %s: %7d' % ( state, inNTuples[state][0].GetEntries() )

# create output dataset
from ROOT import RooDataSet
outputData = RooDataSet( args.protoDataName, args.protoDataName, protoSet )

# draw events from input datasets and combine with proto data
import sys
import P2VV.RooFitWrappers
from ROOT import TRandom3
randGen = TRandom3(args.jobID)
for pEvIt, pSet in enumerate(protoData) :
    if pEvIt % 10000 == 0 :
        print 'P2VV - INFO: generateFromFile(): processing event %d' % pEvIt
        sys.stdout.flush()

    # get state of conditional observables in proto data
    state = [ ]
    for cat in splitCats :
        if splitPars[cat][0] == 'cat' :
            val = pSet.getCatIndex( splitPars[cat][1] )
            bin = splitPars[cat][3][val]
        else :
            val = pSet.getRealValue( splitPars[cat][1] )
            bin = splitPars[cat][2][-1]
            for thresh in splitPars[cat][3] :
                if val < thresh[0] :
                    bin = thresh[1]
                    break
        state.append(bin)

    # draw event from input n-tuple
    input = inNTuples[ tuple(state) ]
    while True :
        nEv = input[0].GetEntries()
        assert nEv - len(input[1]) > 0
        ev = int( float(nEv) * ( 1. - randGen.Rndm() ) )
        if ev in input[1] : continue
        input[1].append(ev)
        break
    input[0].GetEntry(ev)
    for obs in obsNames :
        protoSet.setRealValue( obs, getattr( input[0], obs ) )
    outputData.add(protoSet)

inputFile.Close()

# write generated dataset to file
print 'P2VV - INFO: generateFromFile(): writing generated data to %s' % args.outputFile
outputData.Print()
from ROOT import TObject
outputFile = TFile.Open( args.outputFile, 'RECREATE' )
outputFile.Add(outputData)
outputFile.Write( args.outputFile, TObject.kOverwrite )
