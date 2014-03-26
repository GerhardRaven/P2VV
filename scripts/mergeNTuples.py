nTupleFilePathsIn = ['/project/bfys/raaij/p2vv/data/2011_tupleB.root',
                     '/project/bfys/raaij/p2vv/data/2012_tupleB.root']
nTupleFilePathOut = '/project/bfys/raaij/p2vv/data/combined_tupleB.root'
nTupleName = 'DecayTree'

from ROOT import TFile, TList
nTupleFilesIn = [ ]
nTuplesIn = TList()
print 'reading n-tuples:'
for filePath in nTupleFilePathsIn :
    nTupleFilesIn.append( TFile.Open(filePath) )
    nTuplesIn.Add( nTupleFilesIn[-1].Get(nTupleName) )
    print '    %d entries in file "%s"' % ( nTuplesIn.Last().GetEntries(), filePath )

from ROOT import TTree
print 'merging n-tuples'
nTupleFileOut = TFile.Open( nTupleFilePathOut, 'RECREATE' )
nTupleOut = TTree.MergeTrees(nTuplesIn)
for nTupleFile in nTupleFilesIn :
    nTupleFile.Close()
    del nTupleFile

from ROOT import TObject
print 'writing merged n-tuple with %d entries to file "%s"' % ( nTupleOut.GetEntries(), nTupleFilePathOut )
nTupleFileOut.Write( nTupleFilePathOut, TObject.kOverwrite )
nTupleFileOut.Close()
