nTupleFilePathsIn = [  '/data/bfys/jleerdam/Bs2Jpsiphi/2011_dv33r6p1_s20r1p1_20130919_tupleB_add.root'
                     , '/data/bfys/jleerdam/Bs2Jpsiphi/2012_dv33r6p1_s20r0p1_20130922_tupleB_add.root' ]
nTupleFilePathOut = '/data/bfys/jleerdam/Bs2Jpsiphi/2011_2012_dv33r6p1_s20_201309_tupleB_add.root'
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
