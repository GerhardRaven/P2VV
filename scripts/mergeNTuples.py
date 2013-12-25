nTupleFilePathsIn = [  '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2011_s20r1p1_dv33r6p1_20131217_tupleB_add.root'
                     , '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2012_s20r0p1_dv33r6p1_20131217_tupleB_add.root' ]
nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2011_2012_s20_dv33r6p1_20131217_tupleB_add.root'
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
