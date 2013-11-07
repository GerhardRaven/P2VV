nTupleFilePathIn  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/2011_dv33r6p1_s20r1p1_20131030_tupleB.root'
nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/2011_dv33r6p1_s20r1p1_20131030_tupleB_add.root'
#nTupleFilePathIn  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2012s20r0p1_dv33r6p1_20131030_tupleB.root'
#nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2012s20r0p1_dv33r6p1_20131030_tupleB_add.root'
nTupleName = 'DecayTree'
runPeriod = 2011

from ROOT import TFile
print 'reading file "%s" and cloning n-tuple "%s"' % ( nTupleFilePathIn, nTupleName )
nTupleFileIn = TFile.Open(nTupleFilePathIn)
nTupleIn = nTupleFileIn.Get(nTupleName)
nTupleFileOut = TFile.Open( nTupleFilePathOut, 'RECREATE' )
nTupleOut = nTupleIn.CloneTree()
nTupleFileIn.Close()
del nTupleFileIn

from P2VV.Load import P2VVLibrary
from ROOT import addRunPeriodToTree
print 'adding run period "%s" to n-tuple' % runPeriod
addRunPeriodToTree( nTupleOut, runPeriod, 'runPeriod' )

if nTupleOut.GetEntries() > 9 :
    print 'first ten run-period entries:'
    for it in range(10) :
        b = nTupleOut.GetEntry(it)
        print nTupleOut.runPeriod

from ROOT import TObject
print 'writing n-tuple to file "%s"' % nTupleFilePathOut
nTupleFileOut.Write( nTupleFilePathOut, TObject.kOverwrite )
nTupleFileOut.Close()