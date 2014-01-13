nTupleFilePathIn  = '/data/bfys/jleerdam/Bs2Jpsiphi/Bs2JpsiPhi_2012_s20r0p1_dv33r6p1_20131217_tupleB.root'
nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2012_s20r0p1_dv33r6p1_20131217_tupleB_add.root'

nTupleName = 'DecayTree'
runPeriod = 2012
KKMassBounds = [ 1008., 1016., 1020., 1024., 1032. ]
KKMassInds = range( len(KKMassBounds) + 1 )
pbkgWeight = '' # 'wMC'

from ROOT import TFile
print 'reading file "%s" and cloning n-tuple "%s"' % ( nTupleFilePathIn, nTupleName )
nTupleFileIn = TFile.Open(nTupleFilePathIn)
nTupleIn = nTupleFileIn.Get(nTupleName)
nTupleFileOut = TFile.Open( nTupleFilePathOut, 'RECREATE' )
nTupleOut = nTupleIn.CloneTree()
nTupleFileIn.Close()
del nTupleFileIn

from P2VV.Load import P2VVLibrary
if runPeriod :
    from ROOT import addIntegerToTree
    print 'adding run period "%s" to n-tuple' % runPeriod
    addIntegerToTree( nTupleOut, runPeriod, 'runPeriod' )

if KKMassBounds :
    from ROOT import addCategoryToTree, std
    print 'adding KK-mass category with indices "%s" and KK-mass boundaries "%s" to n-tuple' % ( KKMassInds, KKMassBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in KKMassBounds : bounds.push_back(bound)
    for ind in KKMassInds : inds.push_back(ind)
    addCategoryToTree( nTupleOut, 'mdau2', 'KKMassCat', bounds, inds )

if pbkgWeight :
    from ROOT import copyFloatInTree
    print 'copying peaking-background weight "%s" to branch "pbkgWeight"' % pbkgWeight
    copyFloatInTree( nTupleOut, pbkgWeight, 'pbkgWeight' )

print 'first tree entries:'
for it in range( min( nTupleOut.GetEntries(), 20 ) ) :
    b = nTupleOut.GetEntry(it)
    print 'runPeriod = %4d   mdau2 = %6.1f   KKMassCat = %1d   pbkgWeight = %.2f'\
          % ( nTupleOut.runPeriod if runPeriod else -1, nTupleOut.mdau2, nTupleOut.KKMassCat if KKMassBounds else -1
             , nTupleOut.pbkgWeight if pbkgWeight else 0. )

from ROOT import TObject
print 'writing n-tuple to file "%s"' % nTupleFilePathOut
nTupleFileOut.Write( nTupleFilePathOut, TObject.kOverwrite )
nTupleFileOut.Close()
