#nTupleFilePathIn  = '/data/bfys/jleerdam/Bs2Jpsiphi/Bs2JpsiPhi_2012_s20r0p1_dv33r6p1_20131217_tupleB.root'
nTupleFilePathIn  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged_presc_20140822.root'
nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged_add_20140822.root'

nTupleName = 'DecayTree'
runPeriod = None # 2012
firstData = None # 94387
prescaleBounds = [ 87218., 94387., 104415., 118881., 133786. ]
prescaleInds = [ -1, 0, 1, 0, 1, -2 ]
KKMassBounds = [ 1008., 1016., 1020., 1024., 1032. ]
KKMassInds = range( len(KKMassBounds) + 1 )
BpTBounds = [ 4600., 6300., 8100., 10800. ]
BpTInds = range( len(BpTBounds) + 1 )
nPVBounds = [ 0.5 + float(it) for it in range(11) ]
nPVInds = range( len(nPVBounds) + 1 )
sWeight = 'sWeights_ipatia'
pbkgWeight = 'wMC'
tagOSBounds = [ ] # [ -531.1, -530.9, 530.9, 531.1 ]
tagSSBounds = [ ] # [ -531.1, -530.9, 530.9, 531.1 ]

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

if firstData :
    from ROOT import addCategoryToTree, std
    print 'adding "first data" flag to n-tuple for data with run number < %d' % firstData
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    bounds.push_back( float(firstData) )
    inds.push_back(1)
    inds.push_back(0)
    addCategoryToTree( nTupleOut, 'runNumber', 'firstData', bounds, inds )

if prescaleBounds :
    from ROOT import addCategoryToTree, std
    print 'adding prescale category with indices "%s" and run-number boundaries "%s" to n-tuple' % ( prescaleInds, prescaleBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in prescaleBounds : bounds.push_back(bound)
    for ind in prescaleInds : inds.push_back(ind)
    addCategoryToTree( nTupleOut, 'runNumber', 'hlt2_prescale', bounds, inds )

if KKMassBounds :
    from ROOT import addCategoryToTree, std
    print 'adding KK-mass category with indices "%s" and KK-mass boundaries "%s" to n-tuple' % ( KKMassInds, KKMassBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in KKMassBounds : bounds.push_back(bound)
    for ind in KKMassInds : inds.push_back(ind)
    addCategoryToTree( nTupleOut, 'mdau2', 'KKMassCat', bounds, inds )

if BpTBounds :
    from ROOT import addCategoryToTree, std
    print 'adding B-p_T category with indices "%s" and p_T boundaries "%s" to n-tuple' % ( BpTInds, BpTBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in BpTBounds : bounds.push_back(bound)
    for ind in BpTInds : inds.push_back(ind)
    addCategoryToTree( nTupleOut, 'B_Pt', 'BpTCat', bounds, inds )

if nPVBounds :
    from ROOT import addCategoryToTree, std
    print 'adding number of PVs category with indices "%s" and boundaries "%s" to n-tuple' % ( nPVInds, nPVBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in nPVBounds : bounds.push_back(bound)
    for ind in nPVInds : inds.push_back(ind)
    addCategoryToTree( nTupleOut, 'nPV', 'nPVCat', bounds, inds )

if sWeight :
    from ROOT import copyFloatInTree
    print 'copying signal sWeight "%s" to branch "sWeight_orig"' % sWeight
    copyFloatInTree( nTupleOut, sWeight, 'sWeight_orig' )

if pbkgWeight :
    from ROOT import copyFloatInTree
    print 'copying peaking-background weight "%s" to branch "pbkgWeight"' % pbkgWeight
    copyFloatInTree( nTupleOut, pbkgWeight, 'pbkgWeight' )

if tagOSBounds :
    from ROOT import addFloatToTree, addCategoryToTree, std
    tagOSInds = [ -999, -1, 0, +1, +999 ]
    catOSInds = [ -999, 1, -999, 1, -999 ]
    print 'adding OS tagging category with indices "%s" and true-ID boundaries "%s" to n-tuple' % ( tagOSInds, tagOSBounds )
    bounds  = std.vector('Double_t')()
    tagInds = std.vector('Int_t')()
    catInds = std.vector('Int_t')()
    for bound in tagOSBounds : bounds.push_back(bound)
    for ind in tagOSInds : tagInds.push_back(ind)
    for ind in catOSInds : catInds.push_back(ind)
    addCategoryToTree( nTupleOut, 'trueid', 'iTagOS', bounds, tagInds )
    addCategoryToTree( nTupleOut, 'trueid', 'tagCatP2VVOS', bounds, catInds )
    addFloatToTree( nTupleOut, 0., 'tagomega_os_cb' )

if tagSSBounds :
    from ROOT import addCategoryToTree, std
    tagSSInds = [ -999, -1, 0, +1, +999 ]
    catSSInds = [ -999, 1, -999, 1, -999 ]
    print 'adding SS tagging category with indices "%s" and true-ID boundaries "%s" to n-tuple' % ( tagSSInds, tagSSBounds )
    bounds  = std.vector('Double_t')()
    tagInds = std.vector('Int_t')()
    catInds = std.vector('Int_t')()
    for bound in tagSSBounds : bounds.push_back(bound)
    for ind in tagSSInds : tagInds.push_back(ind)
    for ind in catSSInds : catInds.push_back(ind)
    addCategoryToTree( nTupleOut, 'trueid', 'iTagSS', bounds, tagInds )
    addCategoryToTree( nTupleOut, 'trueid', 'tagCatP2VVSS', bounds, catInds )
    addFloatToTree( nTupleOut, 0., 'tagomega_ss_nn' )

print 'first tree entries:'
for it in range( min( nTupleOut.GetEntries(), 20 ) ) :
    b = nTupleOut.GetEntry(it)
    print 'runPeriod = %4d   firstData = %d   hlt2_prescale = %d   mdau2 = %6.1f   KKMassCat = %1d   BpTCat = %1d   nPVCat = %1d   sWeight = %.2f   pbkgWeight = %.2f'\
          % ( nTupleOut.runPeriod if runPeriod else -1, nTupleOut.firstData if firstData else -1
             , nTupleOut.hlt2_prescale if prescaleBounds else -1, nTupleOut.mdau2, nTupleOut.KKMassCat if KKMassBounds else -1
             , nTupleOut.BpTCat if BpTBounds else -1, nTupleOut.nPVCat if nPVBounds else -1
             , nTupleOut.sWeight_orig if sWeight else 0., nTupleOut.pbkgWeight if pbkgWeight else 0. )
    print '   trueid = %d   iTagOS = %d   iTagSS = %d   tagCatOS = %d   tagCatSS = %d   etaOS = %.1f   etaSS = %.1f'\
          % ( nTupleOut.trueid if tagOSBounds or tagSSBounds else 0
             , nTupleOut.iTagOS if tagOSBounds else 0, nTupleOut.iTagSS if tagSSBounds else 0
             , nTupleOut.tagCatP2VVOS if tagOSBounds else -1, nTupleOut.tagCatP2VVSS if tagSSBounds else -1
             , nTupleOut.tagomega_os_cb if tagOSBounds else -1., nTupleOut.tagomega_ss_nn if tagSSBounds else -1. )

from ROOT import TObject
print 'writing n-tuple to file "%s"' % nTupleFilePathOut
nTupleFileOut.Write( nTupleFilePathOut, TObject.kOverwrite )
nTupleFileOut.Close()
