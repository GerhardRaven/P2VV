from ROOT import TFile
nTupleFile = TFile.Open('/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2011_2012_s20_dv33r6p1_20140213_tupleB_selTrig_TOS.root')
nTuple = nTupleFile.Get('DecayTree')

runs = [ ]
runDict = { }
for ev in nTuple :
    num = ev.runNumber
    runs.append(num)
    if not num in runDict :
        runDict[num] = dict( nCand = 0, nHLT1UB = 0, nHLT1ExclB = 0, nHLT2UB = 0, nHLT2B = 0, nHLT2UBB = 0, magPol = ev.polarity )

    runDict[num]['nCand'] += 1
    if ev.hlt1_unbiased == 1 :
        runDict[num]['nHLT1UB'] += 1
    elif ev.hlt1_excl_biased == 1 :
        runDict[num]['nHLT1ExclB'] += 1
    if ev.hlt2_unbiased == 1 :
        runDict[num]['nHLT2UB'] += 1
        if ev.hlt2_biased == 1 :
            runDict[num]['nHLT2B'] += 1
            runDict[num]['nHLT2UBB'] += 1
    elif ev.hlt2_biased == 1 :
        runDict[num]['nHLT2B'] += 1
    assert ev.polarity == runDict[num]['magPol']
runs.sort()
print '%d runs read from file %s' % ( len(runDict), nTupleFile.GetName() )
print 'found run numbers in the range [%d, %d]' % ( runs[0], runs[-1] )
