brNames = [  ( 'runNumber', 'i' ), ( 'eventNumber', 'l' ), ( 'runPeriod', 'I' ), ( 'polarity', 'I' ), ( 'nPV', 'D' )
           , ( 'sWeights_ipatia', 'F' ), ( 'wMC', 'F' )
           , ( 'hlt1_unbiased', 'I' ), ( 'hlt1_excl_biased', 'I' ), ( 'hlt2_unbiased', 'I' ), ( 'hlt2_biased', 'I' )
           , ( 'triggerDecisionUnbiased', 'I' ), ( 'triggerDecisionBiasedExcl', 'I' )
           , ( 'mass', 'D' ), ( 'mdau1', 'D' ), ( 'mdau2', 'D' )
           , ( 'time', 'D' ), ( 'sigmat', 'D' )
           , ( 'helcosthetaK', 'D' ), ( 'helcosthetaL', 'D' ), ( 'helphi', 'D' )
           , ( 'tagdecision_os_cb', 'I' ), ( 'tagdecision_ss_nn', 'I' ), ( 'sel', 'I' ), ( 'sel_cleantail', 'I' )
           , ( 'tagomega_os_cb', 'D' ), ( 'tagomega_ss_nn', 'D' )
           , ( 'B_P', 'D' ), ( 'B_Pt', 'D' ), ( 'B_eta', 'D' ), ( 'B_phi', 'D' )
           , ( 'Kplus_PX', 'D' ),   ( 'Kplus_PY', 'D' ),   ( 'Kplus_PZ', 'D' ),   ( 'Kplus_LOKI_ETA', 'D' )
           , ( 'Kminus_PX', 'D' ),  ( 'Kminus_PY', 'D' ),  ( 'Kminus_PZ', 'D' ),  ( 'Kminus_LOKI_ETA', 'D' )
           , ( 'muplus_PX', 'D' ),  ( 'muplus_PY', 'D' ),  ( 'muplus_PZ', 'D' ),  ( 'muplus_LOKI_ETA', 'D' )
           , ( 'muminus_PX', 'D' ), ( 'muminus_PY', 'D' ), ( 'muminus_PZ', 'D' ), ( 'muminus_LOKI_ETA', 'D' )
           , ( 'muplus_track_chi2ndof', 'D' ), ( 'muminus_track_chi2ndof', 'D' )
           , ( 'Kplus_track_chi2ndof', 'D' ), ( 'Kminus_track_chi2ndof', 'D' )
           , ( 'B_s0_ENDVERTEX_X', 'D' ), ( 'B_s0_ENDVERTEX_Y', 'D' ), ( 'B_s0_ENDVERTEX_Z', 'D' )
          ]
nTupleFilePathIn = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged.root'
nTupleFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged_presc_20140822.root'

from ROOT import TFile, TTree
nTupleFileIn = TFile.Open(nTupleFilePathIn)
nTupleIn = nTupleFileIn.Get('DecayTree')
nTupleFileOut = TFile.Open( nTupleFilePathOut, 'RECREATE' )
nTupleOut = TTree( 'DecayTree', 'DecayTree' )

from array import array
brDict = { }
typeDict = dict( i = 'L', l = 'L', I = 'l', F = 'f', D = 'd' )
for br in brNames :
    brDict[ br[0] ] = array( typeDict[ br[1] ], [ 0 ] )
    nTupleIn.SetBranchAddress( br[0], brDict[ br[0] ] )
    nTupleOut.Branch( br[0], brDict[ br[0] ], '%s/%s' % ( br[0], br[1] ) )

print 'reading events from file %s: %d entries' % ( nTupleFileIn.GetName(), nTupleIn.GetEntries() )
preScale = 0.2
runNumbers = { 2011 : { 0 : 87219, 1 : 95929 }, 2012 : { 0 : 111183, 1 : 119956 } }
lumiFracs = { 2011 : 0.4 / ( 0.4 + 0.7 ), 2012 : 0.6 / ( 0.6 + 1.5 ) }
from ROOT import TRandom2
rand = TRandom2(12345)
for ev in nTupleIn :
    if brDict['wMC'][0] < 0. :
        assert brDict['hlt2_unbiased'][0] == 1
        if rand.Rndm() > lumiFracs[ brDict['runPeriod'][0] ] :
            brDict['runNumber'][0] = runNumbers[ brDict['runPeriod'][0] ][1]
            if brDict['hlt2_biased'][0] == 1 :
                if rand.Rndm() > preScale : brDict['hlt2_unbiased'][0] = 0
            else :
                brDict['wMC'][0] *= preScale
        else :
            brDict['runNumber'][0] = runNumbers[ brDict['runPeriod'][0] ][0]

    nTupleOut.Fill()

print 'writing events to file %s: %d entries' % ( nTupleFileOut.GetName(), nTupleOut.GetEntries() )
from ROOT import TObject
nTupleFileOut.Write( nTupleFilePathOut, TObject.kOverwrite )
nTupleFileOut.Close()
nTupleFileIn.Close()
