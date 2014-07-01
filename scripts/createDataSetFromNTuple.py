protoFilePathIn = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets2011Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
nTupleFilePathIn = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/DS2011/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_sel.root'
dataFilePathOut = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/DS2011/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_%03d.root'
runPeriods = [ 2011 ]
KKMassBins = [ 990., 1050. ]
removeVars = [ 'tagCatP2VV' ]
selection = 'hlt2_biased==1'
numDataSets = 40
weightName = ''

import P2VV.RooFitWrappers
from ROOT import TFile, RooArgSet
protoFile = TFile.Open(protoFilePathIn)
nTupleFileIn = TFile.Open(nTupleFilePathIn)
protoData = protoFile.Get('JpsiKK_sigSWeight')
nTupleIn = nTupleFileIn.Get('DecayTree')

from ROOT import RooRealVar, RooCategory
obsSet = RooArgSet( protoData.get() )
if runPeriods :
    rp = obsSet.find('runPeriod')
    if rp : obsSet.remove(rp)
    rp = RooCategory( 'runPeriod', 'runPeriod' )
    for per in runPeriods : rp.defineType( 'p%d' % per, per )
    obsSet.add(rp)

if KKMassBins :
    KKCat = obsSet.find('KKMassCat')
    if KKCat : obsSet.remove(KKCat)
    KKCat = RooCategory( 'KKMassCat', 'KKMassCat' )
    for ind in range( len(KKMassBins) - 1 ) : KKCat.defineType( 'bin%d' % ind, ind )
    obsSet.add(KKCat)

    from array import array
    KKBinsArray = array( 'd', KKMassBins )
    from ROOT import RooBinning
    KKBinning = RooBinning( len(KKMassBins) - 1, KKBinsArray, 'KKMassBinning' )
    obsSet.find('mdau2').setBinning( KKBinning, 'KKMassBinning' ) 

if removeVars :
    for var in removeVars : obsSet.remove( obsSet.find(var) )

if weightName :
    wVar = RooRealVar( weightName, weightName, 1. )
    obsSet.add(wVar)

import os, subprocess, tempfile
tmp = tempfile.NamedTemporaryFile()
tmpFileName = os.path.realpath(tmp.name)
tmp.close()

from ROOT import RooDataSet, TObject, gROOT
nEvTot = 0
nEvDS = nTupleIn.GetEntries() / numDataSets
for it in range(numDataSets) :
    startEv = it * nEvDS
    nEv = nEvDS if it < numDataSets - 1 else nTupleIn.GetEntries() - ( numDataSets - 1 ) * nEvDS
    tmpFile = TFile.Open( tmpFileName, 'RECREATE' )
    nTupleSplit = nTupleIn.CopyTree( selection, '', nEv, startEv )
    if weightName :
        dataOut = RooDataSet( protoData.GetName(), protoData.GetTitle(), obsSet, Import = nTupleSplit, WeightVar = ( weightName, True ) )
    else :
        dataOut = RooDataSet( protoData.GetName(), protoData.GetTitle(), obsSet, Import = nTupleSplit )
    tmpFile.Close()

    print 'produced dataset:'
    dataOut.Print()
    nEvTot += dataOut.numEntries()
    dataFileOut = TFile.Open( dataFilePathOut % it, 'RECREATE' )
    dataFileOut.Append(dataOut)
    dataFileOut.Write( dataFilePathOut % it, TObject.kOverwrite )
    dataFileOut.Close()
os.remove(tmpFileName)
print 'total number of events: %d' % nEvTot
