cut = 'hlt2_biased==1'
dataName = 'DecayTree'
dataFilePathIn  = 'fitNTuple_peakBkg_2011_2012_Reco14_TOS_20140822.root'
dataFilePathOut = 'fitNTuple_peakBkg_2011_2012_Reco14_TOS_HLT2B_20140822.root'

from ROOT import TFile
dataFile = TFile.Open(dataFilePathIn)
data = dataFile.Get(dataName)
print 'read n-tuple with %d entries from "%s"' % ( data.GetEntries(), dataFile.GetName() )

from ROOT import TObject
newDataFile = TFile.Open( dataFilePathOut, 'RECREATE' )
if cut :
    newData = data.CopyTree(cut)
else :
    newData = data.CloneTree()
print 'write n-tuple with %d entries to "%s"' % ( newData.GetEntries(), newDataFile.GetName() )
newDataFile.Write( dataFilePathOut, TObject.kOverwrite )
