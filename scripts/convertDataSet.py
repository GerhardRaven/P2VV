dataName = 'JpsiKK_sigSWeight'
cut = 'hlt2_biased==1'
removeObs = [ 'wMC', 'firstData', 'hlt2_prescale', 'polarity', 'mdau1', 'tagCatP2VV' ]
dataFilePathIn  = 'P2VVDataSets20112012Reco14_noMC_DGMass_6KKMassBins_2TagCats.root'
dataFilePathOut = 'P2VVDataSets20112012Reco14_noMC_DGMass_6KKMassBins_2TagCats_HLT2B.root'

from ROOT import TFile
dataFile = TFile.Open(dataFilePathIn)
data = dataFile.Get(dataName)
print 'read dataset from file "%s"' % dataFile.GetName()
data.Print()

import P2VV.RooFitWrappers
from ROOT import TObject, RooFit, RooDataSet, RooArgSet
newArgSet = RooArgSet( data.get() )
for name in removeObs : newArgSet.remove( newArgSet.find(name) )
newData = RooDataSet( dataName, dataName, newArgSet, RooFit.Import(data), RooFit.Cut(cut) )
newData.Print()
newDataFile = TFile.Open( dataFilePathOut, 'RECREATE' )
newDataFile.Add(newData)
print 'write dataset to file "%s"' % newDataFile.GetName()
newDataFile.Write( dataFilePathOut, TObject.kOverwrite )
