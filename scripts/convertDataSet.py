dataNames = [ 'JpsiKK', 'JpsiKK_sigSWeight', 'JpsiKK_cbkgSWeight' ]
cut = 'hlt2_biased==1'
removeObs = [ 'wMC', 'mdau1', 'tagCatP2VV' ] #, 'polarity', 'hlt2_prescale', 'nPVCat', 'BpTCat' ]
dataFilePathIn  = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_4TagCats.root'
dataFilePathOut = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_4TagCats_HLT2B.root'

import P2VV.RooFitWrappers
from ROOT import TObject, TFile, RooFit, RooDataSet, RooArgSet
dataFile = TFile.Open(dataFilePathIn)
newDataFile = TFile.Open( dataFilePathOut, 'RECREATE' )
newData = [ ]
print 'read datasets from file "%s"' % dataFile.GetName()
for dataName in dataNames :
    print 'reading dataset "%s"' % dataName
    data = dataFile.Get(dataName)
    data.Print()

    newArgSet = RooArgSet( data.get() )
    for name in removeObs : newArgSet.remove( newArgSet.find(name) )
    newData.append( RooDataSet( dataName, dataName, newArgSet, RooFit.Import(data), RooFit.Cut(cut) ) )
    newData[-1].Print()
    newDataFile.Add( newData[-1] )

print 'write dataset to file "%s"' % newDataFile.GetName()
newDataFile.Write( dataFilePathOut, TObject.kOverwrite )
