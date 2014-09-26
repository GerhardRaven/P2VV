runPeriod = 2012
cut = 'hlt2_biased==1'
if runPeriod : cut += ( ' && ' if cut else '' ) + 'runPeriod==%d' % runPeriod
dataNames = [ 'JpsiKK_sigSWeight' ] #[ 'JpsiKK', 'JpsiKK_sigSWeight', 'JpsiKK_cbkgSWeight' ]
removeObs = [ 'wMC', 'mdau1', 'tagCatP2VV' ] #, 'polarity', 'hlt2_prescale', 'nPVCat', 'BpTCat' ]
dataFilePathIn  = 'P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root'
dataFilePathOut = 'P2VVDataSets2012Reco14_I2Mass_6KKMassBins_2TagCats_20140309.root'

import P2VV.RooFitWrappers
from ROOT import TObject, TFile, RooFit, RooDataSet, RooArgSet, RooCategory
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
    if runPeriod :
        newArgSet.remove( newArgSet.find('runPeriod') )
        rp = RooCategory( 'runPeriod', 'runPeriod' )
        rp.defineType( 'p%d' % runPeriod, runPeriod )
        newArgSet.add(rp)

    newData.append( RooDataSet( dataName, dataName, newArgSet, RooFit.Import(data), RooFit.Cut(cut) ) )
    newData[-1].Print()
    newDataFile.Add( newData[-1] )

print 'write dataset to file "%s"' % newDataFile.GetName()
newDataFile.Write( dataFilePathOut, TObject.kOverwrite )
