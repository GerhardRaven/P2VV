dataName = 'JpsiKK_sigSWeight'
dataFilePathIn  = 'P2VVDataSets20112012Reco14_noMC_DGMass_6KKMassBins_2TagCats.root'
dataFilePathOut = 'P2VVDataSets20112012Reco14_noMC_DGMass_6KKMassBins_2TagCats_HLT2B.root'

from ROOT import TFile
dataFile = TFile.Open(dataFilePathIn)
data = dataFile.Get(dataName)

import P2VV.RooFitWrappers
from ROOT import TObject, RooFit, RooDataSet
newData = RooDataSet( dataName, dataName, data.get(), RooFit.Import(data), RooFit.Cut('hlt2_biased==1') )
newDataFile = TFile.Open( dataFilePathOut, 'RECREATE' )
newDataFile.Add(newData)
newDataFile.Write( dataFilePathOut, TObject.kOverwrite )
