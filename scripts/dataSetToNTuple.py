dsFilePath = 'P2VVDataSets20112012Reco14_noMC_I2Mass_6KKMassBins_2TagCats_HLT2B.root'
ntupleFilePath = 'fitNTuple_2011_2012_Reco14_20140128.root'

from ROOT import TFile
dsFile = TFile.Open(dsFilePath)
dataSet = dsFile.Get('JpsiKK_sigSWeight')
dsFile.Close()

import P2VV.RooFitWrappers
from ROOT import TObject
ntupleFile = TFile.Open( ntupleFilePath, 'RECREATE' )
ntuple = dataSet.buildTree( Name = 'DecayTree', Title = 'DecayTree', WeightName = 'sWeights_ipatia', RooFitFormat = False )
ntupleFile.Write( ntupleFilePath, TObject.kOverwrite )
ntupleFile.Close()
