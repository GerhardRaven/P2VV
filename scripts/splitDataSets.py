from ROOT import TFile, RooArgSet
dataFile = TFile.Open('P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_kinematics.root')
data = dataFile.Get('JpsiKK_sigSWeight')
dataFile.Close()
varSet = RooArgSet( data.get() )
for var in [ 'wMC', 'hlt2_prescale', 'polarity', 'tagCatP2VV', 'nPVCat', 'BpTCat', 'runPeriod' ] :
    varSet.remove( varSet.find(var) )

from ROOT import RooCategory
varSet11 = RooArgSet(varSet)
varSet12 = RooArgSet(varSet)
runPeriod11 = RooCategory( 'runPeriod', 'runPeriod' )
runPeriod12 = RooCategory( 'runPeriod', 'runPeriod' )
runPeriod11.defineType( 'p2011', 2011 )
runPeriod12.defineType( 'p2012', 2012 )
varSet11.add(runPeriod11)
varSet12.add(runPeriod12)
from ROOT import RooFit, RooDataSet
data11 = RooDataSet( data.GetName(), data.GetTitle(), varSet11, RooFit.Import(data), RooFit.Cut('hlt2_biased==1 && runPeriod==2011') )
data12 = RooDataSet( data.GetName(), data.GetTitle(), varSet12, RooFit.Import(data), RooFit.Cut('hlt2_biased==1 && runPeriod==2012') )

from ROOT import TObject
dataFile11 = TFile.Open( 'P2VVDataSets2011Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root', 'RECREATE' )
dataFile11.Add(data11)
dataFile11.Write( 'P2VVDataSets2011Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root', TObject.kOverwrite )
dataFile11.Close()

dataFile12 = TFile.Open( 'P2VVDataSets2012Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root', 'RECREATE' )
dataFile12.Add(data12)
dataFile12.Write( 'P2VVDataSets2012Reco14_I2Mass_6KKMassBins_2TagCats_kinematics_HLT2B.root', TObject.kOverwrite )
dataFile12.Close()
