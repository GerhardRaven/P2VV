filePathIn  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root'
filePathOut = 'P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats.root'
dataSetName = 'JpsiKK_sigSWeight'
weightName  = 'sWeights_ipatia'
runPeriods  = [ 2011, 2012 ]

from P2VV import RooFitDecorators
from ROOT import TFile
print 'reading data set "%s" from file "%s"' % ( dataSetName, filePathIn )
fileIn = TFile.Open(filePathIn)
dsIn = fileIn.Get(dataSetName)
obsSetIn = dsIn.get()
dsIn.Print()

runPeriodCat = obsSetIn.find('runPeriod')
for period in runPeriods :
    assert runPeriodCat.isValidIndex(period)
if len(runPeriods) < runPeriodCat.numTypes() :
    from ROOT import RooCategory
    runPeriodCatOut = RooCategory( runPeriodCat.GetName(), runPeriodCat.GetTitle() )
    for catType in runPeriodCat :
        if catType.getVal() in runPeriods : runPeriodCatOut.defineType( catType.GetName(), catType.getVal() )

    from ROOT import RooArgSet
    obsSetOut = RooArgSet(obsSetIn)
    obsSetOut.remove(runPeriodCat)
    obsSetOut.add(runPeriodCatOut)
else :
    obsSetOut = obsSetIn

cuts = '||'.join( 'runPeriod==%d' % period for period in runPeriods )

from ROOT import RooDataSet, RooFit
print 'creating new data set for run periods %s' % runPeriods
print 'applying cuts "%s"' % cuts
dsOut = RooDataSet( dsIn.GetName(), dsIn.GetTitle(), obsSetOut, Import = dsIn, Cut = cuts, WeightVar = ( weightName, True ) )
dsOut.Print()

print 'writing data set to file "%s"' % filePathOut
fileOut = TFile.Open( filePathOut, 'RECREATE' )
fileOut.Add(dsOut)
fileOut.Write()
fileOut.Close()
fileIn.Close()
