from ROOT import TFile
f = TFile("/stuff/PhD/p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats.root")
data = f.Get('JpsiKK_sigSWeight')
obs = data.get()
time = obs.find('time')

from P2VV.Utilities.General import make_binning
binning = make_binning(data, time, 15)

from ROOT import RooBinning
from array import array
rbins = RooBinning(len(binning) - 1, array('d', binning))

from P2VV.Parameterizations import TimeAcceptance
bc = TimeAcceptance.BinCounter(data, time, rbins)

bc.add_bins('hlt1', {'runPeriod' : 'p2011', 'hlt1_excl_biased_dec' : 'exclB'})
bc.add_bins('hlt1', {'runPeriod' : 'p2011', 'hlt1_excl_biased_dec' : 'notExclB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2011', 'hlt2_biased' : 'B', 'hlt2_unbiased' : 'UB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2011', 'hlt2_biased' : 'B', 'hlt2_unbiased' : 'notUB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2011', 'hlt2_biased' : 'notB', 'hlt2_unbiased' : 'UB'})

bc.add_bins('hlt1', {'runPeriod' : 'p2012', 'hlt1_excl_biased_dec' : 'exclB'})
bc.add_bins('hlt1', {'runPeriod' : 'p2012', 'hlt1_excl_biased_dec' : 'notExclB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2012', 'hlt2_biased' : 'B', 'hlt2_unbiased' : 'UB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2012', 'hlt2_biased' : 'B', 'hlt2_unbiased' : 'notUB'})
bc.add_bins('hlt2', {'runPeriod' : 'p2012', 'hlt2_biased' : 'notB', 'hlt2_unbiased' : 'UB'})

bc.run()
