from RooFitWrappers import RealVar
from RooFitWrappers import RooObject
from ROOT import TFile
from ROOT import RooDataSet
from ROOT import RooWorkspace
from ROOT import RooArgSet

ws = RooObject()
ws.setWorkspace(RooWorkspace("workspace"))

m = RealVar('data_m', Observable = True, Unit = 'MeV/c^2', MinMax=(5000, 6000))
t = RealVar('data_t', Observable = True, MinMax = (-1, 14), Unit = 'ps', Value = 0)

args = RooArgSet(m._target_(), t._target_())
d = RooDataSet('data', 'data', args)
f = TFile.Open('data.root', 'recreate')
f.WriteTObject(d, 'data')
f.Close()
