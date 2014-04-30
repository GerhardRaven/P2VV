from P2VV.RooFitWrappers import *
from P2VV.Load import LHCbStyle
from ROOT import TFile

obj = RooObject( workspace = 'w')
w = obj.ws()

t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (0.3, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
st = RealVar('sigmat',Title = '#sigma_{t}', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))

# add 20 bins for caching the normalization integral
st.setBins(50, 'cache')

from math import pi
cpsi = RealVar('helcosthetaK', Title = 'cpsi', Observable = True, MinMax = (-1, 1))
ctheta = RealVar('helcosthetaL', Title = 'ctheta', Observable = True, MinMax = (-1, 1))
phi = RealVar('helphi', Title = 'helphi', Observable = True, MinMax = (-pi, pi))
angles = [cpsi, ctheta, phi]

f = TFile("/project/bfys/raaij/p2vv/data/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_HLT2B_20140309.root")
data = f.Get("JpsiKK_sigSWeight")
data = data.reduce('runPeriod == runPeriod::p2012')
from ROOT import TCanvas
from collections import defaultdict
__store = defaultdict(list)
for obs, n_bins, x_title, y_title in zip([t, cpsi, ctheta, phi],
                                [137, 100, 100, 63],
                                ['t [ps]', 'cos #psi', 'cos #theta', '#phi [rad]'],
                                ['Candidates / (0.1 ps)', 'Candidates / 0.2' , 'Candidates / 0.2', 'Candidates / (0.1 rad)']):
    cn = obs.GetName() + '_canvas'
    canvas = TCanvas(cn, cn, 600, 400)
    __store['canvases'].append(canvas)
    frame = obs.frame(Bins = n_bins)
    frame.GetXaxis().SetTitle(x_title)
    frame.GetYaxis().SetTitle(y_title)    
    __store['frames'].append(frame)
    data.plotOn(frame)
    frame.Draw()

psi_cache = TFile("/project/bfys/raaij/p2vv/data/Bs2JpsiPhi_2012_Prescaled_cache.root")
psi_data = psi_cache.Get("1bin_15500.00fs_simple/2027465761870101697/sdata/full_sig_sdata")

n_sig = data.sumEntries()
n_psi = psi_data.sumEntries()

canvas = TCanvas('st_canvas', 'st_canvas', 600, 400)
frame = st.frame(Range = (0, 0.07))
from ROOT import kGreen, kBlue
data.plotOn(frame, MarkerColor = kBlue, Rescale = 1 / n_sig)
psi_data.plotOn(frame, MarkerColor = kGreen, Rescale = 1 / n_psi)
frame.Draw()
frame.GetYaxis().SetRangeUser(0, 0.045)
frame.GetYaxis().SetTitle("Candidates / (0.7 fs)")
frame.GetXaxis().SetTitle("#sigma_{t} [ps]")
