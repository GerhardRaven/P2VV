from P2VV.Load import LHCbStyle
from P2VV.RooFitWrappers import *
from ROOT import TFile

f = TFile("dilution.root")
data = f.Get("data")
diff = data.get().find("diff")

obj = RooObject( workspace = 'w')
obj.ws().put(data)

mean = RealVar('mean', Value = -0.001, MinMax = (-1, 1))
sigma = RealVar('sigma', Value = 0.0011, MinMax = (0.0001, 0.1))
from ROOT import RooGaussian
gauss = Pdf(Name = 'gauss', Type = RooGaussian, Parameters = [diff, mean, sigma])

result = gauss.fitTo(data, Minimizer = 'Minuit2', Save = True, Optimize = 2)

frame = diff.frame(Range = (-0.016, -0.004))
data.plotOn(frame, Binning = (40, -0.016, -0.004))
gauss.plotOn(frame)
frame.GetXaxis().SetTitle('D_{FT} - D')
frame.Draw()
