from P2VV.RooFitWrappers import *
from ROOT import TFile

obj = RooObject( workspace = 'w')
w = obj.ws()

f = TFile.Open("/stuff/PhD/p2vv/data/MC11a.root")
sdata = f.Get("data")
w.put(sdata)

t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = t_minmax)
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550),
             Ranges =  { 'leftsideband'  : ( None, 5330 )
                         , 'signal'        : ( 5330, 5410 )
                         , 'rightsideband' : ( 5410, None ) 
                         } )
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
t_diff = RealVar('time_diff', Title = 'time diff', Unit = 'ps', Observable = True, MinMax = (-10, 10))
observables = [t, m, mpsi, st, t_diff]


from ROOT import RooCBShape as CrystalBall

mean_one = RealVar('mean_one', Value = 0, MinMax = (-1, 1))
sigma_one = RealVar('sigma_one', Value = 0.03, MinMax = (0.1, 1))
gauss_one = Pdf(Name = 'gauss_one', Type = Gaussian, Parameters = [t_diff, mean_one, sigma_one])


mean_two = RealVar('mean_two', Value = 0, MinMax = (-1, 1))
sigma_two = RealVar('sigma_two', Value = 0.8, MinMax = (0.001, 2))
alpha_two = RealVar('alpha_two', Value = 2, MinMax = (0.1, 50))
n_two = RealVar('n_two', Value = 2, MinMax = (0.1, 50), Constant = True)
cb_two = Pdf(Name = 'cb_two', Type = CrystalBall, Parameters = [t_diff, mean, sigma_two,
                                                                alpha_two, n_two])
c_exp = RealVar('c_exp', Value = -0.0001, MinMax = (-1, -1e-6))
exp = Pdf(Name = 'exp', Type = Exponential, Parameters = [t_diff, c_exp])

frac_cb = RealVar('frac_cb', Value = 0.2, MinMax = (1e-6, 0.9999))
gauss_cb = SumPdf(Name = 'gauss_cb', PDFs = [gauss_one, cb_two], Yields = {'cb_two' : frac_cb})

gauss_comp = Component('gauss_comp', [gauss_pdf], Yield = (5e5, 1e5, 1e6))
exp_comp = Component('exp_comp', [exp], Yield = (1e3, 1, 1e5))

diff_pdf = buildPdf(Name = 'diff_pdf', Components = [gauss_comp, exp_comp], Observables = [t_diff])
diff_pdf.Print('t')

diff_result = diff_pdf.fitTo(sdata, **fitOpts)

frame = t_diff.Frame(Range = (-0.5, 0.5))
sdata.plotOn(frame)
diff_pdf.plotOn(frame)

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 500, 500)
frame.Draw()
