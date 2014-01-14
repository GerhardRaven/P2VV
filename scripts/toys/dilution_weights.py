from P2VV.RooFitWrappers import *
obj = RooObject( workspace = 'w')

mpsi = RealVar('mpsi', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3025, 3165))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1, 1))

from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
dln = DoubleLogNormal(st, frac_ln2 = dict(Value = 0.312508), k1 = dict(Value = 0.801757),
                      k2 = dict(Value = 1.37584), median = dict(Value = 0.0309409))
ln = dln.pdf()

# Resolution models
from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as Multi_TimeResolution
tres_args = dict(time = t, sigmat = st, PerEventError = True, Cache = True)
tres_1 = Multi_TimeResolution(Name = 'tres', ParNamePrefix = 'one',
                              ScaleFactors = [(2, 2.00), (1, 1.174)],
                              Fractions = [(2, 0.143)], **tres_args)
tres_2 = TimeResolution(Name = 'tres', ParNamePrefix = 'two',
                        timeResSigmaSF = dict(Value = 4.), **tres_args)

# Gaussians for Time
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
g1 = Prompt_Peak(t, tres_1.model(), Name = 'g1')

g2 = Prompt_Peak(t, tres_2.model(), Name = 'g2')

# Mass shapes
from P2VV.Parameterizations.MassPDFs import Background_PsiMass as PsiBkgPdf
bkg_m = PsiBkgPdf(mpsi, Name = 'bkg_mpsi')

from P2VV.Parameterizations.MassPDFs import DoubleCB_Psi_Mass as PsiMassPdf
sig_m = PsiMassPdf(mpsi, Name = 'psi_m', mpsi_alpha_1 = dict(Value = 2, Constant = True),
                   mpsi_sigma_sf = dict(Value = 2.3), mpsi_frac = dict(Value = 0.8))

one = Component('one', [g1.pdf(), sig_m.pdf(), ln], Yield = (25000, 100, 100000))
two = Component('two', [g2.pdf(), bkg_m.pdf(), ln], Yield = (25000, 100, 100000))

pdf = buildPdf(Name = 'pdf', Components = [one, two], Observables = [mpsi, t, st])
mass_pdf = buildPdf(Name = 'mass_pdf', Components = [one, two], Observables = [mpsi])

from P2VV.Utilities.Resolution import SplitSigmat
split = SplitSigmat('', st)
sigmat_cat = split.split_cats()[0]

from P2VV import Dilution

dilutions = []
da = RealVar('da', Observable = True, MinMax = (0.01, 1.1))
dft = RealVar('dft', Observable = True, MinMax = (0.01, 1.1))
from ROOT import RooDataSet
from ROOT import RooArgSet
test_data = RooDataSet("test_data", "test_data", RooArgSet([da, dft]))

from multiprocessing import Process
from multiprocessing import Queue

class Calculator(Process):
    def __init__(self, mass_pdf, pdf, sigmat_cat, t, st, n, n_event = 100000):
        Process.__init__(self)
        self.__pdf = pdf
        self.__mass_pdf = mass_pdf
        self.__st_cat = sigmat_cat
        self.__t = t
        self.__st = st
        self.__n = n
        from ROOT import RooFit
        self.__spec = pdf.prepareMultiGen(RooArgSet(mpsi, t, st), RooFit.NumEvents(n_event))
        self.__queue = Queue()

    def run(self):
        fitOpts = dict(NumCPU = 1, Timer = 1, Save = True, Minimizer = 'Minuit2',
                       Optimize = 2, Offset = True)
        i = 0
        while i < self.__n:
            data = self.__pdf.generate(self.__spec)
            for i in range(3):
                mass_result = self.__mass_pdf.fitTo(data, **fitOpts)
                if mass_result.status() == 0:
                    i += 1
                    break
            else:
                continue
            from P2VV.Utilities.SWeights import SData
            sData = SData(Pdf = self.__mass_pdf, Data = data, Name = 'MassSPlot')
            sdata = sData.data('one')
            st_cat = sdata.addColumn(self.__st_cat._target_())
            d_ft = Dilution.dilution_bins(sdata, self.__t, self.__st, st_cat, t_range = 2)
            d_a = Dilution.signal_dilution_dg(sdata, self.__st, 1.2, 0.2, 2)
            self.__queue.put((d_a, d_ft))
        self.__queue.put('done')

    def queue(self):
        return self.__queue

n_p = 1
calculators = []
n_toys = 10
n_t = n_p * [n_toys / n_p]
for i in range(n_toys % n_p):
    n_t[i] += 1

for n in n_t:
    c = Calculator(mass_pdf, pdf, sigmat_cat, t, st, n)
    calculators.append(c)
    c.start()

args = RooArgSet(da, dft)
while len(calculators):
    for i, calculator in enumerate(calculators):
        msg = calculator.queue().get()
        if msg == 'done':
            calculator.join()
            calculators.pop(i)
            continue
        else:
            d_a, d_ft = msg
        da.setVal(d_a[0])
        dft.setVal(d_ft[0])
        dft.setError(d_ft[1])
        test_data.add(args)
    if test_data.numEntries() % 100 == 0:
        print 'completed ', test_data.numEntries()
        
diff = FormulaVar(Name = 'diff', Formula = '@0 - @1', Arguments = (dft, da), data = test_data)
from ROOT import TFile
f = TFile("dilution.root", "recreate")
f.WriteTObject(test_data, "data")
f.Close()
