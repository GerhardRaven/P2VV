from P2VV.RooFitWrappers import *
obj = RooObject( workspace = 'w')

st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.01, 0.07))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1, 1))

from P2VV.Parameterizations.SigmatPDFs import DoubleLogNormal
dln = DoubleLogNormal(st, frac_ln2 = dict(Value = 0.312508), k1 = dict(Value = 0.801757),
                      k2 = dict(Value = 1.37584), median = dict(Value = 0.0309409))
ln = dln.pdf()

# Resolution models
from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
tres_args = dict(time = t, sigmat = st, PerEventError = True, Cache = True)
tres_1 = TimeResolution(Name = 'tres', ParNamePrefix = 'one',
                        timeResSigmaSF = dict(Value = 1.2), **tres_args)
tres_2 = TimeResolution(Name = 'tres', ParNamePrefix = 'two',
                        timeResSigmaSF = dict(Value = 2.), **tres_args)

# Gaussians
from P2VV.Parameterizations.TimePDFs import Prompt_Peak
g1 = Prompt_Peak(t, tres_1.model(), Name = 'g1')
g2 = Prompt_Peak(t, tres_2.model(), Name = 'g2')

one = Component('one', [g1.pdf(), ln])
frac2 = RealVar('frac2', Value = 0.2, MinMax = (0.01, 0.99))
two = Component('two', [g2.pdf(), ln], Yield = (frac2,))

pdf_one = buildPdf(Name = 'pdf_one', Components = [one], Observables = [t, st])
pdf = buildPdf(Name = 'pdf', Components = [one, two], Observables = [t, st])

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
    def __init__(self, pdf, sigmat_cat, t, st, n, n_event = 100000):
        Process.__init__(self)
        self.__pdf = pdf
        self.__st_cat = sigmat_cat
        self.__t = t
        self.__st = st
        self.__n = n
        from ROOT import RooFit
        self.__spec = pdf.prepareMultiGen(RooArgSet(t, st), RooFit.NumEvents(n_event))
        self.__queue = Queue()

    def run(self):
        for i in range(self.__n):
            data = pdf.generate(self.__spec)
            st_cat = data.addColumn(self.__st_cat._target_())
            d_ft = Dilution.dilution_bins(data, self.__t, self.__st, st_cat, t_range = 2)
            d_a = Dilution.signal_dilution_dg(data, self.__st, 1.2, 0.2, 2)
            self.__queue.put((d_a, d_ft))
        self.__queue.put('done')

    def queue(self):
        return self.__queue

n_p = 4
calculators = []
n_toys = 5000
n_t = n_p * [n_toys / n_p]
for i in range(n_toys % n_p):
    n_t[i] += 1

for n in n_t:
    c = Calculator(pdf, sigmat_cat, t, st, n)
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
