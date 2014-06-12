#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog data_type')
(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)

from P2VV.Utilities.Resolution import input_data, prefix
if args[0] not in input_data.keys():
    print parser.print_usage()
    print "Possible samples are: %s" % ' '.join(input_data.keys())
    sys.exit(-2)

from P2VV.RooFitWrappers import *
from P2VV.Load import P2VVLibrary
from P2VV.Load import LHCbStyle

from ROOT import RooMsgService

obj = RooObject( workspace = 'w')
w = obj.ws()

from math import pi
t  = RealVar('time', Title = 't', Unit='ps', Observable = True, MinMax=(-5, 14))
m  = RealVar('mass', Title = 'B mass', Unit = 'MeV', Observable = True, MinMax = (5200, 5550))
mpsi = RealVar('mdau1', Title = 'J/psi mass', Unit = 'MeV', Observable = True, MinMax = (3030, 3150))
mphi = RealVar('mdau2', Title = 'phi mass', Unit = 'MeV', Observable = True, MinMax = (990, 1050))
st = RealVar('sigmat',Title = '#sigma_{t}', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
nPV = RealVar('nPV', Title = 'N_{PV}', Observable = True, MinMax = (0, 10))

from math import pi
cpsi = RealVar('helcosthetaK', Title = 'cpsi', Observable = True, MinMax = (-1, 1))
ctheta = RealVar('helcosthetaL', Title = 'ctheta', Observable = True, MinMax = (-1, 1))
phi = RealVar('helphi', Title = 'helphi', Observable = True, MinMax = (-pi, pi))
angles = [cpsi, ctheta, phi]

# add 20 bins for caching the normalization integral
for i in [ st ] : i.setBins( 20 , 'cache' )

# Categories
excl_biased = Category('triggerDecisionBiasedExcl', States = {'ExclBiased' : 1, 'Unbiased' : 0})
hlt1_biased = Category('hlt1_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt1_unbiased = Category('hlt1_unbiased_dec', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
hlt2_biased = Category('hlt2_biased', States = {'biased' : 1, 'not_biased' : 0}, Observable = True)
hlt2_unbiased = Category('hlt2_unbiased', States = {'unbiased' : 1, 'not_unbiased' : 0}, Observable = True)
selected = Category('sel', States = {'Selected' : 1, 'NotSelected' : 0})

observables = [t, m, mpsi, mphi, st, excl_biased, selected, nPV,
               hlt1_biased, hlt1_unbiased, hlt2_biased, hlt2_unbiased] + angles
               
project_vars = [st, excl_biased]

# now build the actual signal PDF...
from ROOT import RooGaussian as Gaussian
from ROOT import RooExponential as Exponential
from ROOT import RooDecay as Decay

# B time PDF
from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
tres_args = dict(time = t, sigmat = st, Cache = True,
                 Parameterise = 'RMS',
                 TimeResSFParam = 'linear', SplitFracs = False,
                 ScaleFactors = [(2, 2.1), (1, 1.26)],
                 Fractions = [(2, 0.2)], SplitMean = False)
tres = TimeResolution(Name = 'tres', **tres_args)

from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time as TimePdf
sig_t = TimePdf(t, tres.model())

# B mass pdf
from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as MassPdf
sig_m = MassPdf(m, Name = 'sig_m', AvSigParameterisation = True)
sig_m = sig_m.pdf()

# background
psi_c = RealVar( 'psi_c',  Unit = '1/MeV', Value = -0.0004, MinMax = (-0.1, -0.0000001))
bkg_m = Pdf(Name = 'bkg_m',  Type = Exponential, Parameters = [m, psi_c])

background = Component('background', (bkg_m,), Yield= (200000,500,5000000) )
signal = Component('signal', (sig_m, sig_t.pdf()), Yield= (200000,500,500000) )

from P2VV.Utilities.DataHandling import readData
tree_name = 'DecayTree'

## cut = 'runNumber > 0 && sel == 1 && (hlt1_biased == 1 || hlt1_unbiased_dec == 1) && hlt2_biased == 1 && '
cut = 'runNumber > 0 && sel == 1 && hlt1_unbiased_dec == 1 && hlt2_unbiased == 1 && '
cut += ' && '.join(['%s < 4' % e for e in ['muplus_track_chi2ndof', 'muminus_track_chi2ndof', 'Kplus_track_chi2ndof', 'Kminus_track_chi2ndof']])

from ROOT import RooAbsData
storage = RooAbsData.getDefaultStorageType()
RooAbsData.setDefaultStorageType(RooAbsData.Tree)

data = readData(input_data[args[0]]['data'], tree_name, ntupleCuts = cut,
                NTuple = True, observables = observables)

from ROOT import TFile
tmp_file = TFile.Open(os.path.join(prefix, 'p2vv/data/tmp.root'), 'recreate')
data.store().tree().SetDirectory(tmp_file)
data_cut = data.reduce("time > 0.3")
data_cut.store().tree().SetDirectory(tmp_file)

## Build PDF
mass_pdf = buildPdf(Components = (background, signal), Observables = (m, ), Name = 'mass_pdf')
mass_pdf.Print("t")

## Fit options
fitOpts = dict(NumCPU = 4, Timer = 1, Save = True,
               Verbose = False, Optimize = 1, Minimizer = 'Minuit2')

# make sweighted dataset. TODO: use mumu mass as well...
from P2VV.Utilities.SWeights import SData

result_cut = mass_pdf.fitTo(data_cut, **fitOpts)

__canvases = []

# Plot mass pdf
def plot_mass(data, prefix = ''):
    nBins = 100
    from ROOT import kDashed, kRed, kGreen, kBlue, kBlack
    from ROOT import TCanvas
    mass_canvas = TCanvas('%smass_canvas' % prefix, '%smass_canvas' % prefix, 600, 530)
    __canvases.append(mass_canvas)
    from P2VV.Utilities.Plotting import plot
    frames = plot(mass_canvas, m, pdf = mass_pdf, data = data,
                  dataOpts = dict(MarkerSize = 0.8, MarkerColor = kBlack, Binning = nBins),
                  pdfOpts  = dict(LineWidth = 2),
                  frameOpts = dict(Title = args[0]),
                  yTitleOffset = 1 / 0.7,
                  plotResidHist = True,
                  components = { 'bkg_*'     : dict( LineColor = kRed,   LineStyle = kDashed )
                                 , 'psi_*'  : dict( LineColor = kGreen, LineStyle = kDashed )
                                 , 'sig_*'     : dict( LineColor = kBlue,  LineStyle = kDashed )})
    
    frames[1].GetXaxis().SetTitle('M_{J/#psi K^{+}K^{-}} [MeV/#font[12]{c}^2]')
    frames[0].GetYaxis().SetTitle('Candidates / (%3.2f MeV/#font[12]{c}^2)' % ((t.getMax() - t.getMin()) / float(nBins)))

plot_mass(data_cut, 'cut_')

pdf_pars = mass_pdf.getParameters(data.get())
for p in sig_m.getParameters(data.get()):
    p = pdf_pars.find(p.GetName())
    p.setConstant(True)

result = mass_pdf.fitTo(data, **fitOpts)
plot_mass(data)

# Calculate sWeights
from P2VV.Utilities.SWeights import SData
sData = SData(Pdf = mass_pdf, Data = data, Name = 'MassSPlot')
sig_sdata = sData.data('signal')
sig_sdata.store().tree().SetDirectory(tmp_file)

bkg_sdata = sData.data('background')
bkg_sdata.store().tree().SetDirectory(tmp_file)

RooAbsData.setDefaultStorageType(storage)

# Load results from the resolution fit
from P2VV.CacheUtils import CacheFiles
directory = '1bin_19000.00fs_simple/7051155956274815792'

cache_files = CacheFiles(*input_data[args[0]]['cache'].rsplit('/', 1))
cache_dir, cache_file = cache_files.getFromCache(directory)
if not cache_dir:
    print 'Cannot find cached directory %s in %s' (cache_file.GetName(), directory)
    sys.exit(-2)

rd = cache_dir.Get('results')
time_result = rd.Get('time_result_double_RMS_Mixing_linear')
assert(time_result)

# PV bins
from array import array
PV_bounds = array('d', [-0.5 + i for i in range(12)])

from ROOT import TPaveText
year_label = TPaveText(0.71, 0.72, 0.89, 0.85, "NDC")
year_label.SetFillColor(0)
year_label.AddText(args[0].split('_')[0][-4:])
year_label.SetBorderSize(0)

from ROOT import TCanvas
canvas = TCanvas('time_canvas', 'time_canvas', 600, 400)
frame = t.frame()
sig_sdata.plotOn(frame)
frame.Draw()
frame.GetYaxis().SetRangeUser(1, 1000)
frame.GetYaxis().SetTitle('Candidates / (%3.2f ps)' % t.getBinning().averageBinWidth())
frame.GetXaxis().SetTitle('t [ps]')
canvas.SetLogy()
year_label.Draw()
from P2VV.Utilities.Resolution import plot_dir
canvas.Print(os.path.join(plot_dir, 'prescaled_sig_time_%s.pdf' % args[0]), EmbedFonts = True)

assert(False)

# Build WPV shape
from P2VV.Parameterizations import WrongPV
reweigh_data = dict(B = sig_sdata, bkg = bkg_sdata)
masses = {'B' : m}
weights = 'B'
wpv_builder = WrongPV.ShapeBuilder(t, masses, UseKeysPdf = True, Weights = weights, Draw = True,
                                   InputFile = input_data[args[0]]['wpv'],
                                   Workspace = input_data[args[0]]['workspace'],
                                   Reweigh = dict(Data = reweigh_data, DataVar = nPV,
                                                  Binning = PV_bounds),
                                   MassResult = result_cut)
wpv_signal = wpv_builder.shape('B')
sig_wpv = Component('sig_wpv', (wpv_signal, m), Yield = (888, 1, 300000))

# Build time PDF
time_pdf = buildPdf(Components = [signal, sig_wpv], Observables = [t], Name = 'time_pdf')


pdf_vars = time_pdf.getVariables()
for p in time_result.floatParsFinal():
    pdf_par = pdf_vars.find(p.GetName())
    if pdf_par and (pdf_par.GetName() not in ['N_sig_wpv']):
        pdf_par.setVal(p.getVal())
        pdf_par.setError(p.getError())
        pdf_par.setConstant(True)
        print 'Set value, error and constant for %s' % pdf_par.GetName()

time_pdf.Print('t')

time_result = time_pdf.fitTo(sig_sdata, SumW2Error = False, **fitOpts)
