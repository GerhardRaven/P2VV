from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

mass = RealVar( 'B_s0_MM',  Title = 'm(J/#psi K^{+}K^{-})', Unit = 'MeV/c^{2}', Observable = True
                 , Value = 5368., MinMax = ( 5200., 5550. ),
                 Ranges = dict(  LeftSideBand  = ( 5200., 5320. )
                                 , Signal        = ( 5320., 5420. )
                                 , RightSideBand = ( 5420., 5550. )
                                 , PeakBkg       = ( 5390., 5440. )
                                 ))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0, 0.1))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1.5, 8))

excl_biased = Category('hlt1_excl_biased_dec', Observable = True,
                       States = { 'excl_biased' : 1, 'unbiased' : 0 } )
prescaled = False

if prescaled:
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhiPrescaled_ntupleAB_20130531.root'
    cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
else:
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhi_ntupleAB_20130531.root'
    cut = '%s==1 && (%s==1 || %s==1) && %s==1 && ' % ('sel', 'hlt1_biased', 'hlt1_unbiased',
                                                      'hlt2_biased')
cut += ' && '.join(['%s < 4' % e for e in ['muplus_TRACK_CHI2NDOF', 'muminus_TRACK_CHI2NDOF', 'Kplus_TRACK_CHI2NDOF', 'Kminus_TRACK_CHI2NDOF']])
cut += ' && sel_cleantail == 1'

for o in [mass, t, st]:
    cut += ' && {0} > {1} && {0} < {2}'.format(o.GetName(), o.getMin(), o.getMax())

from P2VV.GeneralUtils import readData
data = readData(filePath = input_file, dataSetName = 'DecayTree', NTuple = True,
                observables = [mass, t, st, excl_biased], Rename = 'JpsiphiData', ntupleCuts = cut)

sigFrac = 0.504
nEvents     = data.sumEntries()
nSignal     = nEvents * sigFrac
nBackground = nEvents * ( 1. - sigFrac )

from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass
sig_m = SignalBMass(Name = 'sig_m', mass = mass)

signal = Component( 'signal', [sig_m.pdf()], Yield = ( nSignal,     0., nEvents ) )

from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass
bkg_m = BackgroundBMass(Name = 'bkg_m', mass = mass)
background = Component( 'background', [bkg_m.pdf()], Yield = ( nBackground, 0., nEvents ) )

mass_pdf = buildPdf( [ signal, background ], Observables = [ mass ], Name = 'mass_pdf')

fitOpts = dict(NumCPU = 4, Optimize = 2, Save = True, Timer = True, Minimizer = 'Minuit2',
               Offset = True)

mass_result = mass_pdf.fitTo(data, **fitOpts)

# categories for splitting the PDF
if not prescaled:
    split_cats = [[excl_biased]]
    # get mass parameters that are split
    split_params = [[par for par in mass_pdf.Parameters() if par.getAttribute('Yield')]]

    # build simultaneous mass PDF
    from P2VV.RooFitWrappers import SimultaneousPdf
    sWeight_mass_pdf = SimultaneousPdf(mass_pdf.GetName() + '_simul'
                                       , MasterPdf       = mass_pdf
                                       , SplitCategories = split_cats
                                       , SplitParameters = split_params)

    # set yields for categories
    split_cat  = sWeight_mass_pdf.indexCat()
    split_vars = sWeight_mass_pdf.getVariables()
    from P2VV.GeneralUtils import getSplitPar
    from math import sqrt
    for state in split_cat:
        sigYield = getSplitPar( 'N_signal', state.GetName(), split_vars )
        bkgYield = getSplitPar( 'N_background', state.GetName(), split_vars )
        
        selStr = '{0} == {0}::{1}'.format(split_cat.GetName(), state.GetName())
        
        nEv    = data.sumEntries()
        nEvBin = data.sumEntries(selStr)
        
        sigYield.setVal( sigYield.getVal() * nEvBin / nEv )
        sigYield.setError( sqrt( sigYield.getVal() ) )
        sigYield.setMin(0.)
        sigYield.setMax(nEvBin)
        bkgYield.setVal( bkgYield.getVal() * nEvBin / nEv )
        bkgYield.setError( sqrt( bkgYield.getVal() ) )
        bkgYield.setMin(0.)
        bkgYield.setMax(nEvBin)
    
    # determine mass parameters in each sub-sample with a fit
    sim_mass_result = sWeight_mass_pdf.fitTo(data, **fitOpts)
    sweight_pdf = sWeight_mass_pdf
else:
    sweight_pdf = mass_pdf

from P2VV.GeneralUtils import SData
sdata = SData( Pdf = sweight_pdf, Data = data, Name = 'mass_sdata')
sig_sdata = sdata.data('signal')
bkg_sdata = sdata.data('background')

from ROOT import TFile
f = TFile.Open(input_file, "update")
tree = f.Get("DecayTree")

from ROOT import addVertexErrors
from ROOT import std
dss = std.list("RooDataSet*")()
dss.push_back(sig_sdata)
dss.push_back(bkg_sdata)
addVertexErrors(tree, dss, cut)

from P2VV.Load import LHCbStyle

from ROOT import TCanvas
canvases = [TCanvas(n, n, 600, 400) for n in ['sv_canvas', 'st_canvas', 'psi_canvas']]
canvases[0].cd()
sv_err = sig_sdata.get().find('sv_err')
frame = sv_err.frame()
frame.GetXaxis().SetTitle('#sigma_{t,SV} [ps]')
frame.GetYaxis().SetTitle('Candidates / fs')
sig_sdata.plotOn(frame)
frame.Draw()

canvases[1].cd()
frame = st.frame()
frame.GetXaxis().SetTitle('#sigma_{t} [ps]')
frame.GetYaxis().SetTitle('Candidates / fs')
sig_sdata.plotOn(frame)
frame.Draw()

canvases[2].cd()
psi_err = sig_sdata.get().find('jpsi_vx_err')
frame = psi_err.frame()
frame.GetXaxis().SetTitle('#sigma_{t,J/#psi} [ps]')
frame.GetYaxis().SetTitle('Candidates / fs')
sig_sdata.plotOn(frame)
frame.Draw()
