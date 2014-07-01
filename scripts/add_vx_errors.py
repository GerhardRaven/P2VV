#!/usr/bin/env python
import optparse
import sys
import os
from math import sqrt

parser = optparse.OptionParser(usage = 'usage: %prog <2011|2012>')
parser.add_option("--dls", dest = "dls", default = False,
                  action = 'store_true', help = 'Add DLS')
parser.add_option("-p", "--prescaled", dest = "prescaled", default = False,
                  action = 'store_true', help = 'Use prescaled sample')

(options, args) = parser.parse_args()

if len(args) != 1:
    print parser.print_usage()
    sys.exit(-2)
elif args[0] not in ['2011', '2012']:
    print parser.print_usage()
    print "Year must be 2011 or 2012"
    sys.exit(-2)
else:
    year = args[0]

from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

mass = RealVar('B_s0_LOKI_MASS_JpsiConstr', Title = 'm(J/#psi K^{+}K^{-})', Unit = 'MeV/c^{2}', Observable = True,
               Value = 5368., MinMax = (5200., 5550.))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0, 0.1))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (-1.5, 8))

excl_biased = Category('hlt1_excl_biased_dec', Observable = True,
                       States = { 'excl_biased' : 1, 'unbiased' : 0 } )
nPV = RealVar('nPV', Title = 'Number of PVs', Observable = True, MinMax = (0, 10))

if options.prescaled:
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhiPrescaled_ntupleAB_20130531.root'
    cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
else:
    if year == '2011':
        input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhi_2011_Reco14_ntupleAB_20130906.root'
    elif year == '2012':
        input_file = '/glusterfs/bfys/users/raaij/NTuples/2012/Bs2JpsiPhi_2012_ntupleAB_20130604.root'
    cut = '%s==1 && (%s==1 || %s==1) && %s==1 && ' % ('sel', 'hlt1_biased', 'hlt1_unbiased',
                                                      'hlt2_biased')
cut += ' && '.join(['%s < 4' % e for e in ['muplus_TRACK_CHI2NDOF', 'muminus_TRACK_CHI2NDOF', 'Kplus_TRACK_CHI2NDOF', 'Kminus_TRACK_CHI2NDOF']])
cut += ' && sel_cleantail == 1'

for o in [mass, t, st]:
    cut += ' && {0} > {1} && {0} < {2}'.format(o.GetName(), o.getMin(), o.getMax())

from P2VV.Utilities.DataHandling import readData
data = readData(filePath = input_file, dataSetName = 'DecayTree', NTuple = True,
                observables = [mass, t, st, excl_biased], Rename = 'JpsiphiData', ntupleCuts = cut)

sigFrac = 0.504
nEvents     = data.sumEntries()
nSignal     = nEvents * sigFrac
nBackground = nEvents * ( 1. - sigFrac )

## from P2VV.Parameterizations.MassPDFs import Ipatia2_Signal_Mass as SignalBMass
## sig_m = SignalBMass(Name = 'sig_m', mass = mass, m_sig_lambda = -2.5, m_sig_zeta = 0.01,
##                     m_sig_alpha_1 = 3.0, m_sig_alpha_2 = 2.5, m_sig_n_1 = 1.5, m_sig_n_2 = 2.0)
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
if not options.prescaled:
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
    from P2VV.Utilities.General import getSplitPar
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

from P2VV.Utilities.SWeights import SData
sdata = SData( Pdf = sweight_pdf, Data = data, Name = 'mass_sdata')
sig_sdata = sdata.data('signal')
bkg_sdata = sdata.data('background')

from ROOT import TFile
f = TFile.Open(input_file, "update")
tree = f.Get("DecayTree")

from ROOT import addVertexErrors, addJpsiDLS
from ROOT import std
dss = std.list("RooDataSet*")()
dss.push_back(sig_sdata)
dss.push_back(bkg_sdata)
if options.dls:
    addJpsiDLS(tree, dss, cut)
else:
    addVertexErrors(tree, dss, cut)

output_file = "vx_errors.root"
option = "update" if os.path.exists(output_file) else "new"
f = TFile(output_file, option)
if not f.Get(year):
    f.mkdir(year)

d = f.Get(year)
for n, ds in {'sig_sdata' : sig_sdata, 'bkg_sdata' : bkg_sdata}.iteritems():    
    d.WriteTObject(ds, n, 'overwrite')

f.Write()
f.Close()
