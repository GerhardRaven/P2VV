from P2VV.RooFitWrappers import *

obj = RooObject( workspace = 'w')
w = obj.ws()

mass = RealVar( 'B_s0_LOKI_MASS_JpsiConstr',  Title = 'm(J/#psi K^{+}K^{-})', Unit = 'MeV/c^{2}', Observable = True
                 , Value = 5368., MinMax = ( 5200., 5550. ),
                 Ranges = dict(  LeftSideBand  = ( 5200., 5320. )
                                 , Signal        = ( 5320., 5420. )
                                 , RightSideBand = ( 5420., 5550. )
                                 , PeakBkg       = ( 5390., 5440. )
                                 ))
st = RealVar('sigmat',Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12))
t  = RealVar('time', Title = 'decay time', Unit='ps', Observable = True, MinMax = (0.3, 14))

excl_biased = Category('hlt1_excl_biased_dec', Observable = True,
                       States = { 'excl_biased' : 1, 'unbiased' : 0 } )
data_type = "2011_signal"

if data_type == "2011_prescaled":
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhiPrescaled_2011_Reco14_ntupleAB_20140430.root'
    cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
elif data_type == "2011_signal":
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2011/Bs2JpsiPhi_2011_Reco14_ntupleAB_20130906.root'
    cut = '%s==1 && (%s==1 || %s==1) && %s==1 && ' % ('sel', 'hlt1_biased', 'hlt1_unbiased_dec', 'hlt2_biased')
elif data_type == "2012_prescaled":
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2012/Bs2JpsiPhiPrescaled_2012_ntupleAB_20130419.root'
    cut = 'sel == 1 && triggerDecisionUnbiasedPrescaled == 1 && '
elif data_type == "2012_signal":
    input_file = '/glusterfs/bfys/users/raaij/NTuples/2012/Bs2JpsiPhi_2012_ntupleAB_20130604.root'
    cut = '%s==1 && (%s==1 || %s==1) && %s==1 && ' % ('sel', 'hlt1_biased', 'hlt1_unbiased_dec', 'hlt2_biased')
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
if data_type.endswith("signal"):
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
from P2VV.Utilities.DataHandling import correctWeights
sdata = SData( Pdf = sweight_pdf, Data = data, Name = 'mass_sdata')

sig_sdata = correctWeights(sdata.data('signal'), splitCatNames = [excl_biased.GetName()],
                           ImportIntoWS = False)
bkg_sdata = correctWeights(sdata.data('background'), splitCatNames = [excl_biased.GetName()],
                           ImportIntoWS = False)

from array import array
weights = array('d')

from ROOT import TFile
f = TFile.Open(input_file, "read")
tree = f.Get("DecayTree")

from ROOT import gDirectory
tree.Draw(">>elist", cut, "entrylist")

output = "temp.root"
output_file = TFile(output, "recreate")
copy_tree = tree.CopyTree(cut)
copy_tree.Write()

j = 0
for i in range(copy_tree.GetEntries()):
    r = sig_sdata.get(j)
    weights.append(sig_sdata.weight())

from ROOT import addSWeightToTree
addSWeightToTree(weights, len(weights), copy_tree, 'sweight')
copy_tree.Write()
output_file.Close()

