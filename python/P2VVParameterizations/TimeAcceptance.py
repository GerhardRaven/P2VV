from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin
from P2VVParameterizations.GeneralUtils import valid_combinations, exclusive_combinations

class TimeAcceptance ( _util_parse_mixin, _util_extConstraints_mixin ) :
    def __init__( self, **kwargs ) : 
        if 'Acceptance' in kwargs:
            self._acceptance = kwargs.pop( 'Acceptance' )
        else:
            raise KeyError('TimeAcceptance: please specify an acceptance')
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw(kwargs)

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def acceptance( self ) : return self._acceptance
    def __mul__(self,rhs) : return self.acceptance() * rhs


class LP2011_TimeAcceptance ( TimeAcceptance ) :
    def __init__(self, **kwargs ) :
        from RooFitWrappers import ConstVar, FormulaVar, BinnedPdf
        self._parseArg( 'time',      kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._a = ConstVar(Name = 'eff_a', Value =  1.45 )
        self._c = ConstVar(Name = 'eff_c', Value = -2.37 )
        self._eff = FormulaVar('eff_shape', "(@0 > 0.) ? (1.0 / (1.0 + (@1 * @0) ** (@2))) : 0.0", [self._time, self._a, self._c])

        from P2VVBinningBuilders import build1DVerticalBinning
        self._binning, self._eff_func = build1DVerticalBinning('time_binning', self._eff, self._time, 0.05, 1.)

        TimeAcceptance.__init__( self, Acceptance = BinnedPdf(Name = 'time_acceptance', Observable = self._time, Function = self._eff, Binning = self._binning))

class Moriond2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        from RooFitWrappers import HistFunc
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       Value = 0., MinMax = (0.2, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histogram = kwargs.pop('Histogram', 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
        acceptance_file = TFile.Open(input_file)
        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        self._hist = acceptance_file.Get(histogram)
        if not self._hist:
            raise ValueError, 'Cannot get acceptance historgram %s from file' % histogram
        TimeAcceptance.__init__( self, Acceptance = HistFunc(self._hist.GetName(), Histogram = self._hist, Observables = [self._time]))
        acceptance_file.Close()

class Paper2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        from RooFitWrappers import HistFunc
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       Value = 0., MinMax = (0.2, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histograms = kwargs.pop('Histograms')
        acceptance_file = TFile.Open(input_file)
        acceptance_name = kwargs.pop('Name', 'Paper2012_Acceptance')

        data = kwargs.pop('Data')
        cuts = ' || '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in histograms.iterkeys()])
        data = data.reduce(cuts)
        total = data.sumEntries()

        exclusive = exclusive_combinations(histograms.keys())
        rel_spec = {}
        ## We assume that the categories have been constructed to be mutually exclusive.
        for comb in exclusive:
            cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
            rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}

        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        self._hists = {}
        for cat, name in histograms.iteritems():
            h = acceptance_file.Get(name)
            if not h:
                raise ValueError, 'Cannot get acceptance historgram %s from file %s' % (histogram, input_file)
            self._hists[cat] = h
        from collections import defaultdict
        bin_spec = defaultdict(dict)
        from array import array
        for (cat, label), hist in self._hists.iteritems():
            xaxis = hist.GetXaxis()
            bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)))
            heights = [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)]
            d = dict(bins = bins, heights = heights)
            bin_spec[cat][label] = d
        ## FIXME: make sure all the bins are set constant if needed
        TimeAcceptance.__init__( self, Acceptance = dict(Bins = bin_spec, Relative = rel_spec, Observable = self._time, ConditionalCategories = True, Name = acceptance_name, FitAcceptance = False))
        acceptance_file.Close()
