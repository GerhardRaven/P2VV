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
        from array import array
        from RooFitWrappers import BinnedPdf
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       Value = 0., MinMax = (0.2, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histogram = kwargs.pop('Histogram', 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
        binning_name = kwargs.pop('BinningName', 'efficiency_binning')
        name = kwargs.pop('Name', 'Moriond2012_Acceptance')

        acceptance_file = TFile.Open(input_file)
        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        self._hist = acceptance_file.Get(histogram)
        if not self._hist:
            raise ValueError, 'Cannot get acceptance historgram %s from file' % histogram
        xaxis = self._hist.GetXaxis()
        bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, self._hist.GetNbinsX() + 2)))
        heights = [self._hist.GetBinContent(i) for i in range(1, self._hist.GetNbinsX() + 1)]

        acceptance_file.Close()

        from RooFitWrappers import RealVar
        heights = [RealVar('%s_bin_%03d' % (name, i + 1), Observable = False, Value = v,
                           Constant = True) for i, v in enumerate(heights)]
        # Add a binning for this category and state
        from ROOT import RooBinning
        self._binning = RooBinning(len(bins) - 1, bins)
        self._binning.SetName(binning_name)
        self._time.setBinning(self._binning, binning_name)

        from RooFitWrappers import BinnedPdf
        self._shape = BinnedPdf(name, Observable = self._time,
                                Binning = binning_name, Coefficients = heights)
        TimeAcceptance.__init__(self, Acceptance = self._shape)

class Paper2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        from RooFitWrappers import HistFunc
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       Value = 0., MinMax = (0.2, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histograms = kwargs.pop('Histograms')
        acceptance_file = TFile.Open(input_file)
        fit = kwargs.pop('Fit')
        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file

        acceptance_name = kwargs.pop('Name', 'Paper2012_FitAcceptance')

        data = kwargs.pop('Data')

        from collections import defaultdict
        levels = defaultdict(list)
        combinations = [(state, label) for state, info in histograms.iteritems() for label in info.iterkeys()]
        for comb in combinations:
            level = comb[0].GetName()[ : 4]
            levels[level].append(comb)
        cuts = ' || '.join(['{0} == {0}::{1}'.format(s.GetName(), l) for s, l in combinations])
        data = data.reduce(cuts)
        total = data.sumEntries()

        valid = valid_combinations(levels.itervalues())
        rel_spec = {}
        for comb in valid:
            cuts = ' && '.join(['{0} == {0}::{1}'.format(state.GetName(), label) for state, label in comb])
            rel_spec[comb] = {'Value' : data.sumEntries(cuts) / total, "Constant" : True}

        from array import array
        bin_spec = defaultdict(dict)

        for cat, label_info in histograms.iteritems():
            for label, info in label_info.iteritems():
                histogram = info.pop('histogram', None)
                if histogram:
                    hist = acceptance_file.Get(histogram)
                    if not hist:
                        raise ValueError, 'Cannot get acceptance histrogram %s from file %s' % (histogram, input_file)
                    xaxis = hist.GetXaxis()
                    bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)))
                    heights = [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)]
                else:
                    bins = array('d', (b for b in info['bins']))
                    heights = array('d', (h for h in info['heights']))
                d = dict(bins = bins, heights = heights)
                if 'average' in info:
                    d['average'] = info['average']
                bin_spec[cat][label] = d
        acceptance_file.Close()
        
        ## FIXME: make sure all the bins are set constant if needed
        TimeAcceptance.__init__( self, Acceptance = dict(Bins = bin_spec, Relative = rel_spec,
                                                         Observable = self._time,
                                                         ConditionalCategories = True,
                                                         Name = acceptance_name,
                                                         FitAcceptance = fit,
                                                         UseSingleBinConstraint = False))
