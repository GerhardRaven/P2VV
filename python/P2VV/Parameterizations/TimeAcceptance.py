from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin
from P2VV.Parameterizations.GeneralUtils import valid_combinations, exclusive_combinations
from P2VV.Parameterizations.TimeResolution import TimeResolution

## Since all time acceptances are implemented in the resolution model, we inherit from there
class TimeAcceptance ( TimeResolution ) :
    def __init__( self, **kwargs ) : 
        if 'Acceptance' in kwargs: self._acceptance = kwargs.pop( 'Acceptance' )
        else: raise KeyError('P2VV - ERROR: TimeAcceptance.__init__(): TimeAcceptance: please specify an acceptance')
        # Workaround for acceptance normalization in cFit...
        # print 'P2VV - WARNING: TimeAcceptance.__init__(): set "NOCacheAndTrack" for %s' % self._acceptance.GetName()

        from P2VV.RooFitWrappers import BinnedPdf, EffResModel
        if type(self._acceptance) == EffResModel and type(self._acceptance.efficiency()) == BinnedPdf:
            # note: constant optimization WILL evaluate RooBinnedPdf as a PDF,
            # and thus normalize it...
            self._acceptance.efficiency().setForceUnitIntegral(True)
            print 'P2VV - WARNING: TimeAcceptance.__init__(): switched setForceUnitIntegral to true for %s' % self._acceptance.GetName()
        print self.acceptance()
        
        TimeResolution.__init__(self, Model = self._acceptance,
                                Conditionals = self._acceptance.ConditionalObservables(),
                                Constraints = self._acceptance.ExternalConstraints(),
                                **kwargs)

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def acceptance( self ) : return self._acceptance

class Moriond2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        from array import array
        from P2VV.RooFitWrappers import BinnedPdf
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       MinMax = (0.3, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histogram = kwargs.pop('Histogram', 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
        binning_name = kwargs.pop('BinningName', 'efficiency_binning')
        name = kwargs.pop('Name', 'Moriond2012_Acceptance')
        model = kwargs.pop('ResolutionModel')

        acceptance_file = TFile.Open(input_file)
        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        self._hist = acceptance_file.Get(histogram)
        if not self._hist:
            raise ValueError, 'Cannot get acceptance histogram %s from file' % histogram
        xaxis = self._hist.GetXaxis()
        bins = array('d', (xaxis.GetBinLowEdge(i) for i in range(1, self._hist.GetNbinsX() + 2)))
        self._binHeights = [self._hist.GetBinContent(i) for i in range(1, self._hist.GetNbinsX() + 1)]

        acceptance_file.Close()
        print 'P2VV - INFO: Moriond2012_TimeAcceptance.__init__(): using time efficiency histogram "%s" from file "%s"'\
              % ( histogram, input_file )

        from P2VV.RooFitWrappers import RealVar
        self._binHeights = [RealVar('%s_bin_%03d' % (name, i + 1), Observable = False, Value = v,
                           Constant = True) for i, v in enumerate(self._binHeights)]
        # Add a binning for this category and state
        from ROOT import RooBinning
        self._binning = RooBinning(len(bins) - 1, bins)
        self._binning.SetName(binning_name)
        self._time.setBinning(self._binning, binning_name)

        from P2VV.RooFitWrappers import BinnedPdf
        self._shape = BinnedPdf(name + '_shape', Observable = self._time,
                                Binning = binning_name, Coefficients = self._binHeights)
        from P2VV.RooFitWrappers import EffResModel
        TimeAcceptance.__init__(self, Acceptance = EffResModel(Name = name, Efficiency = self._shape,
                                                               ResolutionModel = model['model'],
                                                               ConditionalObservables = model.conditionalObservables(),
                                                               ExternalConstraints = model.externalConstraints()))

class Paper2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       MinMax = (0.3, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histograms = kwargs.pop('Histograms')
        original = kwargs.pop('Original')
        acceptance_file = TFile.Open(input_file)
        fit = kwargs.pop('Fit')
        model = kwargs.pop('ResolutionModel')
        binHeightMinMax = kwargs.pop('BinHeightMinMax', None)

        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        print 'P2VV - INFO: Paper2012_TimeAcceptance.__init__(): using time efficiency histograms file "%s"' % input_file

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

        valid = valid_combinations(levels.values())
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
                    print 'P2VV - INFO: Paper2012_TimeAcceptance.__init__(): using time efficiency histogram "%s"' % histogram
                else:
                    bins = array('d', (b for b in info['bins']))
                    heights = array('d', (h for h in info['heights']))
                d = dict(bins = bins, heights = heights)
                if 'average' in info:
                    d['average'] = info['average']
                bin_spec[cat][label] = d
        acceptance_file.Close()

        from P2VV.RooFitWrappers import MultiHistEfficiencyModel
        ## FIXME: make sure all the bins are set constant if needed
        mhe = MultiHistEfficiencyModel(Bins = bin_spec, Relative = rel_spec, Observable = self._time,
                                       ConditionalCategories = True, Name = acceptance_name,
                                       FitAcceptance = fit, UseSingleBinConstraint = False,
                                       ResolutionModel = model['model'], Original = original,
                                       ConditionalObservables = model.conditionalObservables(),
                                       ExternalConstraints = model.externalConstraints(),
                                       BinHeightMinMax = binHeightMinMax)
        TimeAcceptance.__init__(self, Acceptance = mhe)
