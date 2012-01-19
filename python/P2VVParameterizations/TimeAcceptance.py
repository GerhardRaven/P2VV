
from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin

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
        TimeAcceptance.__init__( self, Acceptance = HistFunc('time_acceptance', Histogram = self._hist, Observables = [self._time]))
