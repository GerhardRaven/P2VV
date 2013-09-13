from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin
class OtherPdf( _util_parse_mixin ) :
    def __init__( self, **kwargs ) :
        self._pdf = kwargs.pop('pdf')
        for (k,v) in kwargs.iteritems() : setattr( self, '_' + k, v )

    def pdf(self) : return self._pdf

class AmorosoPdf( OtherPdf ) :
    def __init__(self, obs,**kwargs) :
        suffix = kwargs.pop('Suffix', obs.GetName())
        self._a     = self._parseArg('a_%s'      % suffix, kwargs, Title = 'a',      Value = 0,     Constant = True)
        self._theta = self._parseArg('theta_%s'  % suffix, kwargs, Title = 'theta',  Value = 1,     Constant = True)
        self._alpha = self._parseArg('alpha_%s'  % suffix, kwargs, Title = 'alpha',  Value = 25,    MinMax = (0.1, 100))
        self._beta  = self._parseArg('beta_%s'   % suffix, kwargs, Title = 'beta',   Value = 0.3,   MinMax = (0.01,1))
        offset      = self._parseArg('offset_%s' % suffix, kwargs, Title = 'offset', Value = 22200, MinMax = (20000, 30000))
        obs_offset = obs
        if offset.getVal() != 0:
            obs_offset = self._parseArg( '%s_offset' % obs.GetName(), kwargs, Formula = '(@0 - @1) > 0 ? (@0 - @1) : 0.001'
                                        , Arguments = [obs, offset], ObjectType = 'FormulaVar' )

        from P2VV.RooFitWrappers import Pdf
        from ROOT import RooAmoroso as Amoroso

        Name = kwargs.pop('Name', "%s_%s" % (self.__class__.__name__, suffix))
        OtherPdf.__init__(self, pdf = Pdf( Name = Name
                                          , Type = Amoroso
                                          , Parameters = (obs_offset, self._a, self._theta, self._alpha, self._beta)
                                           )
                          )
