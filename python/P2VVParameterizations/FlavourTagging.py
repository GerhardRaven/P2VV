###########################################################################################################################################
## P2VVParameterizations.FlavourTagging: Flavour tagging parameters                                                                      ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin


class TaggingParams ( _util_parse_mixin ):
    def __init__( self, **kwargs ) :
        if 'Dilutions' in kwargs : self._dilutions = kwargs.pop('Dilutions')
        else : raise KeyError('TaggingParams: no key word argument \"Dilutions\" found')

        if 'ADilWTags' in kwargs : self._ADilWTags = kwargs.pop('ADilWTags')
        else : raise KeyError('TaggingParams: no key word argument \"ADilWTags\" found')

        if 'CEvenOdds' in kwargs : self._CEvenOdds = kwargs.pop('CEvenOdds')
        else : raise KeyError('TaggingParams: no key word argument \"CEvenOdds\" found')

    def __getitem__( self, kw ) :
        if kw == 'dilution'                  : return self._dilutions[0]
        if kw == 'ADilWTag'                  : return self._ADilWTags[0]
        if kw in [ 'avgCEven',  'avgCOdd' ]  : return self._CEvenOdds[0][kw]
        if kw in [ 'avgCEvens', 'avgCOdds' ] : return [ CEvenOdd[kw] for CEvenOdd in self._CEvenOdds ]

        return getattr( self, '_' + kw )


class WTagsCoefAsyms_TaggingParams( TaggingParams ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import FormulaVar

        self._parseArg( 'wTag',    kwargs, Title = 'B wrong tag probability',    Value = 0.25, MinMax = ( 0., 0.5 ) )
        self._parseArg( 'wTagBar', kwargs, Title = 'Bbar wrong tag probability', Value = 0.25, MinMax = ( 0., 0.5 ) )

        if 'CEvenOdd' in kwargs :
            CEvenOdd = kwargs.pop('CEvenOdd')
        else :
            from P2VVParameterizations.BBbarAsymmetries import Coefficients_CEvenOdd
            if 'AProd' in kwargs and 'ANorm' in kwargs :
                AProdVal = kwargs.pop('AProd')
                ANormVal = kwargs.pop('ANorm')
                avgCOdd  = ( AProdVal + ANormVal ) / ( 1. + AProdVal * ANormVal )
            elif 'AvgCEven' in kwargs and 'AvgCOdd' in kwargs :
                avgCEven = kwargs.pop('AvgCEven')
                avgCOdd  = kwargs.pop('AvgCOdd') / avgCEven
            CEvenOdd = Coefficients_CEvenOdd( avgCEven = 1., avgCOdd = avgCOdd )

        self._check_extraneous_kw( kwargs )
        TaggingParams.__init__( self
                               , Dilutions = [ FormulaVar(  'tagDilution', '1. - @0 - @1'
                                                          , [ self._wTag, self._wTagBar ], Title = 'Average tagging dilution' ) ]
                               , ADilWTags = [ FormulaVar(  'ADilWTag',    '(@0 - @1) / (1. - @0 - @1)'
                                                          , [ self._wTag, self._wTagBar ], Title = 'Dilution/wrong tag asymmetry' ) ]
                               , CEvenOdds = [ CEvenOdd ]
                              )

