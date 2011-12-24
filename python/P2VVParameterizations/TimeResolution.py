###########################################################################################################################################
## P2VVParameterizations.TimeResolution: Time resolution models                                                                          ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin


class TimeResolution ( _util_parse_mixin ) :
    def __init__( self, model ) : self._model = model
    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def model( self ) : return self._model


class Truth_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel
        self._parseArg( 'time',         kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )

        from ROOT import RooTruthModel as TruthModel
        TimeResolution.__init__( self, ResolutionModel( Name = 'timeResModelTruth'
                                                      , Type = TruthModel
                                                      , Parameters = [ self._time ] )
                               )

class Gaussian_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel
        self._parseArg( 'time',         kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMu',    kwargs, Title = 'Decay time resolution mean',  Value = 0.,  )
        self._parseArg( 'timeResSigma', kwargs, Title = 'Decay time resolution width', Value = 0.05 )

        from ROOT import RooGaussModel as GaussModel
        TimeResolution.__init__( self, ResolutionModel( Name = 'timeResModelGauss'
                                                      , Type = GaussModel
                                                      , Parameters  = [ self._time, self._timeResMu, self._timeResSigma ] )
                               )

class LP2011_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel, AddModel, ConstVar, RealVar
        self._parseArg( 'time',      kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMu', kwargs, Value = -0.0027, Constant = True )
        self._parseArg( 'timeResSF', kwargs, Value = 1.0, MinMax = ( 0.5, 5. ) )

        sigmas = [ ( 3, 0.513  ), ( 2, 0.0853 ), ( 1, 0.0434 ) ]
        fracs  = [ ( 3, 0.0017 ), ( 2, 0.165 ) ]
        self._timeResSigmas = [ ConstVar( Name = 'timeResSigma%s' % num, Value = val ) for num, val in sigmas ]
        self._timeResFracs  = [ ConstVar( Name = 'timeResFrac%s'  % num, Value = val ) for num, val in fracs  ]

        from ROOT import RooGaussModel as GaussModel
        TimeResolution.__init__( self, AddModel( 'timeResModelLP2011'
                                               , [ ResolutionModel( Name = 'timeResLP2011_%s' % numVal[0]
                                                                  , Type = GaussModel
                                                                  , Parameters  = [ self._time, self._timeResMu, sigma, self._timeResSF ]
                                                                  )\
                                                   for ( numVal, sigma ) in zip( sigmas, self._timeResSigmas ) ]
                                               , self._timeResFracs )
                               )

