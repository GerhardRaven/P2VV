###########################################################################################################################################
## P2VVParameterizations.TimeResolution: Time resolution models                                                                          ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin


class TimeResolution ( _util_parse_mixin ) :
    def __init__( self, model, constraints = None ) : 
        self._model = model
        # constraints should be a type which can be put on rhs of a list.__iadd__ ( += ), i.e.
        # c = list()
        # c += TimeResolution.ExternalConstraints()
        self._constraints = constraints
    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def model( self ) : return self._model
    def ExternalConstraints( self ) : return self._constraints


class Truth_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )

        self._check_extraneous_kw( kwargs )
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

        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        TimeResolution.__init__( self, ResolutionModel( Name = 'timeResModelGauss'
                                                      , Type = GaussModel
                                                      , Parameters  = [ self._time, self._timeResMu, self._timeResSigma ] )
                               )

class LP2011_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel, AddModel, ConstVar, RealVar
        self._parseArg( 'time',      kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMu', kwargs, Value = -0.0027 )
        self._parseArg( 'timeResSF', kwargs, Value = 1.0, MinMax = ( 0.5, 5. ) )

        sigmas = [ ( 3, 0.513  ), ( 2, 0.0853 ), ( 1, 0.0434 ) ]
        fracs  = [ ( 3, 0.0017 ), ( 2, 0.165 ) ]
        self._timeResSigmas = [ ConstVar( Name = 'timeResSigma%s' % num, Value = val ) for num, val in sigmas ]
        self._timeResFracs  = [ ConstVar( Name = 'timeResFrac%s'  % num, Value = val ) for num, val in fracs  ]

        #TODO: add a kw which indicates
        #        1) leave parameters fixed   
        #        2) leave parameters free, and define an external constraint
        #      note that 2 does not imply that one has to actually use the constraint... 
        #      in case 2, it would be nice if one could specify the values of the constraint parameters...
        from ROOT import RooGaussian as Gaussian
        from RooFitWrappers import Pdf
        constraint = Pdf(  Name = self._timeResSF.GetName() + '_constraint', Type = Gaussian
                         , Parameters = [  self._timeResSF
                                         , ConstVar( Name = 'tres_SF_constraint_mean',  Value = 1.00 )
                                         , ConstVar( Name = 'tres_SF_constraint_sigma', Value = 0.04 )
                                        ]
                        )

        Name =  kwargs.pop( 'Name', 'timeResModelLP2011' )
        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        TimeResolution.__init__( self, AddModel( Name
                                               , [ ResolutionModel( Name = 'timeResLP2011_%s' % numVal[0]
                                                                  , Type = GaussModel
                                                                  , Parameters  = [ self._time, self._timeResMu, sigma, self._timeResSF ]
                                                                  )\
                                                   for ( numVal, sigma ) in zip( sigmas, self._timeResSigmas ) ]
                                               , self._timeResFracs )
                                 , constraints = [constraint]
                               )

class Gamma_Sigmat( _util_parse_mixin ) :
    def pdf(self) :
        return self._pdf
    def __init__( self, **kwargs ) :
        from ROOT import RooGamma as Gamma
        from RooFitWrappers import Pdf
        self._parseArg( 'st',    kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.001,0.5) )
        self._parseArg( 'st_mu', kwargs,   Value = 0, Constant = True )
        self._parseArg( 'st_sig_gamma', kwargs,  MinMax = (5,15) )
        self._parseArg( 'st_sig_beta',  kwargs,  MinMax = (0.0001,0.01) , Value = 0.0025 )
        self._pdf = Pdf( Name = kwargs.pop('Name','sig_st')
                       , Type = Gamma
                       , Parameters = ( self._st, self._st_sig_gamma, self._st_sig_beta, self._st_mu ) 
                       )

