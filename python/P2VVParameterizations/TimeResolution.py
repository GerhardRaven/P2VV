###########################################################################################################################################
## P2VVParameterizations.TimeResolution: Time resolution models                                                                          ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VVParameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin

# initial values for time resolution parameters in fit
timeResSigmaSFVal = 1.4
timeResSigmaSFErr = 0.06

timeResSigmaSFConstrVal = 1.45
timeResSigmaSFConstrErr = 0.06


class TimeResolution ( _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin ):
    def __init__( self, **kwargs ) : 
        if 'Model' in kwargs : self._model = kwargs.pop( 'Model' )
        else :                 raise KeyError('TimeResolution: please specify a resolution model')
        if 'Name' in kwargs: self._Name = kwargs.pop('Name')

        # cache integrals as a function of observables
        cache = kwargs.pop('Cache', True)
        from ROOT import RooAbsReal, RooArgSet
        realObs = RooArgSet( [ o._var for o in self._model.Observables() if isinstance(o._var,RooAbsReal) and o != self._time  ]  )
        if cache and len(realObs) :
            print 'invoking %s.parameterizeIntegral(%s)' % ( self._model.GetName(),[o.GetName() for o in realObs] )
            self._model.setParameterizeIntegral( realObs )
            for o in realObs :
                if not o.hasBinning('cache') : 
                    print 'adding cache binning to %s' % o.GetName()
                    o.setBins( 20 , 'cache' )


        _util_conditionalObs_mixin.__init__( self, kwargs )
        _util_extConstraints_mixin.__init__( self, kwargs )
        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def model( self ) : return self._model


class Truth_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )

        self._check_extraneous_kw( kwargs )
        from ROOT import RooTruthModel as TruthModel
        TimeResolution.__init__( self, Model = ResolutionModel( Name = 'timeResModelTruth'
                                                               , Type = TruthModel
                                                               , Parameters = [ self._time ]
                                                              )
                               )

class Gaussian_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        scaleBias = kwargs.pop('BiasScaleFactor', True)
        pee = kwargs.pop('PerEventError', False)
        use_offset = kwargs.pop('TimeResSFOffset', False)
        extraArgs = {}
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        if pee :
            self._parseArg( 'sigmat',  kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
            self._parseArg( 'bias',    kwargs, Title = 'Decay time biased',         Value = -0.17,    MinMax = (-1, 1)  )
            self._parseArg( 'sigmaSF', kwargs, Title = 'Decay time scale factor',   Value = 1.46,     MinMax = (0.1, 2.5) )
            if use_offset:
                self._offset = self._parseArg('offset', kwargs, Value = 0.01, Error = 0.001,
                                              MinMax = (0.000001, 1))
                from RooFitWrappers import LinearVar
                linear_var = LinearVar(Name = 'sigmaSF_linear', Observable = self._sigmat,
                                       Slope = self._sigmaSF, Offset = self._offset)
                params = [ self._time, self._bias, linear_var ]
            elif scaleBias :
                params = [ self._time, self._bias, self._sigmat, self._sigmaSF, self._sigmaSF ]
            else :
                self._parseArg( 'biasSF', kwargs, Title = 'Decay time bias scale factor', Value = 1, Constant = True )
                params = [ self._time, self._bias, self._sigmat, self._biasSF,  self._sigmaSF ]
            extraArgs['ConditionalObservables'] = [ self._sigmat ]

        else :
            self._parseArg( 'timeResMu', kwargs, Title = 'Decay time resolution mean', Value = -0.00017, MinMax = (-1, 1)  )
            self._parseArg( 'timeResSigma', kwargs, Title = 'Decay time resolution width', Value = 0.05,  MinMax = (0.0001, 2.5) )
            params = [ self._time, self._timeResMu, self._timeResSigma ]

        name = kwargs.pop('Name', 'Gaussian_TimeResolution')
        cache = kwargs.pop('Cache', True)
        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        from RooFitWrappers import ResolutionModel
        TimeResolution.__init__(  self
                                , Name = name
                                , Model = ResolutionModel(  Name = 'timeResModelGauss'
                                                          , Type = GaussModel
                                                          , Parameters  = params
                                                          , **extraArgs)
                                , Cache = cache
                               )

class LP2011_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        namePF = kwargs.pop( 'ResolutionNamePrefix', '' )

        from RooFitWrappers import ResolutionModel, AddModel, ConstVar, RealVar
        from ROOT import RooNumber
        self._parseArg('time',      kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ))
        self._timeResMu = self._parseArg('%stimeResMu' % namePF, kwargs, Value = -0.0027)
        self._timeResSF = self._parseArg('%stimeResSF' % namePF, kwargs, Value = 1.0, Error = 0.04, MinMax = ( 0.1, 5. ) )

        sigmas = [ ( 3, 0.513  ), ( 2, 0.0853 ), ( 1, 0.0434 ) ]
        fracs  = [ ( 3, 0.0017 ), ( 2, 0.165 ) ]
        self._timeResSigmas = [ ConstVar( Name = 'timeResSigma%s' % num, Value = val ) for num, val in sigmas ]
        self._timeResFracs  = [ ConstVar( Name = 'timeResFrac%s'  % num, Value = val ) for num, val in fracs  ]

        constraints = []
        timeResSFConstr = kwargs.pop( 'timeResSFConstraint', None )
        if type(timeResSFConstr) == str and timeResSFConstr == 'fixed' :
            self._timeResSF.setConstant(True)
        elif timeResSFConstr :
            from ROOT import RooGaussian as Gaussian
            from RooFitWrappers import Pdf
            constraints.append( Pdf(  Name = self._timeResSF.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResSF
                                                    , ConstVar(  Name = '%stres_SF_constraint_mean' % namePF
                                                               , Value = self._timeResSF.getVal()
                                                              )
                                                    , ConstVar(  Name = '%stres_SF_constraint_sigma' % namePF
                                                               , Value = self._timeResSF.getError()
                                                              )
                                                   ]
                                   )
                             )

        self._check_extraneous_kw( kwargs )
        Name = kwargs.pop('Name', 'timeResModelLP2011')
        from ROOT import RooGaussModel as GaussModel
        models = [ ResolutionModel(  Name = 'timeResLP2011_%s' % numVal[0]
                                     , Type = GaussModel
                                     , Parameters = [self._time, self._timeResMu, sigma, self._timeResSF]
                                     , ExternalConstraints = constraints
                                     )\
                   for ( numVal, sigma ) in zip( sigmas, self._timeResSigmas )
                   ]
        TimeResolution.__init__(self, Name = Name
                                , Model = AddModel(Name, Models = models, Fractions = self._timeResFracs)
                                , Constraints = constraints
                                )

class Multi_Gauss_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        namePF = kwargs.pop( 'ResolutionNamePrefix', '' )

        from RooFitWrappers import ResolutionModel, AddModel, ConstVar
        from RooFitWrappers import RealVar, FormulaVar, LinearVar
        from ROOT import RooNumber
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ))
        self._parseArg('sigmat', kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
        self._timeResMu = self._parseArg('%stimeResMu' % namePF, kwargs, Value = -0.0027, MinMax = (-2, 2))
        self._timeResMuSF = self._parseArg('%stimeResMuSF' % namePF, kwargs, Value = 1.0, Constant = True)
        self._timeResSigmaOffset = self._parseArg( '%stimeResSigmaOffset' % namePF, kwargs, Value = 0.01, Error = 0.001, MinMax = ( 0.00001, 1 ) )

        sigmasSFs = kwargs.pop('ScaleFactors', [(2, 3), (1, 1)])
        fracs     = kwargs.pop('Fractions', [(2, 0.165)])
        assert(len(sigmasSFs) - 1 == len(fracs))

        cache = kwargs.pop('Cache', True)
        pee = kwargs.pop('PerEventError', False)
        param_rms = kwargs.pop('ParamRMS', False)
        ## Can only reparameterize 2 Gaussians right now
        if param_rms:
            assert(len(sigmasSFs) == 2)

        self._timeResSigmasSFs = [ RealVar( Name = 'timeResSigmaSF_%s' % num, Value = val, MinMax = (0.001, 20) ) for num, val in sigmasSFs ]
        self._timeResFracs  = [ RealVar( Name = 'timeResFrac%s'  % num, Value = val, MinMax = (0.001, 0.99) ) for num, val in fracs  ]

        use_offset = kwargs.pop('TimeResSFOffset', False)
        assert(not (use_offset and param_rms))
        Name = kwargs.pop('Name', 'timeResModelMG')
        if param_rms: 
            from ROOT import RooNumber
            RooInf = RooNumber.infinity()
            from math import sqrt
            self._rms = RealVar('timeResRMS', Value = sqrt((1 - fracs[0][1]) * sigmasSFs[1][1]
                                                           + fracs[0][1] * sigmasSFs[0][1]),
                                MinMax = (0, RooInf))
            self._timeResSigmasSFs[1] = FormulaVar(Name + '_RMS', 'sqrt(1 / (1 - @0) * (@1 * @1 - @0 * @2 * @2))',
                                                   (self._timeResFracs[0], self._rms, self._timeResSigmasSFs[0]))
        elif use_offset:
            assert(pee)
            self._timeResSigmasSFs = [LinearVar(Name = sf.GetName() + '_linear',
                                                Observable = self._sigmat,
                                                Slope = sf, Offset = self._timeResSigmaOffset)
                                      for sf in self._timeResSigmasSFs]

        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel

        models = []
        for ( numVal, sigmaSF ) in zip( sigmasSFs, self._timeResSigmasSFs ):
            if use_offset:
                params = [ self._time, self._timeResMu, sigmaSF ]
            elif pee:
                params = [ self._time, self._timeResMu, self._sigmat, self._timeResMuSF, sigmaSF ]
            else:
                params = [ self._time, self._timeResMu, sigmaSF ]
            model = ResolutionModel(  Name = 'timeResGauss_%s' % numVal[0]
                                      , Type = GaussModel
                                      , Parameters = params
                                      , ConditionalObservables = [ self._sigmat ] if pee else []
                                      )
            models.append(model)
        TimeResolution.__init__(self, Name = Name
                                , Model = AddModel(Name, Models = models, Fractions = self._timeResFracs)
                                , Cache = cache)

class Moriond2012_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from RooFitWrappers import ResolutionModel, AddModel, ConstVar, RealVar
        from ROOT import RooNumber
        self._parseArg( 'time',           kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMean',    kwargs, Value = 0., ObjectType = 'ConstVar' )
        self._parseArg( 'timeResSigma',   kwargs, Title = 'Decay time error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
        self._parseArg( 'timeResMeanSF',  kwargs, Value = 1., ObjectType = 'ConstVar' )
        self._parseArg( 'timeResSigmaSF', kwargs, Value = 1.45, Error = 0.06, MinMax = ( 0.1, 5. ) )

        constraints = []
        timeResSFConstr = kwargs.pop( 'timeResSFConstraint', None )
        if type(timeResSFConstr) == str and timeResSFConstr == 'fixed' and isinstance( self._timeResSigmaSF, RealVar ) :
            self._timeResSigmaSF.setConstant(True)

        elif timeResSFConstr and isinstance( self._timeResSigmaSF, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from RooFitWrappers import Pdf
            constraints.append( Pdf(  Name = self._timeResSigmaSF.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResSigmaSF
                                                    , ConstVar( Name = 'tres_SF_constraint_mean'
                                                               ,  Value = self._timeResSigmaSF.getVal() )
                                                    , ConstVar( Name = 'tres_SF_constraint_sigma'
                                                               , Value = self._timeResSigmaSF.getError() )
                                                   ]
                                   )
                              )

        Name =  kwargs.pop('Name', 'timeResModelMoriond2012')
        cache = kwargs.pop('Cache', True)
        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        TimeResolution.__init__(  self
                                , Model =  ResolutionModel(  Name = Name
                                                           , Type = GaussModel
                                                           , Parameters = [  self._time
                                                                           , self._timeResMean, self._timeResSigma
                                                                           , self._timeResMeanSF, self._timeResSigmaSF
                                                                          ]
                                                           , ConditionalObservables = [ self._timeResSigma ]
                                                           , ExternalConstraints = constraints
                                                          )
                                , Conditional = self._timeResSigma
                                , Constraints = constraints
                                , Cache = cache
                               )

        from ROOT import RooArgSet

class Paper2012_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        useOffset = kwargs.pop( 'timeResSFOffset', False )

        from RooFitWrappers import ResolutionModel, AddModel, ConstVar, RealVar, LinearVar
        from ROOT import RooNumber
        self._parseArg( 'time',           kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMean',    kwargs, Value = 0., Error = 0.1, MinMax = ( -2., 2. ), Constant = True )
        self._parseArg( 'timeResSigma',   kwargs, Title = 'Decay time error', Unit = 'ps', Observable = True, MinMax = ( 0.0, 0.2 ) )
        self._parseArg( 'timeResMeanSF',  kwargs, timeResMeanSF = self._timeResSigma )
        self._parseArg( 'timeResSigmaSF', kwargs, Value = timeResSigmaSFVal, Error = timeResSigmaSFErr, MinMax = ( 0.8, 2.1 ) )
        if useOffset :
            self._parseArg( 'timeResSigmaOffset', kwargs, Value = 0.0065, Error = 0.001, MinMax = ( 0.00001, 0.1 ) )

        constraints = []
        timeResMeanConstr = kwargs.pop( 'timeResMeanConstraint', None )
        if type(timeResMeanConstr) == str and timeResMeanConstr == 'fixed' and isinstance( self._timeResMean, RealVar ) :
            self._timeResMean.setConstant(True)

        elif timeResMeanConstr and isinstance( self._timeResMean, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from RooFitWrappers import Pdf
            constraints.append( Pdf(  Name = self._timeResMean.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResMean
                                                    , ConstVar( Name = 'tresMean_constraint_mean'
                                                               , Value = self._timeResMean.getVal() )
                                                    , ConstVar( Name = 'tresMean_constraint_sigma'
                                                               , Value = self._timeResMean.getError() )
                                                   ]
                                   )
                             )

        timeResSFConstr = kwargs.pop( 'timeResSFConstraint', None )
        if type(timeResSFConstr) == str and timeResSFConstr == 'fixed' and isinstance( self._timeResSigmaSF, RealVar ) :
            self._timeResSigmaSF.setVal(timeResSigmaSFConstrVal)
            self._timeResSigmaSF.setError(timeResSigmaSFConstrErr)
            self._timeResSigmaSF.setConstant(True)

        elif timeResSFConstr and isinstance( self._timeResSigmaSF, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from RooFitWrappers import Pdf
            constraints.append( Pdf(  Name = self._timeResSigmaSF.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResSigmaSF
                                                    , ConstVar( Name = 'tres_SF_constraint_mean'
                                                               ,  Value = timeResSigmaSFConstrVal )
                                                    , ConstVar( Name = 'tres_SF_constraint_sigma'
                                                               , Value = timeResSigmaSFConstrErr )
                                                   ]
                                   )
                              )

        if useOffset:
            self._timeResSigmaLinear = LinearVar(Name = 'timeResSigmaLinear',
                                                 Observable = self._timeResSigma,
                                                 Slope = self._timeResSigmaSF,
                                                 Offset = self._timeResSigmaOffset)
            constraints.append(Pdf(Name = self._timeResSigmaOffset.GetName() + '_constraint', Type = Gaussian
                                   , Parameters = [self._timeResSigmaOffset
                                                   , ConstVar(Name = 'tres_offset_constraint_mean'
                                                              , Value = self._timeResSigmaOffset.getVal())
                                                   , ConstVar(Name = 'tres_offset_constraint_sigma'
                                                              , Value = self._timeResSigmaOffset.getError())
                                                   ]
                                   )
                               )
        
        Name =  kwargs.pop( 'Name', 'timeResModelPaper2012' )
        cache = kwargs.pop( 'Cache', True )
        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        if useOffset:
            parameters = [self._time, self._timeResMean, self._timeResSigmaLinear]
        else:
            parameters = [self._time, self._timeResMean, self._timeResSigma,
                          self._timeResMeanSF, self._timeResSigmaSF]
        TimeResolution.__init__(  self
                                , Model =  ResolutionModel(  Name = Name
                                                           , Type = GaussModel
                                                           , Parameters = parameters
                                                           , ConditionalObservables = [ self._timeResSigma ]
                                                           , ExternalConstraints = constraints
                                                          )
                                , Conditional = self._timeResSigma
                                , Constraints = constraints
                                , Cache = cache
                               )

        from ROOT import RooArgSet

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

