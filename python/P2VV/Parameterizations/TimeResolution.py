###########################################################################################################################################
## P2VVParameterizations.TimeResolution: Time resolution models                                                                          ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin

# initial values for time resolution parameters in fit
timeResSigmaSFVal = 1.45
timeResSigmaSFErr = 0.06

timeResSigmaSFConstrVal = 1.45
timeResSigmaSFConstrErr = 0.06

class TimeResolution ( _util_parse_mixin, _util_extConstraints_mixin, _util_conditionalObs_mixin ):
    def __init__( self, **kwargs ) : 
        if 'Model' not in kwargs : raise KeyError('TimeResolution: please specify a resolution model')
        self._model = kwargs.pop( 'Model' )
        if 'Name'  in kwargs : self._Name  = kwargs.pop('Name')

        if not hasattr( self, '_time' ) :
            from P2VV.RooFitWrappers import RealVar
            self._time = RealVar( Name = self._model.convVar().GetName() )

        params = kwargs.pop( 'Parameters', None )
        if params != None :
            self._params = [ var for var in params ]
            self._model._target_()._parameters = [ var for var in params ]

        # cache integrals as a function of observables
        from ROOT import RooArgSet, RooAbsCategory
        obs = RooArgSet( o for o in self._model.Observables() if o != self._time   )
        if kwargs.pop('Cache', True) and len(obs) :
            print 'P2VV - INFO: invoking %s.parameterizeIntegral(%s)' % ( self._model.GetName(),[o.GetName() for o in obs] )
            self._model.setParameterizeIntegral( obs )
            for o in filter( lambda i : not ( isinstance(i,RooAbsCategory) or i.hasBinning('cache') ), obs ):
                print 'P2VV::TimeResolution: adding cache binning wth 20 bins to %s' % o.GetName()
                o.setBins( 20 , 'cache' )

        _util_conditionalObs_mixin.__init__( self, kwargs )
        self.addConditionals( self._model.ConditionalObservables() )

        _util_extConstraints_mixin.__init__( self, kwargs )
        self.addConstraints( self._model.ExternalConstraints() )

        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def model( self ) : return self._model


class Truth_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        from P2VV.RooFitWrappers import ResolutionModel
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )

        self._check_extraneous_kw( kwargs )
        from ROOT import RooTruthModel as TruthModel
        TimeResolution.__init__( self, Model = ResolutionModel( Name = namePF + 'timeResModelTruth'
                                                               , Type = TruthModel
                                                               , Parameters = [ self._time ]
                                                              )
                               )

class Gaussian_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        namePF = self.getNamePrefix(kwargs)
        scaleBias = kwargs.pop('BiasScaleFactor', True)
        pee = kwargs.pop('PerEventError', False)
        sf_param = kwargs.pop('TimeResSFParam', '')
        self.__split_mean = kwargs.pop('SplitMean', False)
        
        assert(sf_param in ['', 'linear', 'quadratic'])

        extraArgs = {}
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._timeResMu = self._parseArg( 'timeResMu', kwargs, Title = 'Decay time resolution mean', Value = -0.004, MinMax = (-1, 1) )
        if pee :
            self._parseArg( 'sigmat',  kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
            self._timeResSigmaSF = self._parseArg( 'timeResSigmaSF', kwargs, Title = 'Decay time scale factor',
                                                   Value = 1.46, MinMax = (0.1, 2.5) )
            if sf_param:
                self._offset = self._parseArg( 'offset', kwargs, Value = 0.01, Error = 0.001, MinMax = ( 0.000001, 1 ) )
            if sf_param == 'linear':
                linear_var = self._parseArg( 'timeResSigmaSF_linear', kwargs, ObjectType = 'LinearVar', ObsVar = self._sigmat
                                            , Slope = self._timeResSigmaSF, Offset = self._offset )
                params = [self._time, self._timeResMu, linear_var]
            elif sf_param == 'quadratic':
                self._quad = self._parseArg( 'quad', kwargs, Value = -0.01, Error = 0.003, MinMax = ( -0.1, 0.1 ) )
                self._offset.setVal(-0.005)
                self._timeResSigmaSF.setVal(0.0019)
                poly_var = self._parseArg( 'timeResSigmaSF_quad', kwargs, ObsVar = self._sigmat, ObjectType = 'PolyVar'
                                          , Coefficients = [ self._offset, self._timeResSigmaSF, self._quad ] )
                params = [self._time, self._timeResMu, poly_var]
            elif scaleBias :
                params = [self._time, self._timeResMu, self._sigmat, self._timeResSigmaSF, self._timeResSigmaSF]
            else :
                self._timeResMuSF = self._parseArg( 'timeResMuSF', kwargs, Title = 'Decay time bias scale factor', Value = 1
                                                   , Constant = True )
                params = [ self._time, self._timeResMu, self._sigmat, self._timeResMuSF,  self._timeResSigmaSF ]
            extraArgs['ConditionalObservables'] = [self._sigmat]
        else :
            self._timeResSigma = self._parseArg( 'timeResSigma', kwargs, Title = 'Decay time resolution width',
                                                 Value = 0.05,  MinMax = (0.0001, 2.5) )
            params = [ self._time, self._timeResMu, self._timeResSigma ]

        from ROOT import RooGaussModel as GaussModel
        from P2VV.RooFitWrappers import ResolutionModel
        TimeResolution.__init__(  self
                                , Name = namePF + kwargs.pop( 'Name', self.__class__.__name__ )
                                , Model = ResolutionModel(  Name = '%stimeResModelGauss' % namePF
                                                          , Type = GaussModel
                                                          , Parameters  = params
                                                          , **extraArgs)
                                , Cache = kwargs.pop('Cache', True)
                               )
        self._check_extraneous_kw( kwargs )

    def splitVars(self):
        sv = [self._timeResSigmaSF]
        if self.__split_mean:
            sv.append(self._timeResMu)
        return sv

class Multi_Gauss_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        Name = kwargs.pop('Name', 'timeResModelMG')
        namePF = self.getNamePrefix(kwargs)
        from P2VV.RooFitWrappers import ResolutionModel, AddModel
        from ROOT import RooNumber
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ))
        self._parseArg('sigmat', kwargs, Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = (0.0001, 0.12) )
        self._timeResMuSF = self._parseArg( 'timeResMuSF', kwargs, Value = 1.0, Constant = True )
        self._timeResSigmaOffset = self._parseArg( 'timeResSigmaOffset', kwargs, Value = 0.01, Error = 0.001, MinMax = ( 0.00001, 1 ) )
        sigmasSFs = kwargs.pop('ScaleFactors', [(2, 3), (1, 1)])
        gexps = kwargs.pop('GExp', dict([(i[0], False) for i in sigmasSFs]))
        fracs     = kwargs.pop('Fractions', [(2, 0.165)])
        split_fracs = kwargs.pop('SplitFracs', True)
        self.__split_mean = kwargs.pop('SplitMean', False)
        self.__simultaneous = kwargs.pop('Simultaneous', False)
        
        sf_param = kwargs.pop('TimeResSFParam', False)
        from P2VV.RooFitWrappers import RealVar
        self.__placeholder = self._parseArg('tres_placeholder', kwargs, Value = 0,
                                            MinMax = (-1e6, 1e6), Constant = True)
        
        assert(len(sigmasSFs) - 1 == len(fracs))

        ## If this is to be fitted in a non-simultaneous fit, the mean sigmat
        ## value should be passed in to get the parameterisation correct.
        self.__mu_param = kwargs.pop('MeanParameterisation', '')
        assert(not (self.__split_mean and self.__mu_param))
        assert(self.__mu_param in ['quadratic', 'linear', ''])
        self._timeResMu = self._parseArg( 'timeResMu', kwargs, Value = -0.0027, MinMax = ( -2, 2 ) )
        if self.__mu_param == 'linear':
            self._mu = self._timeResMu
            self._mu_offset = self._parseArg( 'timeResMu_offset', kwargs, Value = 0, MinMax = (-2, 2))
            self._mu_slope = self._parseArg( 'timeResMu_slope', kwargs, Value = 0, MinMax = (-2, 2))
            if self.__simultaneous:
                formula = '@0 + @1 * @2'
                args = [self._mu_offset, self._mu_slope, self.__placeholder]
            else:
                formula = '@0 + @1 * (@2 - @3)'
                args = [self._mu_offset, self._mu_slope, self._sigmat, self.__placeholder]
            self._timeResMu = self._parseArg('timeResMu_linear', kwargs, Formula = formula,
                                             ObjectType = 'FormulaVar', Arguments = args)
        elif self.__mu_param == 'quadratic':
            self._mu = self._timeResMu
            self._mu_offset = self._parseArg( 'timeResMu_offset', kwargs, Value = -0.001723, MinMax = (-2, 2))
            self._mu_slope = self._parseArg( 'timeResMu_slope', kwargs, Value = -0.00431, MinMax = (-2, 2))
            self._mu_quad = self._parseArg( 'timeResMu_quad', kwargs, Value = -0.00380, MinMax = (-2, 2))
            if self.__simultaneous:
                formula = '@1 + @2 * @0 + @3 * @0 * @0'
                args = [self.__placeholder, self._mu_offset, self._mu_slope, self._mu_quad]
            else:
                formula = '@2 + @3 * (@0 - @1) + @4 * (@0 - @1) * (@0 - @1)'
                args = [self._sigmat, self.__placeholder, self._mu_offset, self._mu_slope, self._mu_quad]
            self._timeResMu = self._parseArg('timeResMu_quadratic', kwargs, Formula = formula,
                                             ObjectType = 'FormulaVar', Arguments = args)
        self._cache = kwargs.pop('Cache', True)
        pee = kwargs.pop('PerEventError', False)
        if not hasattr(pee, '__iter__'):
            pee = [pee for i in range(len(sigmasSFs))]
        param = kwargs.pop('Parameterise', False)
        assert(param in [False, 'RMS', 'Comb'])
        self._timeResSigmasSFs = [ self._parseArg( 'timeResSigmaSF_%s' % num, kwargs, Value = val, MinMax = (0.001, 20) )\
                                  for num, val in sigmasSFs ]
        self._timeResFracs     = [ self._parseArg( 'timeResFrac%s' % num, kwargs, Value = val, MinMax = (0.0001, 0.99))\
                                  for num, val in fracs ]

        from ROOT import RooNumber
        RooInf = RooNumber.infinity()
        if param == 'RMS': 
            if sf_param:
                self._parseArg( 'sf_mean_offset', kwargs, Value = 1.5, MinMax = (-20, 20) )
                self._parseArg( 'sf_mean_slope', kwargs, Value = -0.1, MinMax = (-20, 20) )
                self._parseArg('sf_sigma_offset', kwargs, Value = 0.4, MinMax = (-20, 20) )
                self._parseArg('sf_sigma_slope', kwargs, Value = -0.05, MinMax = (-20, 20) )
                if self.__simultaneous:
                    formula = '@0 + @1 * @2'
                    args = {'mean'  : [self._sf_mean_offset, self._sf_mean_slope, self.__placeholder],
                            'sigma' : [self._sf_sigma_offset, self._sf_sigma_slope, self.__placeholder]}
                else:
                    formula = '@2 + @3 * (@0 - @1)'
                    args = {'mean'  : [self._sigmat, self.__placeholder, self._sf_mean_offset, self._sf_mean_slope],
                            'sigma' : [self._sigmat, self.__placeholder, self._sf_sigma_offset, self._sf_sigma_slope]}
                self._sf_mean = self._parseArg('timeResSFMean_linear', kwargs,
                                                Formula = formula, ObjectType = 'FormulaVar',
                                                Arguments = args['mean'])
                self._sf_sigma = self._parseArg('timeResSFSigma_linear', kwargs,
                                                Formula = formula, ObjectType = 'FormulaVar',
                                                Arguments = args['sigma'])
            else:
                from math import sqrt
                self._sf_mean = self._parseArg('timeResSFMean', kwargs
                                               , Value = ( ( 1. - fracs[-1][1] ) * sigmasSFs[-1][1] + fracs[-1][1] * sigmasSFs[-2][1] )
                                               , MinMax = (0.5, 5 ))
                self._sf_sigma = self._parseArg( 'timeResSFSigma', kwargs, Value = sqrt((1 - fracs[-1][1]) * sigmasSFs[-1][1] * sigmasSFs[-1][1] +
                                                                                        fracs[-1][1] * sigmasSFs[-2][1] * sigmasSFs[-2][1] - self._sf_mean.getVal() ** 2),
                                                 MinMax = ( 0.01, 2 ))
            self._timeResSigmasSFs[-1] = self._parseArg( Name + '_SF1', kwargs, Formula = '- sqrt(@0 / (1 - @0)) * @1 + @2'
                                                        , Arguments = (self._timeResFracs[-1], self._sf_sigma, self._sf_mean)
                                                        , ObjectType = 'FormulaVar')
            self._timeResSigmasSFs[-2] = self._parseArg( Name + '_SF2', kwargs, Formula = 'sqrt((1 - @0) / @0) * @1 + @2'
                                                        , Arguments = (self._timeResFracs[-1], self._sf_sigma, self._sf_mean)
                                                        , ObjectType = 'FormulaVar')
            if sf_param:
                self._realVars = []
            else:
                self._realVars = [self._sf_mean, self._sf_sigma]
                if split_fracs:
                    self._realVars += self._timeResFracs
        elif param == 'Comb':
            if sf_param:
                self._parseArg( 'sf_mean_slope', kwargs, Value = -3.41081, MinMax = (-20, 20) )
                self._parseArg( 'sf_mean_offset', kwargs, Value = 1.43297, MinMax = (-20, 20) )
                self._parseArg( 'sf2_slope', kwargs, Value = -1, MinMax = (-20, 20) )
                self._parseArg( 'sf2_offset', kwargs, Value = 0, MinMax = (-20, 20) )
                if self.__simultaneous:
                    formula = '@0 + @1 * @2'
                    args = {'mean'  : [self._sf_mean_offset, self._sf_mean_slope, self.__placeholder],
                            'sigma' : [self._sf_sigma_offset, self._sf_sigma_slope, self.__placeholder]}
                else:
                    formula = '@2 + @3 * (@0 - @1)'
                    args = {'mean'  : [self._sigmat, self.__placeholder, self._sf_mean_offset, self._sf_mean_slope],
                            'sf2' : [self._sigmat, self.__placeholder, self._sf2_offset, self._sf2_slope]}
                self._sf2_original = self._timeResSigmasSFs[-2]
                self._timeResSigmasSFs[-2] = self._parseArg( self._sf2_original.GetName() + '_linear', kwargs,
                                                             Formula = formula, ObjectType = 'FormulaVar',
                                                             Arguments = args['sf2'])
                self._sf_mean = self._parseArg( 'timeResSFMean_linear', kwargs,
                                                Formula = formula, ObjectType = 'FormulaVar', 
                                                Arguments = [ self._sf_mean_offset, self._sf_mean_slope, self.__placeholder] )
            else:
                self._sf_mean = self._parseArg('timeResSFMean', kwargs
                                               , Value = ((1. - fracs[-1][1]) * sigmasSFs[-1][1] + fracs[-1][1] * sigmasSFs[-2][1])
                                               , MinMax = (0.5, 5))
            self._timeResSigmasSFs[-1] = self._parseArg(Name + '_SFMean', kwargs, Formula = '(1 / (1 - @0)) * (@1 - @0 * @2)',
                                                        Arguments = ( self._timeResFracs[-1], self._sf_mean, self._timeResSigmasSFs[-2]),
                                                        ObjectType = 'FormulaVar')
            if sf_param:
                self._realVars = []
            else:
                self._realVars = [self._sf_mean, self._timeResSigmasSFs[:-1]]
                if split_fracs:
                    self._realVars += self._timeResFracs
        else:
            self._realVars = [sf for sf in self._timeResSigmasSFs]
            if split_fracs:
                self._realVars += self._timeResFracs

        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        from ROOT import RooGExpModel as GExpModel

        models = []
        for ( numVal, pee ), sigmaSF, in zip( zip(sigmasSFs, pee), self._timeResSigmasSFs):
            gexp = gexps[numVal[0]]
            if not pee:
                params = [ self._time, self._timeResMu, sigmaSF ]
                if gexp:
                    rlife = self._parseArg( 'rlife_%d' % numVal[0], kwargs, Value = 0.1, MinMax = ( 0.0001, 10 ) )
                    params += [rlife]
            else:
                if gexp:
                    rlife = self._parseArg( 'rlife_%d' % numVal[0], kwargs, Value = 0.1, MinMax = ( 0.0001, 10 ) )
                    rlife_sf = self._parseArg( 'rlife_sf', kwargs, Value = 1, ObjectType = 'ConstVar' )
                    self._realVars += [rlife]
                    params = [ self._time, self._timeResMu, self._sigmat, rlife, self._timeResMuSF, sigmaSF, self._timeResMuSF, 'false', 'Normal' ]
                else:
                    params = [ self._time, self._timeResMu, self._sigmat, self._timeResMuSF, sigmaSF ]
                
            model_type = GExpModel if gexp else GaussModel
            model = ResolutionModel(  Name = '%stimeResGauss_%s' % (namePF, numVal[0])
                                      , Type = model_type
                                      , Parameters = params
                                      , ConditionalObservables = [ self._sigmat ] if pee else [])
            models.append(model)
        TimeResolution.__init__(self, Name = Name
                                , Model = AddModel(Name, Models = models, Fractions = self._timeResFracs)
                                , Cache = self._cache)
    def splitVars(self):
        sv = self._realVars[:]
        if self.__split_mean: sv.append(self._timeResMu)
        return sv

    def sigmatPlaceHolder(self):
        return self.__placeholder

class Core_Rest_TimeResolution( TimeResolution ):
    def __init__(self, **kwargs):
        self._core_model = kwargs.pop('CoreModel')
        Name = kwargs.pop('Name', 'timeResModelCR')
        namePF = self.getNamePrefix(kwargs)

        self._sigmat = self._core_model._sigmat
        self._time = self._core_model._time
        ## self._timeResMu = self._parseArg( 'timeResRestMu', kwargs, Value = 0, MinMax = (-1, 1))
        self._timeResMu = self._core_model._timeResMu
        self._timeResMuSF = self._core_model._timeResMuSF

        self._left_sigmaSF = self._parseArg( 'timeResRestLSSF', kwargs, Value = 1.5, MinMax = (0.01, 50))
        self._right_sigmaSF = self._parseArg( 'timeResRestRSSF', kwargs, Value = 1.5, MinMax = (0.01, 50))
        self._left_rlifeSF = self._parseArg( 'timeResRestLTSF', kwargs, Value = 6, MinMax = (0.01, 50))
        self._right_rlifeSF = self._parseArg( 'timeResRestRTSF', kwargs, Value = 2, MinMax = (0.01, 50))

        self._frac_left = self._parseArg('timeResFracLeft', kwargs, Value = 0.2, MinMax = (0.01, 0.99))
        
        self._timeResFracRest = self._parseArg('timeResFracRest', kwargs, Value = 0.05, MinMax = (0.0001, 0.99))

        from P2VV.RooFitWrappers import ResolutionModel, AddModel
        from ROOT import RooGExpModel
        gexps = []
        for side, pars in (('right', [self._right_sigmaSF, self._right_rlifeSF, 'false', 'Flipped']),
                           ('left',  [self._left_sigmaSF, self._left_rlifeSF, 'false', 'Normal'])):
            gexp = ResolutionModel(Name = '%sgexp_%s' % (namePF, side), Type = RooGExpModel,
                                     Parameters = [self._time, self._timeResMu, self._sigmat, self._sigmat, self._timeResMuSF] + pars,
                                     ConditionalObservables = [self._sigmat])
            gexps.append(gexp)
        gexps = AddModel("%stimeResModelRest" % namePF, Models = gexps, Fractions = [self._frac_left])
        TimeResolution.__init__(self, Name = Name,
                                Model = AddModel(Name, Models = [gexps, self._core_model.model()],
                                                 Fractions = [self._timeResFracRest]),
                                Cache = self._core_model._cache)

    def splitVars(self):
        sv = self._core_model.splitVars()
        sv.append(self._sigmaSF)
        return sv

    def sigmatPlaceHolder(self):
        return self._core_model.sigmatPlaceHolder()
        
class Paper2012_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import ResolutionModel, AddModel
        from P2VV.RooFitWrappers import RealVar
        from ROOT import RooNumber
        Name = kwargs.pop('Name', 'timeResModelPaper2012')
        namePF = self.getNamePrefix(kwargs)
        self._parseArg( 'time',           kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMean',    kwargs, Value = 0., Error = 0.1, MinMax = ( -2., 2. ), Constant = True )
        self._parseArg( 'timeResSigma',   kwargs, Title = 'Decay time error', Unit = 'ps', Observable = True, MinMax = ( 0.0, 0.2 ) )
        self._parseArg( 'timeResMeanSF',  kwargs, Value = 1., Constant = True )
        self._parseArg( 'timeResSigmaSF', kwargs, Value = timeResSigmaSFVal, Error = timeResSigmaSFErr, MinMax = ( 0.8, 2.1 ) )
        
        constraints = set()
        timeResMeanConstr = kwargs.pop( 'timeResMeanConstraint', None )
        if type(timeResMeanConstr) == str and timeResMeanConstr == 'fixed' and isinstance( self._timeResMean, RealVar ) :
            self._timeResMean.setConstant(True)

        elif timeResMeanConstr and isinstance( self._timeResMean, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from P2VV.RooFitWrappers import Pdf
            constraints.add( Pdf(  Name = self._timeResMean.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResMean
                                                    , self._parseArg( 'tresMean_constraint_mean', kwargs
                                                                     , Value = self._timeResMean.getVal(), ObjectType = 'ConstVar' )
                                                    , self._parseArg( 'tresMean_constraint_sigma', kwargs
                                                                     , Value = self._timeResMean.getError(), ObjectType = 'ConstVar' )
                                                   ]
                                   )
                             )
        self.__timeResSFConstr = kwargs.pop( 'timeResSFConstraint', None )

        parameters, pc = self.__buildModel()
        constraints |= pc

        models = []
        for i, params in enumerate(parameters):
            from ROOT import RooGaussModel as GaussModel
            model = ResolutionModel(Name = namePF + 'timeResGauss_%d' % (len(parameters) - i)
                                    , Type = GaussModel
                                    , Parameters = params
                                    , ConditionalObservables = [self._timeResSigma]
                                    , ExternalConstraints = constraints)
            models.append(model)

        cache = kwargs.pop( 'Cache', True )
        self._check_extraneous_kw( kwargs )

        if len(models) == 1:
            TimeResolution.__init__(self, Name = Name, Model = models[0], Cache = cache)
        else:
            TimeResolution.__init__(self, Name = Name, Cache = cache
                                    , Model = AddModel(namePF + Name, Models = models, Fractions = [self._timeResSigmaFrac2]))

    def __buildModel(self):
        from P2VV.RooFitWrappers import RealVar
        constraints = set()
        if type(self.__timeResSFConstr) == str and self.__timeResSFConstr == 'fixed' and isinstance( self._timeResSigmaSF, RealVar ) :
            self._timeResSigmaSF.setVal(timeResSigmaSFConstrVal)
            self._timeResSigmaSF.setError(timeResSigmaSFConstrErr)
            self._timeResSigmaSF.setConstant(True)

        elif self.__timeResSFConstr and isinstance( self._timeResSigmaSF, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from P2VV.RooFitWrappers import Pdf
            constraints.add(Pdf(Name = self._timeResSigmaSF.GetName() + '_constraint', Type = Gaussian
                                , Parameters = [self._timeResSigmaSF
                                                , self._parseArg( 'tres_SF_constraint_mean', { }, ObjectType = 'ConstVar'
                                                                  , Value = self._timeResSigmaSF.getVal() )
                                                , self._parseArg( 'tres_SF_constraint_sigma', { }, ObjectType = 'ConstVar'
                                                                  , Value = self._timeResSigmaSF.getError() )]))
        parameters = [self._time, self._timeResMean, self._timeResSigma,
                      self._timeResMeanSF, self._timeResSigmaSF]
        return [parameters], constraints

    
class Gamma_Sigmat( _util_parse_mixin ) :
    def pdf(self) :
        return self._pdf
    def __init__( self, **kwargs ) :
        from ROOT import RooGamma as Gamma
        from P2VV.RooFitWrappers import Pdf
        self._parseArg( 'st',    kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.001,0.5) )
        self._parseArg( 'st_mu', kwargs,   Value = 0, Constant = True )
        self._parseArg( 'st_sig_gamma', kwargs,  MinMax = (5,15) )
        self._parseArg( 'st_sig_beta',  kwargs,  MinMax = (0.0001,0.01) , Value = 0.0025 )
        self._pdf = Pdf( Name = kwargs.pop('Name','sig_st')
                       , Type = Gamma
                       , Parameters = ( self._st, self._st_sig_gamma, self._st_sig_beta, self._st_mu ) 
                       )

class GExp_Gauss(TimeResolution):
    def __init__(self, time, sigmat, **kwargs):
        namePF = self.getNamePrefix(kwargs)
        mean = self._parseArg( 'peak_mean', kwargs, Value = 0, MinMax = ( -0.5, 0.5 ) )
        mean_sf = self._parseArg( 'mean_sf', kwargs, Value = 1, Constant = True )
        gexp_rlife = self._parseArg( 'rlife', kwargs, Value = 0.1, MinMax = ( 0.001, 10 ) )
        gexp_sigma_sf = self._parseArg( 'gexp_sigma_sf', kwargs, Value = 1., MinMax = ( 0.001, 10 ) )
        gexp_rlife_sf = self._parseArg( 'gexp_rlife_sf', kwargs, Value = 1., MinMax = ( 0.001, 10 ), Constant = True )
        gauss_sigma_sf = self._parseArg( 'gauss_sigma_sf', kwargs, Value = 1., MinMax = ( 0.5, 100 ) )
        gexp_frac = self._parseArg( 'gexp_frac', kwargs, Value = 0.2, MinMax = ( 0.0001, 0.9999 ) )

        from ROOT import RooGExpModel
        params = [time, mean, sigmat, gexp_rlife, mean_sf, gexp_sigma_sf, gexp_rlife_sf, 'false', 'Normal']
        gexp_model = ResolutionModel(Name = "%sgexp_model" % namePF, Type = RooGExpModel,
                                     Parameters = params, ConditionalObservables = [sigmat])

        from ROOT import RooGaussModel as GaussModel
        gauss_model = ResolutionModel(Name = "%sgauss_model" % namePF, Type = GaussModel,
                                      ConditionalObservables = [sigmat],
                                      Parameters = [time, mean, sigmat, mean_sf, gauss_sigma_sf])

        TimeResolution.__init__(self, Name = namePF + kwargs.pop( 'Name', self.__class__.__name__ ),
                                Model = AddModel("%sadd_model" % namePF, Fractions = [gexp_frac],
                                                 Models = [gexp_model, gauss_model],
                                                 ConditionalObservables = [sigmat]), Cache = cache)


def rms(sf1, f, sf2):
   return (1 - f) * sf1 * sf1 + f * sf2 * sf2
