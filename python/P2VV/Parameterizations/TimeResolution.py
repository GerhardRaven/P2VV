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
        namePF = self.getNamePrefix(kwargs)
        from P2VV.RooFitWrappers import ResolutionModel, AddModel
        from ROOT import RooNumber
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ))
        self._parseArg('sigmat', kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
        self._timeResMu = self._parseArg( 'timeResMu', kwargs, Value = -0.0027, MinMax = ( -2, 2 ) )
        self._timeResMuSF = self._parseArg( 'timeResMuSF', kwargs, Value = 1.0, Constant = True )
        self._timeResSigmaOffset = self._parseArg( 'timeResSigmaOffset', kwargs, Value = 0.01, Error = 0.001, MinMax = ( 0.00001, 1 ) )
        sigmasSFs = kwargs.pop('ScaleFactors', [(2, 3), (1, 1)])
        gexps = kwargs.pop('GExp', dict([(i[0], False) for i in sigmasSFs]))
        fracs     = kwargs.pop('Fractions', [(2, 0.165)])
        split_fracs = kwargs.pop('SplitFracs', True)
        self.__split_mean = kwargs.pop('SplitMean', False)

        sf_param = kwargs.pop('TimeResSFParam', False)
        if sf_param:
            self._parseArg('sigmat_mean', kwargs, Value = 3.47648e-02, Constant = True)
            self.__st_placeholder = self._parseArg( self._sigmat.GetName() + '_placeholder', kwargs, Value = 0.033
                                                   , MinMax = (-0.2, 0.2), Constant = True )
        else:
            self.__st_placeholder = None
        
        assert(len(sigmasSFs) - 1 == len(fracs))

        cache = kwargs.pop('Cache', True)
        pee = kwargs.pop('PerEventError', False)
        if not hasattr(pee, '__iter__'):
            pee = [pee for i in range(len(sigmasSFs))]
        param = kwargs.pop('Parameterise', False)
        assert(param in [False, 'RMS', 'Comb'])
        ## Can only reparameterize 2 Gaussians right now
        if param:
            assert(len(sigmasSFs) == 2)

        self._timeResSigmasSFs = [ self._parseArg( 'timeResSigmaSF_%s' % num, kwargs, Value = val, MinMax = (0.001, 20) )\
                                  for num, val in sigmasSFs ]
        self._timeResFracs     = [ self._parseArg( 'timeResFrac%s' % num, kwargs, Value = val, MinMax = (0.0001, 0.99) )\
                                  for num, val in fracs ]

        Name = kwargs.pop('Name', 'timeResModelMG')

        from ROOT import RooNumber
        RooInf = RooNumber.infinity()
        if param == 'RMS': 
            from math import sqrt
            self._rms = self._parseArg( 'timeResRMS', kwargs, Value = sqrt( (1 - fracs[0][1]) * sigmasSFs[1][1]\
                                       + fracs[0][1] * sigmasSFs[0][1]), MinMax = ( 0.8, 5 ) )
            self._timeResSigmasSFs[1] = self._parseArg( Name + '_RMS', kwargs, Formula = 'sqrt(1 / (1 - @0) * (@1 * @1 - @0 * @2 * @2))'
                                                       , Arguments = ( self._timeResFracs[0], self._rms, self._timeResSigmasSFs[0] )
                                                       , ObjectType = 'FormulaVar' )
            self._realVars = [self._rms, self._timeResSigmasSFs[0]]
            if split_fracs:
                self._realVars += self._timeResFracs
        elif param == 'Comb':
            if sf_param:
                self._parseArg( 'sf2_slope', kwargs, Value = -4.08319, MinMax = (-20, 20) )
                self._parseArg( 'sf2_offset', kwargs, Value = 2.03079, MinMax = (-20, 20) )
                self._sf2_original = self._timeResSigmasSFs[0]
                self._timeResSigmasSFs[0] = self._parseArg( self._sf2_original.GetName() + '_linear', kwargs, ObjectType = 'PolyVar'
                                                           , ObsVar = self.__st_placeholder
                                                           , Coefficients = [ self._sf2_offset, self._sf2_slope ] )
                self._parseArg( 'sfc_slope', kwargs, Value = -3.41081, MinMax = (-20, 20) )
                self._parseArg( 'sfc_offset', kwargs, Value = 1.43297, MinMax = (-20, 20) )
                self._comb = self._parseArg( 'timeResComb_linear', kwargs, ObjectType = 'PolyVar', ObsVar = self.__st_placeholder
                                            , Coefficients = [ self._sfc_offset, self._sfc_slope ] )
            else:
                self._comb = self._parseArg(  'timeResComb', kwargs
                                            , Value = ( ( 1. - fracs[0][1] ) * sigmasSFs[1][1] + fracs[0][1] * sigmasSFs[0][1] )
                                            , MinMax = (0.5, 5 )
                                           )
                self._timeResSigmasSFs[1] = self._parseArg( Name + '_Comb', kwargs, Formula = '(1 / (1 - @0)) * (@1 - @0 * @2)'
                                                           , Arguments = ( self._timeResFracs[0], self._comb, self._timeResSigmasSFs[0] )
                                                           , ObjectType = 'FormulaVar' )
            if sf_param:
                self._realVars = []
            else:
                self._realVars = [self._comb, self._timeResSigmasSFs[0]]
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
                    params = [ self._time, self._timeResMu, self._sigmat, rlife, self._timeResMuSF, sigmaSF, sigmaSF, 'false', 'Normal' ]
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
                                , Cache = cache)
    def splitVars(self):
        sv = self._realVars[:]
        if self.__split_mean: sv.append(self._timeResMu)
        return sv

    def sigmatPlaceHolder(self):
        return self.__st_placeholder
        
class Paper2012_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        sfModel = kwargs.pop( 'timeResSFModel', '' )
        nGauss = kwargs.pop('nGauss', 1)
        assert(sfModel in ['', 'linear', 'quadratic'])
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
        if nGauss == 1:
            if sfModel == 'linear':
                self._parseArg( 'timeResSigmaOffset', kwargs, Value = 0.0065, Error = 0.001, MinMax = ( -0.02, 0.1 ) )
            if sfModel == 'quadratic':
                self._parseArg( 'timeResSigmaOffset', kwargs, Value = -0.004, Error = 0.001, MinMax = ( -0.02, 0.1 ) )
                self._parseArg( 'timeResSigmaSF2', kwargs, Value = -0.11, Error = 0.01, MinMax = ( -1, 1 ) )
        elif nGauss == 2:
            self._parseArg('timeResComb', kwargs, Value = 1.45, Error = 4.18e-03, MinMax = (0.8, 5))
            self._parseArg('timeResSigmaFrac2', kwargs, Value = 0.11, Error = 0.01, MinMax = ( -1, 1 ) )
            self._parseArg('timeResSigmaSF2', kwargs, Value = -0.11, Error = 0.01, MinMax = ( -1, 1 ) )
        
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

        self.__cov = kwargs.pop('Covariance', None)
        if self.__cov: self.__cov = dict([(tuple(sorted(e)), v) for e, v in self.__cov.items()])

        parameters, pc = self.__buildModel(nGauss, sfModel)
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

    def __getCovarianceMatrix(self, parameters, covariance):
        from ROOT import TMatrixTSym
        cm = TMatrixTSym('double')(len(parameters))
        from itertools import combinations
        def cwr(pars):
            return [c for c in combinations(pars, 2)] + [(p, p) for p in pars]
        indeces = cwr(range(len(parameters)))
        if covariance:
            keys = dict(zip(indeces, [tuple(sorted(e)) for e in cwr([p.GetName() for p in parameters])]))
        for i, j in indeces:
            if covariance:
                cm[i][j] = covariance[keys[(i, j)]]
                if i != j : cm[j][i] = covariance[keys[(i, j)]]
            else:
                if i == j :
                    cm[i][j] = 1
                else :
                    cm[i][j] = 0
                    cm[j][i] = 0

        return cm

    def __buildModel(self, nGauss, sfModel):
        if nGauss == 1:
            return self.__buildSingleGauss(sfModel)
        elif nGauss == 2:
            return self.__buildDoubleGauss(sfModel)

    def __buildDoubleGauss(self, sfModel):
        sigmaSF1 = self._parseArg( 'timeResSigmaSF1', { }, Formula = '(1 / (1 - @0)) * (@1 - @0 * @2)'
                                  , Arguments = ( self._timeResSigmaFrac2, self._timeResComb, self._timeResSigmaSF2 )
                                  , ObjectType = 'FormulaVar' )
        parameters = [[self._time, self._timeResMean, self._timeResSigma, self._timeResMeanSF, self._timeResSigmaSF2],
                      [self._time, self._timeResMean, self._timeResSigma, self._timeResMeanSF, sigmaSF1]]
        constr_pars = [p for p in [self._timeResComb, self._timeResSigmaFrac2, self._timeResSigmaSF2] if not p.isConstant()]
        constraints = set()
        if constr_pars:
            means = [self._parseArg( v.GetName() + '_constraint_mean', { }, Value = v.getVal(), ObjectType = 'ConstVar' )
                     for v in constr_pars]
            cm = self.__getCovarianceMatrix(constr_pars, self.__cov)
            from P2VV.RooFitWrappers import MultiVarGaussian
            constraints.add(MultiVarGaussian(Name = self.getNamePrefix( { } ) + 'time_resolution_constraint',
                                                Parameters = constr_pars,
                                                CentralValues = means, Covariance = cm))
        return parameters, constraints
        
    def __buildSingleGauss(self, sfModel):
        from P2VV.RooFitWrappers import RealVar
        constraints = set()
        if sfModel == '':
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
        elif sfModel == 'linear':
            parameters = [self._timeResSigmaOffset, self._timeResSigmaSF]
            calib_offset = self._parseArg( 'calibration_offset', { }, Value = 0.03366, ObjectType = 'ConstVar' )
            self._timeResSigmaLinear = self._parseArg( 'timeResSigmaLinear', { }, Formula = '@0 + @1 + @2 * (@3 - @0)'
                                                      , Arguments = [ calib_offset ] + parameters + [ self._timeResSigma ]
                                                      , ObjectType = 'FormulaVar' )
            parameters = [p for p in parameters if not p.isConstant()]
            if parameters:
                means = [self._parseArg( v.GetName() + '_constraint_mean', { }, Value = v.getVal(), ObjectType = 'ConstVar' )
                         for v in parameters]
                cm = self.__getCovarianceMatrix(parameters, self.__cov)
                from P2VV.RooFitWrappers import MultiVarGaussian
                constraints.add(MultiVarGaussian(Name = self.getNamePrefix( { } ) + 'time_resolution_constraint',
                                                    Parameters = parameters,
                                                    CentralValues = means, Covariance = cm))
            parameters = [self._time, self._timeResMean, self._timeResSigmaLinear]
        elif sfModel == 'quadratic':
            parameters = [self._timeResSigmaOffset, self._timeResSigmaSF, self._timeResSigmaSF2]
            calib_offset = self._parseArg( 'calibration_offset', { }, Value = 0.03366, ObjectType = 'ConstVar')
            self._timeResSigmaQuad = self._parseArg( 'timeResSigmaQuad', { }, Formula = '@0 + @1 + @2 * (@4 - @0) + @3 * (@4 - @0)^2'
                                                    , Arguments = [ calib_offset ] + parameters + [ self._timeResSigma ]
                                                    , ObjectType = 'FormulaVar' )
            parameters = [p for p in parameters if not p.isConstant()]
            if parameters:
                means = [self._parseArg( v.GetName() + '_constraint_mean', { }, Value = v.getVal(), ObjectType = 'ConstVar')
                         for v in parameters]
                cm = self.__getCovarianceMatrix(parameters, self.__cov)
                from P2VV.RooFitWrappers import MultiVarGaussian
                constraints.add(MultiVarGaussian(Name = self.getNamePrefix( { } ) + 'time_resolution_constraint',
                                                    Parameters = parameters,
                                                    CentralValues = means, Covariance = cm))
            parameters = [self._time, self._timeResMean, self._timeResSigmaQuad]
        else:
            print 'Error, invalid model for scale factor %s' % sfModel
            assert(False)
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
