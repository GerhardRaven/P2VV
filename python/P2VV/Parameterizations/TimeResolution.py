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
        for obs in self._model.ConditionalObservables() :
            if obs not in self.conditionalObservables() : self.addConditional(obs)

        _util_extConstraints_mixin.__init__( self, kwargs )
        for constr in self._model.ExternalConstraints() :
            if constr not in self.externalConstraints() : self.addConstraint(constr)

        self._check_extraneous_kw( kwargs )

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def model( self ) : return self._model


class Truth_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        from P2VV.RooFitWrappers import ResolutionModel
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
        sf_model = kwargs.pop('TimeResSFModel', '')
        assert(sf_model in ['', 'linear', 'quadratic'])

        extraArgs = {}
        self._parseArg( 'time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ) )
        self._parseArg( 'timeResMu', kwargs, Title = 'Decay time resolution mean', Value = -0.004, MinMax = (-1, 1)  )
        if pee :
            self._parseArg( 'sigmat',  kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
            self._parseArg( 'sigmaSF', kwargs, Title = 'Decay time scale factor',   Value = 1.46,     MinMax = (0.1, 2.5) )
            if sf_model:
                self._offset = self._parseArg('offset', kwargs, Value = 0.01, Error = 0.001,
                                              MinMax = (0.000001, 1))                
            if sf_model == 'linear':
                from P2VV.RooFitWrappers import LinearVar
                linear_var = LinearVar(Name = 'sigmaSF_linear', Observable = self._sigmat,
                                       Slope = self._sigmaSF, Offset = self._offset)
                params = [self._time, self._timeResMu, linear_var]
            elif sf_model == 'quadratic':
                from P2VV.RooFitWrappers import PolyVar
                self._quad = self._parseArg('quad', kwargs, Value = -0.01, Error = 0.003, MinMax = (-0.1, 0.1))
                self._offset.setVal(-0.005)
                self._sigmaSF.setVal(0.0019)
                poly_var = PolyVar(Name = 'sigmaSF_quad', Observable = self._sigmat, Coefficients = [self._offset, self._sigmaSF, self._quad])
                params = [self._time, self._timeResMu, poly_var]
            elif scaleBias :
                params = [self._time, self._timeResMu, self._sigmat, self._sigmaSF, self._sigmaSF]
            else :
                self._parseArg( 'biasSF', kwargs, Title = 'Decay time bias scale factor', Value = 1, Constant = True )
                params = [ self._time, self._timeResMu, self._sigmat, self._biasSF,  self._sigmaSF ]
            extraArgs['ConditionalObservables'] = [self._sigmat]
        else :
            self._parseArg( 'timeResSigma', kwargs, Title = 'Decay time resolution width', Value = 0.05,  MinMax = (0.0001, 2.5) )
            params = [ self._time, self._timeResMu, self._timeResSigma ]

        name = kwargs.pop('Name', 'Gaussian_TimeResolution')
        cache = kwargs.pop('Cache', True)
        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel
        from P2VV.RooFitWrappers import ResolutionModel
        TimeResolution.__init__(  self
                                , Name = name
                                , Model = ResolutionModel(  Name = 'timeResModelGauss'
                                                          , Type = GaussModel
                                                          , Parameters  = params
                                                          , **extraArgs)
                                , Cache = cache
                               )

    def splitVars(self):
        return [self._sigmaSF]

class Multi_Gauss_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        namePF = kwargs.pop( 'ResolutionNamePrefix', '' )

        from P2VV.RooFitWrappers import ResolutionModel, AddModel, ConstVar
        from P2VV.RooFitWrappers import RealVar, FormulaVar, LinearVar
        from ROOT import RooNumber
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0., MinMax = ( -0.5, 5. ))
        self._parseArg('sigmat', kwargs, Title = 'per-event decaytime error', Unit = 'ps', Observable = True, MinMax = (0.0,0.2) )
        self._timeResMu = self._parseArg('%stimeResMu' % namePF, kwargs, Value = -0.0027, MinMax = (-2, 2))
        self._timeResMuSF = self._parseArg('%stimeResMuSF' % namePF, kwargs, Value = 1.0, Constant = True)
        self._timeResSigmaOffset = self._parseArg( '%stimeResSigmaOffset' % namePF, kwargs, Value = 0.01, Error = 0.001, MinMax = ( 0.00001, 1 ) )

        sigmasSFs = kwargs.pop('ScaleFactors', [(2, 3), (1, 1)])
        fracs     = kwargs.pop('Fractions', [(2, 0.165)])
        split_fracs = kwargs.pop('SplitFracs', True)
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

        self._timeResSigmasSFs = [ RealVar( Name = 'timeResSigmaSF_%s' % num, Value = val, MinMax = (0.001, 10) ) for num, val in sigmasSFs ]
        self._timeResFracs  = [ RealVar( Name = 'timeResFrac%s'  % num, Value = val, MinMax = (0.0001, 0.99) ) for num, val in fracs  ]

        use_offset = kwargs.pop('TimeResSFOffset', False)
        assert(not (use_offset and param))
        Name = kwargs.pop('Name', 'timeResModelMG')

        from ROOT import RooNumber
        RooInf = RooNumber.infinity()
        if param == 'RMS': 
            from math import sqrt
            self._rms = RealVar('timeResRMS', Value = sqrt((1 - fracs[0][1]) * sigmasSFs[1][1]
                                                           + fracs[0][1] * sigmasSFs[0][1]),
                                MinMax = (0.8, 5))
            self._timeResSigmasSFs[1] = FormulaVar(Name + '_RMS', 'sqrt(1 / (1 - @0) * (@1 * @1 - @0 * @2 * @2))',
                                                   (self._timeResFracs[0], self._rms, self._timeResSigmasSFs[0]))
            self._realVars = [self._rms, self._timeResSigmasSFs[0]]
            if split_fracs:
                self._realVars += self._timeResFracs
        elif param == 'Comb':
            self._comb = RealVar('timeResComb', Value = ((1 - fracs[0][1]) * sigmasSFs[1][1]
                                                         + fracs[0][1] * sigmasSFs[0][1]),
                                 MinMax = (0.8, 5))
            self._timeResSigmasSFs[1] = FormulaVar(Name + '_Comb', '(1 / (1 - @0)) * (@1 - @0 * @2)',
                                                   (self._timeResFracs[0], self._comb, self._timeResSigmasSFs[0]))
            self._realVars = [self._comb, self._timeResSigmasSFs[0]]
            if split_fracs:
                self._realVars += self._timeResFracs

        elif use_offset:
            sfs = []
            for sf, pe in zip(self._timeResSigmasSFs, pee):
                if pe:
                    sfo = LinearVar(Name = sf.GetName() + '_linear',
                                    Observable = self._sigmat,
                                    Slope = sf, Offset = self._timeResSigmaOffset)
                else:
                    sfo = sf
                sfs.append(sfo)
            self._realVars = [sf for sf in self._timeResSigmasSFs]
            if split_fracs:
                self._realVars += self._timeResFracs
            self._timeResSigmasSFs = sfs
        else:
            self._realVars = [sf for sf in self._timeResSigmasSFs]
            if split_fracs:
                self._realVars += self._timeResFracs

        self._check_extraneous_kw( kwargs )
        from ROOT import RooGaussModel as GaussModel

        models = []
        for ( numVal, pee ), sigmaSF in zip( zip(sigmasSFs, pee), self._timeResSigmasSFs ):
            if use_offset or not pee:
                params = [ self._time, self._timeResMu, sigmaSF ]
            else:
                params = [ self._time, self._timeResMu, self._sigmat, self._timeResMuSF, sigmaSF ]
            model = ResolutionModel(  Name = 'timeResGauss_%s' % numVal[0]
                                      , Type = GaussModel
                                      , Parameters = params
                                      , ConditionalObservables = [ self._sigmat ] if pee else [])
            models.append(model)
        TimeResolution.__init__(self, Name = Name
                                , Model = AddModel(Name, Models = models, Fractions = self._timeResFracs)
                                , Cache = cache)
    def splitVars(self):
        return self._realVars

class Paper2012_TimeResolution ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        sfModel = kwargs.pop( 'timeResSFModel', '' )
        nGauss = kwargs.pop('nGauss', 1)
        assert(sfModel in ['', 'linear', 'quadratic'])
        from P2VV.RooFitWrappers import ResolutionModel, AddModel
        from P2VV.RooFitWrappers import ConstVar, RealVar, FormulaVar
        from ROOT import RooNumber
        Name = kwargs.pop('Name', 'timeResModelPaper2012')
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
        
        constraints = []
        timeResMeanConstr = kwargs.pop( 'timeResMeanConstraint', None )
        if type(timeResMeanConstr) == str and timeResMeanConstr == 'fixed' and isinstance( self._timeResMean, RealVar ) :
            self._timeResMean.setConstant(True)

        elif timeResMeanConstr and isinstance( self._timeResMean, RealVar ) :
            from ROOT import RooGaussian as Gaussian
            from P2VV.RooFitWrappers import Pdf
            constraints.append( Pdf(  Name = self._timeResMean.GetName() + '_constraint', Type = Gaussian
                                    , Parameters = [  self._timeResMean
                                                    , ConstVar( Name = 'tresMean_constraint_mean'
                                                               , Value = self._timeResMean.getVal() )
                                                    , ConstVar( Name = 'tresMean_constraint_sigma'
                                                               , Value = self._timeResMean.getError() )
                                                   ]
                                   )
                             )
        self.__timeResSFConstr = kwargs.pop( 'timeResSFConstraint', None )

        self.__cov = kwargs.pop('Covariance', None)
        if self.__cov: self.__cov = dict([(tuple(sorted(e)), v) for e, v in self.__cov.items()])

        parameters, pc = self.__buildModel(nGauss, sfModel)
        constraints += pc

        models = []
        for i, params in enumerate(parameters):
            from ROOT import RooGaussModel as GaussModel
            model = ResolutionModel(Name = 'timeResGauss_%d' % (len(parameters) - i)
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
                                    , Model = AddModel(Name, Models = models,
                                                       Fractions = [self._timeResSigmaFrac2]))

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
        from P2VV.RooFitWrappers import ConstVar, RealVar, FormulaVar
        sigmaSF1 = FormulaVar('timeResSigmaSF1', '(1 / (1 - @0)) * (@1 - @0 * @2)',
                              (self._timeResSigmaFrac2, self._timeResComb, self._timeResSigmaSF2))
        parameters = [[self._time, self._timeResMean, self._timeResSigma, self._timeResMeanSF, self._timeResSigmaSF2],
                      [self._time, self._timeResMean, self._timeResSigma, self._timeResMeanSF, sigmaSF1]]
        constr_pars = [p for p in [self._timeResComb, self._timeResSigmaFrac2, self._timeResSigmaSF2] if not p.isConstant()]
        constraints = []
        if constr_pars:
            means = [ConstVar(Name = v.GetName() + '_constraint_mean', Value = v.getVal())
                     for v in constr_pars]
            cm = self.__getCovarianceMatrix(constr_pars, self.__cov)
            from P2VV.RooFitWrappers import MultiVarGaussian
            constraints.append(MultiVarGaussian(Name = 'time_resolution_constraint',
                                                Parameters = constr_pars,
                                                CentralValues = means, Covariance = cm))
        return parameters, constraints
        
    def __buildSingleGauss(self, sfModel):
        from P2VV.RooFitWrappers import ConstVar, RealVar, FormulaVar
        constraints = []
        if sfModel == '':
            if type(self.__timeResSFConstr) == str and self.__timeResSFConstr == 'fixed' and isinstance( self._timeResSigmaSF, RealVar ) :
                self._timeResSigmaSF.setVal(timeResSigmaSFConstrVal)
                self._timeResSigmaSF.setError(timeResSigmaSFConstrErr)
                self._timeResSigmaSF.setConstant(True)

            elif self.__timeResSFConstr and isinstance( self._timeResSigmaSF, RealVar ) :
                from ROOT import RooGaussian as Gaussian
                from P2VV.RooFitWrappers import Pdf
                constraints.append(Pdf(Name = self._timeResSigmaSF.GetName() + '_constraint', Type = Gaussian
                                       , Parameters = [self._timeResSigmaSF
                                                       , ConstVar(Name = 'tres_SF_constraint_mean'
                                                                  , Value = self._timeResSigmaSF.getVal() )
                                                       , ConstVar(Name = 'tres_SF_constraint_sigma'
                                                                  , Value = self._timeResSigmaSF.getError() )]))
            parameters = [self._time, self._timeResMean, self._timeResSigma,
                          self._timeResMeanSF, self._timeResSigmaSF]
        elif sfModel == 'linear':
            parameters = [self._timeResSigmaOffset, self._timeResSigmaSF]
            calib_offset = ConstVar(Name = 'calibration_offset', Value = 0.03366)
            self._timeResSigmaLinear = FormulaVar('timeResSigmaLinear',
                                                  '@0 + @1 + @2 * (@3 - @0)',
                                                  [calib_offset] + parameters + [self._timeResSigma])
            parameters = [p for p in parameters if not p.isConstant()]
            if parameters:
                means = [ConstVar(Name = v.GetName() + '_constraint_mean', Value = v.getVal())
                         for v in parameters]
                cm = self.__getCovarianceMatrix(parameters, self.__cov)
                from P2VV.RooFitWrappers import MultiVarGaussian
                constraints.append(MultiVarGaussian(Name = 'time_resolution_constraint',
                                                    Parameters = parameters,
                                                    CentralValues = means, Covariance = cm))
            parameters = [self._time, self._timeResMean, self._timeResSigmaLinear]
        elif sfModel == 'quadratic':
            parameters = [self._timeResSigmaOffset, self._timeResSigmaSF, self._timeResSigmaSF2]
            calib_offset = ConstVar(Name = 'calibration_offset', Value = 0.03366)
            self._timeResSigmaQuad = FormulaVar('timeResSigmaQuad',
                                                '@0 + @1 + @2 * (@4 - @0) + @3 * (@4 - @0)^2',
                                                [calib_offset] + parameters + [self._timeResSigma])
            parameters = [p for p in parameters if not p.isConstant()]
            if parameters:
                means = [ConstVar(Name = v.GetName() + '_constraint_mean', Value = v.getVal())
                         for v in parameters]
                cm = self.__getCovarianceMatrix(parameters, self.__cov)
                from P2VV.RooFitWrappers import MultiVarGaussian
                constraints.append(MultiVarGaussian(Name = 'time_resolution_constraint',
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

