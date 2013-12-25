from P2VV.Parameterizations.GeneralUtils import _util_parse_mixin, _util_extConstraints_mixin
from P2VV.Parameterizations.GeneralUtils import valid_combinations, exclusive_combinations
from P2VV.Parameterizations.TimeResolution import TimeResolution



## TODO: move this to RooFitWrappers.CubicSplineFun...
def fitAverageExpToHist(hist,knots, tau) :
    # chisquared fit of spline defined by knots through histogram,
    # do a chisquared based on the average value of the spline in the
    # bin, after multiplying with an exponent with lifetime 1/gamma
    # needs the histogram (and knots) to start at positive values,
    gamma = 1.0/tau
    print 'TimeAcceptacne::fitAverageToHist got gamma: ', gamma
    print 'TimeAcceptance::fitAverageToHist got %s knots : %s'% (len(knots), knots)

    from ROOT import TVectorD, TMatrixD, RooCubicSplineKnot
    x = TVectorD( len(knots) )
    for k in range(len(knots)) : x[k] = knots[k]

    rcsk = RooCubicSplineKnot( x.GetMatrixArray(), x.GetNoElements() )

    coefficients= TVectorD()
    covMatrix=TMatrixD ()
    chisq = rcsk.expIntegral( hist, gamma, coefficients, covMatrix )
    assert chisq>=0 # negative chisq is indication of failure...

    print 'average efficiency fit chisq / dof = %s / %s'% (chisq,hist.GetNbinsX()-coefficients.GetNoElements())
    #coefficients.Print()
    return coefficients



## Since all time acceptances are implemented in the resolution model, we inherit from there
class TimeAcceptance ( TimeResolution ) :
    def __init__( self, **kwargs ) :
        if 'Acceptance' in kwargs: self._acceptance = kwargs.pop( 'Acceptance' )
        else: raise KeyError('P2VV - ERROR: TimeAcceptance.__init__(): TimeAcceptance: please specify an acceptance')
        # Workaround for acceptance normalization in cFit...
        # print 'P2VV - WARNING: TimeAcceptance.__init__(): set "NOCacheAndTrack" for %s' % self._acceptance.GetName()

        from P2VV.RooFitWrappers import BinnedPdf
        if type(self._acceptance.efficiency()) == BinnedPdf:
            # note: constant optimization WILL evaluate RooBinnedPdf as a PDF,
            # and thus normalize it... which in this case we don't want to happen...
            self._acceptance.efficiency().setForceUnitIntegral(True)
            print 'P2VV - WARNING: TimeAcceptance.__init__(): switched setForceUnitIntegral to true for %s' % self._acceptance.GetName()

        TimeResolution.__init__(self, Model = self._acceptance,
                                Conditionals = self._acceptance.ConditionalObservables(),
                                Constraints = self._acceptance.ExternalConstraints(),
                                **kwargs)

    def acceptance( self ) : return self._acceptance

class Moriond2012_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True, MinMax = (0.3, 14))
        input_file = kwargs.pop('Input', 'acceptance.root')
        histogram = kwargs.pop('Histogram', 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins')
        name = kwargs.pop('Name', 'Moriond2012_Acceptance')
        namePF = self.getNamePrefix(kwargs)
        model = kwargs.pop('ResolutionModel')
        parameterization = kwargs.pop('Parameterization','BinnedFun')

        from ROOT import TFile
        with TFile.Open(input_file) as acceptance_file :
            if not acceptance_file: raise ValueError, "Cannot open ROOT file %s" % input_file
            self._hist = acceptance_file.Get(histogram)
            self._hist.SetDirectory(0) # disconnect self._hist from file... otherwise it is deleted when file is closed
        if not self._hist: raise ValueError, 'Cannot get acceptance histogram %s from file' % histogram


        if parameterization == 'Spline' :
            from P2VV.RooFitWrappers import CubicSplineGaussModel, CubicSplineFun

            _hist = self._hist
            nbins = _hist.GetNbinsX()
            knots = [ _hist.GetBinLowEdge(1+i) for i in range( 0, nbins+1)  ]
            knots = knots[0:-1:2]
            rhe = _hist.GetBinLowEdge(nbins)+_hist.GetBinWidth(nbins)
            knots.append(rhe)
            self._coefficients = fitAverageExpToHist( _hist,knots,1.5)

            self._shape = CubicSplineFun(  Name = namePF + name +'_shape'
                                         , Observable = self._time
                                         , Knots = knots
                                         , Coefficients = [ self._parseArg( name + '_shape_%s' % it, kwargs, Value = self._coefficients(it)
                                                                           , ObjectType = 'ConstVar' )\
                                                           for it in range( self._coefficients.GetNoElements() ) ]
                                        )

            if False :
                fr = self._time.frame()
                self._shape.plotOn(fr)
                fr.addTH1( _hist.Clone() )
                x = lambda i : self._shape.u(i)
                from ROOT import TLine
                lines = [ TLine( x(i),0,x(i),0.3 )  for i in range( self._shape.knotSize() ) ]
                for line in lines : fr.addObject( line.Clone() )
                fr.Draw()
                import code
                code.interact(local=locals())

            TimeAcceptance.__init__(self, Acceptance = CubicSplineGaussModel(Name = namePF + name,
                                                                   Efficiency = self._shape,
                                                                   ResolutionModel = model['model'],
                                                                   ConditionalObservables = model.ConditionalObservables(),
                                                                   ExternalConstraints = model.ExternalConstraints()))
        elif parameterization == 'BinnedFun' :
            from P2VV.RooFitWrappers import CubicSplineGaussModel
            self._shape = self._parseArg( name + '_shape', kwargs, ObsVar = self._time, Histogram = self._hist, ObjectType = 'BinnedFun' )
            TimeAcceptance.__init__(self, Acceptance = CubicSplineGaussModel(Name = namePF + name,
                                                                   Efficiency = self._shape,
                                                                   ResolutionModel = model['model'],
                                                                   ConditionalObservables = model.ConditionalObservables(),
                                                                   ExternalConstraints = model.ExternalConstraints()))

        elif parameterization == 'BinnedPdf' :
            from P2VV.RooFitWrappers import BinnedPdf
            self._shape = BinnedPdf( namePF + name + '_shape', Observable = self._time, Histogram = self._hist )
            from P2VV.RooFitWrappers import EffResModel
            TimeAcceptance.__init__(self, Acceptance = EffResModel(Name = namePF + name,
                                                                   Efficiency = self._shape,
                                                                   ResolutionModel = model['model'],
                                                                   ConditionalObservables = model.ConditionalObservables(),
                                                                   ExternalConstraints = model.ExternalConstraints()))

        else :
            raise ValueError, 'Requested unknown Moriond2012_TimeAcceptance Parameterization %s'%parameterization

        self._check_extraneous_kw( kwargs )
        print 'P2VV - INFO: Moriond2012_TimeAcceptance.__init__(): using time efficiency histogram "%s" from file "%s"'\
              % ( histogram, input_file )

class Paper2012_csg_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       MinMax = (0.3, 14))
        self._histograms = kwargs.pop('Histograms')
        input_file = kwargs.pop('Input', 'acceptance.root')
        self.__model = kwargs.pop('ResolutionModel')
        name = kwargs.pop('Name', 'Paper2012_BinnedFunAcceptance')
        fit = kwargs.pop('Fit', False)
        namePF = self.getNamePrefix(kwargs)
        
        with TFile.Open(input_file) as acceptance_file :
            if not acceptance_file:
                raise ValueError, "Cannot open histogram file %s" % input_file
            print 'P2VV - INFO: using time efficiency histograms file "%s"' % input_file
            # transform histograms in map of cat state -> histo name
            # assume only one category for now (compositing could be implemented later)
            for (cat, v) in self._histograms.iteritems():
                for (s, info) in v.iteritems() :
                    hist = info.get('histogram', None)
                    if not hist:
                        continue
                    hist = acceptance_file.Get(info['histogram'])
                    if not hist : raise ValueError, ('Failed to get histrogram %s from file %s' \
                                                     % (info['histogram'], input_file))
                    hist.SetDirectory(0) # disconnect self._hist from file... otherwise it is deleted when file is closed
                    info['histogram'] = hist

        parameterization = kwargs.pop('Parameterization','CubicSplineGaussModel')
        assert parameterization in [ 'CubicSplineGaussModel','EffResModel' ]
        if parameterization == 'CubicSplineGaussModel' :
            from P2VV.RooFitWrappers import CubicSplineGaussModel as AcceptanceModel
        else :
            print 'WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING WARNING WARNING WARNING WARNING'
            print 'WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING WARNING WARNING WARNING WARNING'
            print 'WARNING                                                                               WARNING'
            print 'WARNING EffResModel in combination with split BinnedFun not (yet) fully functional... WARNING'
            print 'WARNING                                                                               WARNING'
            print 'WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING WARNING WARNING WARNING WARNING'
            print 'WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING WARNING WARNING WARNING WARNING'
            from P2VV.RooFitWrappers import EffResModel as AcceptanceModel
        self._shape = self._parseArg(name + '_shape', kwargs, ParNamePrefix = namePF, Fit = fit,
                                     ObsVar = self._time, Histograms = self._histograms, ObjectType = 'BinnedFun')
        acceptance = AcceptanceModel(Name = namePF + name, ParNamePrefix = namePF,
                                     Efficiency = self._shape, ResolutionModel = self.__model['model'],
                                     ConditionalObservables = self.__model.ConditionalObservables() | set( self._histograms.iterkeys()),
                                     ExternalConstraints = self.__model.ExternalConstraints())
        TimeAcceptance.__init__(self, Acceptance = acceptance, Cache = kwargs.pop('Cache', True))
        self._check_extraneous_kw( kwargs )

    def build_constraints(self, original, values):
        # We're fitting and using the average constraint, first
        # create a shape and then the constraint.
        binning = self._shape.base_binning()
        constraints = []
        for (prefix, cat, state), parameters in self._shape.coefficients().iteritems():
            if prefix:
                shape_name = '_'.join([prefix, cat, state, 'shape'])
            else:
                shape_name = '_'.join([cat, state, 'shape'])
            from P2VV.RooFitWrappers import BinnedPdf
            shape = BinnedPdf(Name = shape_name, Observable = self._time, Coefficients = parameters,
                              Binning = binning if type(binning) == str else binning.GetName())
            shape.setForceUnitIntegral(True) # insure BinnedPdf is used as function, not PDF, i.e. skip normalization!

            from P2VV.RooFitWrappers import EffResModel
            res_model = self.__model.model()
            eff_model = EffResModel(Name = shape_name + '_effres', Efficiency = shape,
                                    ResolutionModel = res_model)

            # Set all observables constant for the shape to work around a
            # RooFit limitation
            from ROOT import RooArgSet
            obs_set = RooArgSet(res_model.Observables().intersection(res_model.ConditionalObservables()))
            observables = eff_model.getObservables(obs_set)
            for o in observables: o.setConstant(True)

            av_name = '_'.join((prefix, cat, state, 'average'))
            from P2VV.RooFitWrappers import RealVar, Pdf
            mean  = RealVar(Name = av_name + '_constraint_mean', Value = values[(cat, state)][0], Constant = True)
            sigma = RealVar(Name = av_name + '_constraint_sigma', Value = values[(cat, state)][1], Constant = True)
            constraints.append(Pdf(Name = av_name + '_constraint', Type = 'RooAvEffConstraint',
                                   Parameters = [original, eff_model, mean, sigma]))
        return constraints

    def shapes(self):
        return [self._shape]

class Paper2012_mer_TimeAcceptance(TimeAcceptance):
    def __init__(self, **kwargs ) :
        from ROOT import TFile
        self._parseArg('time', kwargs, Title = 'Decay time', Unit = 'ps', Observable = True,
                       MinMax = (0.3, 14))
        histograms = kwargs.pop('Histograms')
        original = kwargs.pop('Original')
        input_file = kwargs.pop('Input', 'acceptance.root')
        acceptance_file = TFile.Open(input_file)
        fit = kwargs.pop('Fit')
        model = kwargs.pop('ResolutionModel')
        binHeightMinMax = kwargs.pop('BinHeightMinMax', (0.001, 0.999))
        spline = kwargs.pop('Spline', False)
        smooth = kwargs.pop('SmoothSpline', 0.1)

        if not acceptance_file:
            raise ValueError, "Cannot open histogram file %s" % input_file
        print 'P2VV - INFO: using time efficiency histograms file "%s"' % input_file

        acceptance_name = kwargs.pop('Name', 'Paper2012_FitAcceptance')
        namePF = self.getNamePrefix(kwargs)

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
        mhe = MultiHistEfficiencyModel( Name = namePF + acceptance_name, Observable = self._time,
                                       Bins = bin_spec, Relative = rel_spec,
                                       ConditionalCategories = True,
                                       FitAcceptance = fit, UseSingleBinConstraint = False,
                                       ResolutionModel = model['model'], Original = original,
                                       ConditionalObservables = model.ConditionalObservables(),
                                       ExternalConstraints = model.ExternalConstraints(),
                                       BinHeightMinMax = binHeightMinMax,
                                       Spline = spline, SmoothSpline = smooth)
        TimeAcceptance.__init__(self, Acceptance = mhe)
        self._check_extraneous_kw( kwargs )
