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

class BinCounter(object):
    def __init__(self, data, time, binning, weightSquared = False):
        self.__data = data
        self.__time = data.get().find(time.GetName())
        self.__binning = binning
        self.__wSq = weightSquared

        from collections import defaultdict
        self.__levels = defaultdict(dict)
        self.__levelsSq = defaultdict(dict)
        self.__categories = defaultdict(set)
        
    def add_bins(self, level, cat_def):
        d = sorted(['{0} == {0}::{1}'.format(c, s) for c, s in cat_def.iteritems()])
        self.__levels[level][tuple(d)] = self.__binning.numBins() * [0.]
        self.__levelsSq[level][tuple(d)] = self.__binning.numBins() * [0.]
        self.__categories[level] |= set(cat_def.iterkeys())
        
    def run(self):
        obs = self.__data.get()
        categories = dict([(c, obs.find(c)) for cats in self.__categories.itervalues() for c in cats])
        assert(all(categories.itervalues()))
        for i in range(self.__data.numEntries()):
            self.__data.get(i)
            for level, cats in self.__categories.iteritems():
                k = tuple(sorted(['{0} == {0}::{1}'.format(c, categories[c].getLabel()) for c in cats]))
                b = self.__binning.binNumber(self.__time.getVal())
                if k not in self.__levels[level]:
                    continue
                self.__levels[level][k][b] += self.__data.weight()
                if self.__wSq :
                    self.__levelsSq[level][k][b] += self.__data.weight()**2

    def get_bins(self, level, cat_def):
        k = tuple(sorted(['{0} == {0}::{1}'.format(c, s) for c, s in cat_def.iteritems()]))
        return self.__levels[level][k]

    def get_binsSq(self, level, cat_def):
        k = tuple(sorted(['{0} == {0}::{1}'.format(c, s) for c, s in cat_def.iteritems()]))
        return self.__levelsSq[level][k]

    def levels(self):
        return self.__levels

    def levelsSq(self):
        return self.__levelsSq

    def categories(self):
        return self.__categories
    
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
        fit = kwargs.pop('Fit', False)

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
            self._coefficients = fitAverageExpToHist( _hist,knots, 1. / 0.66091 )

            if fit : coefArgs = dict( ObjectType = 'RealVar', MinMax = ( -50, 50. ) )
            else :   coefArgs = dict( ObjectType = 'ConstVar' )
            coefVars = [ self._parseArg( name + '_shape_%03d' % it, kwargs, Value = self._coefficients(it), **coefArgs )\
                        for it in range( self._coefficients.GetNoElements() ) ]
            if fit : coefVars[0].setConstant(True)
            self._shape = CubicSplineFun( Name = namePF + name +'_shape', Observable = self._time, Knots = knots, Coefficients = coefVars )

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
        rbo = kwargs.pop('RandomBinOrder', False)
        self.__bc = None
        
        print 'P2VV - INFO: Paper2012_csg_TimeAcceptance: using time efficiency histograms file "%s"' % input_file
        # transform histograms in map of cat state -> histo name
        # assume only one category for now (compositing could be implemented later)
        for (cat, v) in self._histograms.iteritems():
            for (s, info) in v.iteritems() :
                hist = info.get('histogram', None)
                if not hist:
                    continue
                with TFile.Open(input_file) as acceptance_file:
                    if not acceptance_file:
                        raise ValueError, "Cannot open histogram file %s" % input_file
                    hist = acceptance_file.Get(info['histogram'])
                    if not hist : raise ValueError, ('Failed to get histogram %s from file %s' \
                                                     % (info['histogram'], input_file))
                    hist.SetDirectory(0) # disconnect self._hist from file... otherwise it is deleted when file is closed
                    info['histogram'] = hist
                    print 'P2VV - INFO: Paper2012_csg_TimeAcceptance: time efficiency category "%s/%s": using histogram "%s"'\
                           % ( cat.GetName(), s, hist.GetName() )

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
                                     ObsVar = self._time, Histograms = self._histograms,
                                     RandomBinOrder = rbo, ObjectType = 'BinnedFun')
        acceptance = AcceptanceModel(Name = namePF + name, ParNamePrefix = namePF,
                                     Efficiency = self._shape, ResolutionModel = self.__model['model'],
                                     ConditionalObservables = self.__model.ConditionalObservables() | set( self._histograms.iterkeys()),
                                     ExternalConstraints = self.__model.ExternalConstraints())
        TimeAcceptance.__init__(self, Acceptance = acceptance, Cache = kwargs.pop('Cache', True))
        self._check_extraneous_kw( kwargs )

    def build_av_constraints(self, original, values):
        # We're fitting and using the average constraint, first
        # create a shape and then the constraint.
        binning = self._shape.base_binning()
        constraints = set()
        for (prefix, cat, state), parameters in self._shape.coefficients().iteritems():
            if len(parameters) == 1:
                parameters[0].setConstant(True)
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
            constraints.add(Pdf(Name = av_name + '_constraint', Type = 'RooAvEffConstraint',
                                Parameters = [original, eff_model, mean, sigma]))
        return constraints

    def build_multinomial_constraints(self, data, observables, extra_cat = None):
        from P2VV.RooFitWrappers import RooObject
        from collections import defaultdict
        from itertools import product

        ## We're fitting and using the multinomial constraint, first count the
        ## number of events of a types N_a, N_b and N_ab in data and use them to
        ## set the values of the parameters close to the minimum. Then construct
        ## the constraint.
        binning = self._shape.base_binning()
        constraints = set()

        ## Reorganise the coefficients returned by the shape for convenience. If
        ## no prefix has been defined and an extra category to split on has been
        ## given, replace the prefixes.
        epsilons = defaultdict(dict)
        coefficients = self._shape.coefficients()
        if extra_cat and all((not k[0] for k in coefficients.iterkeys())):
            states = [(c.GetName(), s.GetName()) for c, s in product([extra_cat], extra_cat)]
            coefficients = dict([(('_'.join(prefix), cat, state), v) for (prefix, ((p, cat, state), v)) in product(states, coefficients.iteritems())])
        for (prefix, cat, state), parameters in coefficients.iteritems():
            level = cat[ : 4]
            ## Check if a constraint for this combination already exists.
            constraint_name = '_'.join([e for e in (prefix, level, 'multinomial') if e])
            if constraint_name in RooObject._ws._rooobjects:
                continue
            epsilons[(prefix, level)][(cat, state)] = parameters

        if not(epsilons):
            return constraints
                
        obs = dict([(o.GetName(), o) for o in observables.itervalues()])
        time = obs['time']

        ## Count the number of events of types N_a, N_b and N_ab in bins in the
        ## data.
        from P2VV.Parameterizations.GeneralUtils import valid_combinations
        self.__bc = BinCounter(data, time, binning)
        for (prefix, level), parameters in epsilons.iteritems():
            valid = valid_combinations([[(obs[c], s) for c, s in parameters.iterkeys()]])
            for e in valid:
                d = {}
                if prefix and '_' in prefix and prefix.split('_')[0] in observables:
                    d.update(dict([prefix.split('_')]))
                for c, l in e:
                    d[c.GetName()] = l
                self.__bc.add_bins(level, d)
        self.__bc.run()

        ## Build the constraint
        from ROOT import RooArgList
        for (prefix, level), parameters in epsilons.iteritems():
            args = {'Epsilons' : [], 'N' : []}
            par_keys = set(parameters.iterkeys())
            ## Calculate the valid combinations for this level, 2 for HLT1 like
            ## 3 for HLT2 like
            valid = valid_combinations([[(obs[c], s) for c, s in par_keys]])
            for e in valid:
                ## Add the cut on the prefix
                d = {}
                if prefix and '_' in prefix and prefix.split('_')[0] in observables:
                    d.update(dict([prefix.split('_')]))
                m = set()
                for c, l in e:
                    d[c.GetName()] = l
                    m.add((c.GetName(), l))
                m &= par_keys
                bins = self.__bc.get_bins(level, d)
                if len(m) == 1:
                    ## A single entry is a signal state, this is either N_a or N_b.
                    ## Insert at the front so it is sure they will be either the first
                    ## or the second argument to the EffConstraint constructor.
                    k = list(m)[0]
                    eps = parameters[k]
                    args['Epsilons'].insert(0, parameters[k])
                    args['N'].insert(0, bins)
                elif len(m) == 2:
                    ## Both entries are a signal state, this must be N_ab
                    args['N'].append(bins)
            ## Calculate values of epsilon from the number of events in the different
            ## categories in each bin 
            ea, eb = args['Epsilons']
            if len(valid) == 2:
                na, nb = args['N']
                if len(ea) > 1 and len(eb) == 1:
                    ## Ensure that if there if either of the lists of epsilons
                    ## contains only a single entry, it is the first one.
                    eb, ea = ea, eb
                    nb, na = na, nb
                ## HLT1 like combination
                for e in ea: e.setConstant(True)
                if len(ea) == 1:
                    eav = len(na) * [ea[0].getVal()]
                else:
                    eav = [e.getVal() for e in ea]
                ebv = [eav[i] * nb[i] / na[i] for i in range(len(na))]
            else:
                ## HLT2 like combination
                na, nb, nab = args['N']
                if len(ea) > 1 and len(eb) == 1:
                    ## Ensure that if there if either of the lists of epsilons
                    ## contains only a single entry, it is the first one.
                    eb, ea = ea, eb
                    nb, na, nab = na, nb, nab
                eav = [nab[i] / (nb[i] + nab[i]) for i in range(len(na))]
                ebv = [nab[i] / (na[i] + nab[i]) for i in range(len(na))]
            if len(ea) == 1:
                ## If there is only one bin, use the average of the
                ## calculated values.
                ea[0].setVal(sum(eav) / float(len(eav)))
            else:
                for i, e in enumerate(eav): ea[i].setVal(e)
            for i, e in enumerate(ebv): eb[i].setVal(e)
            ## Build constraint 
            from P2VV.RooFitWrappers import EffConstraint
            constraint_name = '_'.join([e for e in (prefix, level, 'multinomial') if e])
            constraints.add(EffConstraint(Name = constraint_name, **args))
        return constraints

    def build_poisson_constraints(self, data, observables, strategy = 1):
        # get acceptance parameters and observables
        from collections import defaultdict
        binning = self._shape.base_binning()
        numBins = binning.numBins()
        coefficients = self._shape.coefficients()
        constraints = set()
        effs = defaultdict(dict)
        obs = dict([(o.GetName(), o) for o in observables.itervalues()])

        # get acceptance-shape parameters
        from P2VV.RooFitWrappers import RooObject
        for (pref, cat, state), pars in coefficients.iteritems():
            level = cat[ : 4]
            if '%s_poisson' % pref in RooObject._ws._rooobjects: continue
            effs[(pref, level)][(cat, state)] = pars
        if not effs: return constraints

        # get trigger categories
        assert len(effs) == 2 and effs.keys()[0][0] == effs.keys()[1][0] and effs.keys()[0][1] != effs.keys()[1][1]
        prefix = effs.keys()[0][0]
        levels = [k[1] for k in effs.keys()]
        trigCats = [[cat for cat in sorted(effs[(prefix, lev)].keys())] for lev in levels]
        if len(trigCats[0]) == 2 and len(trigCats[1]) == 2 and trigCats[1][0][0] == trigCats[1][1][0]:
            levels = [levels[1], levels[0]]
            trigCats = [trigCats[1], trigCats[0]]
        else:
            assert len(trigCats[0]) == 2 and len(trigCats[1]) == 2 and trigCats[0][0][0] == trigCats[0][1][0]
        assert all(cat[0] in obs for cats in trigCats for cat in cats)
        if len(effs[(prefix, levels[0])][trigCats[0][0]]) == numBins and len(effs[(prefix, levels[0])][trigCats[0][1]]) == 1:
            trigCats[0] = [trigCats[0][1], trigCats[0][0]]

        # get valid combinations of trigger categories and sums of (squared) weights for each combination
        from P2VV.Parameterizations.GeneralUtils import valid_combinations
        self.__bc = BinCounter(data, obs['time'], binning, True)
        validStates = valid_combinations([[(obs[c], s) for c, s in effs[(prefix, levels[0])].iterkeys()],
                                         [(obs[c], s) for c, s in effs[(prefix, levels[1])].iterkeys()]])
        assert len(validStates) == 6
        catStates = {}
        combCats = ['1A2A', '1A2B', '1A2AB', '1B2A', '1B2B', '1B2AB']
        for states in validStates:
            stDict = {}
            if prefix and '_' in prefix and prefix.split('_')[0] in observables:
                stDict.update(dict([prefix.split('_')]))
            for cat, state in states:
                stDict[cat.GetName()] = state
            name = combCats[3 * int(stDict[trigCats[0][1][0]] == trigCats[0][1][1])\
                   + 1 * int(stDict[trigCats[1][0][0]] == trigCats[1][0][1]) + 2 * int(stDict[trigCats[1][1][0]] == trigCats[1][1][1]) - 1]
            catStates[name] = stDict
            self.__bc.add_bins('all', stDict)
        assert all(st for st in catStates.itervalues())
        self.__bc.run()

        # create a dictionary of constraint arguments
        constrArgs = dict(Name = '%s_poisson' % prefix, Strategy = strategy, NumBins = numBins, Parameters = {},
                          SumW = dict([(n, self.__bc.get_bins('all', catStates[n])) for n in combCats]),
                          SumWSq = dict([(n, self.__bc.get_binsSq('all', catStates[n])) for n in combCats]))

        from P2VV.RooFitWrappers import RealVar
        constrArgs['Parameters']['nu'] = [RealVar(Name = '%s_yield_bin_%03d' % (prefix, bin+1), Value = 0., MinMax = (0., 1.))\
                                          for bin in range(numBins)]
        self._shape.setYields(constrArgs['Parameters']['nu'])
        constrArgs['Parameters']['eps1A'] = effs[(prefix, levels[0])][trigCats[0][0]]
        constrArgs['Parameters']['eps1B'] = effs[(prefix, levels[0])][trigCats[0][1]]
        constrArgs['Parameters']['eps2A'] = effs[(prefix, levels[1])][trigCats[1][0]]
        constrArgs['Parameters']['eps2B'] = effs[(prefix, levels[1])][trigCats[1][1]]
        assert all(len(pars) in [1, numBins] for pars in constrArgs['Parameters'].itervalues())
        for par in constrArgs['Parameters']['eps1A']:
            if hasattr(par, 'setConstant') : par.setConstant(True)

        # set estimates for the values of the acceptance parameters
        from math import sqrt
        eps1ACommon = constrArgs['Parameters']['eps1A'][0].getVal() if len(constrArgs['Parameters']['eps1A']) == 1 else -1.
        if any(len(constrArgs['Parameters'][k]) == 1 for k in ['eps1B', 'eps2A', 'eps2B']):
            # set values common to all bins
            N_1A2A  = sum(constrArgs['SumW']['1A2A'])
            N_1A2B  = sum(constrArgs['SumW']['1A2B'])
            N_1A2AB = sum(constrArgs['SumW']['1A2AB'])
            N_1B2A  = sum(constrArgs['SumW']['1B2A'])
            N_1B2B  = sum(constrArgs['SumW']['1B2B'])
            N_1B2AB = sum(constrArgs['SumW']['1B2AB'])
            if len(constrArgs['Parameters']['eps1B']) == 1 :
                if eps1ACommon > 0.:
                    eps1A = eps1ACommon
                else:
                    sumWTot = 0.
                    eps1A = 0.
                    for bIt in range(numBins):
                        sumW = sum(constrArgs['SumW'][cat][bIt] for cat in ['1A2A', '1A2B', '1A2AB', '1B2A', '1B2B', '1B2AB'])
                        sumWTot += sumW
                        eps1A += sumW * constrArgs['Parameters']['eps1A'][bIt].getVal()
                    eps1A /= sumWTot

                eps1B = (N_1B2A + N_1B2AB) / (N_1A2A + N_1A2AB) * eps1A
                constrArgs['Parameters']['eps1B'][0].setVal(eps1B)
            if len(constrArgs['Parameters']['eps2A']) == 1:
                wSq_2B = sum(constrArgs['SumWSq']['1A2B']) + sum(constrArgs['SumWSq']['1B2B'])
                if abs(wSq_2B) > 1.e-5:
                    eps2A = (N_1A2AB + N_1B2AB) / (N_1A2B + N_1B2B + N_1A2AB + N_1B2AB)
                    constrArgs['Parameters']['eps2A'][0].setVal(eps2A)
                else:
                    assert abs(N_1A2B) <= 1.e-5 and abs(N_1B2B) <= 1.e-5
                    if constrArgs['Parameters']['eps2A'][0].getMax() < 1.: constrArgs['Parameters']['eps2A'][0].setMax(1.)
                    constrArgs['Parameters']['eps2A'][0].setVal(1.)
                    constrArgs['Parameters']['eps2A'][0].setConstant()
            if len(constrArgs['Parameters']['eps2B']) == 1:
                wSq_2A = sum(constrArgs['SumWSq']['1A2A']) + sum(constrArgs['SumWSq']['1B2A'])
                if abs(wSq_2A) > 1.e-5:
                    eps2B = (N_1A2AB + N_1B2AB) / (N_1A2A + N_1B2A + N_1A2AB + N_1B2AB)
                    constrArgs['Parameters']['eps2B'][0].setVal(eps2B)
                else:
                    assert abs(N_1A2A) <= 1.e-5 and abs(N_1B2A) <= 1.e-5
                    if constrArgs['Parameters']['eps2B'][0].getMax() < 1.: constrArgs['Parameters']['eps2B'][0].setMax(1.)
                    constrArgs['Parameters']['eps2B'][0].setVal(1.)
                    constrArgs['Parameters']['eps2B'][0].setConstant()

        for bIt in range(numBins):
            # set values for one bin
            N_1A2A  = constrArgs['SumW']['1A2A'][bIt]
            N_1A2B  = constrArgs['SumW']['1A2B'][bIt]
            N_1A2AB = constrArgs['SumW']['1A2AB'][bIt]
            N_1B2A  = constrArgs['SumW']['1B2A'][bIt]
            N_1B2B  = constrArgs['SumW']['1B2B'][bIt]
            N_1B2AB = constrArgs['SumW']['1B2AB'][bIt]
            eps1A = eps1ACommon if eps1ACommon > 0. else constrArgs['Parameters']['eps1A'][bIt].getVal()
            nu = (N_1A2A + N_1A2AB) * (N_1A2B + N_1A2AB) / N_1A2AB / eps1A
            constrArgs['Parameters']['nu'][bIt].setRange((max(0., nu - 10. * sqrt(nu)), nu + 10. * sqrt(nu)))
            constrArgs['Parameters']['nu'][bIt].setVal(nu)
            if len(constrArgs['Parameters']['eps1B']) == numBins:
                eps1B = (N_1B2A + N_1B2AB) / (N_1A2A + N_1A2AB) * eps1A
                constrArgs['Parameters']['eps1B'][bIt].setVal(eps1B)
            if len(constrArgs['Parameters']['eps2A']) == numBins:
                wSq_2B = constrArgs['SumWSq']['1A2B'][bIt] + constrArgs['SumWSq']['1B2B'][bIt]
                if abs(wSq_2B) > 1.e-5:
                    eps2A = (N_1A2AB + N_1B2AB) / (N_1A2B + N_1B2B + N_1A2AB + N_1B2AB)
                    constrArgs['Parameters']['eps2A'][bIt].setVal(eps2A)
                else:
                    assert abs(N_1A2B) <= 1.e-5 and abs(N_1B2B) <= 1.e-5
                    if constrArgs['Parameters']['eps2A'][bIt].getMax() < 1.: constrArgs['Parameters']['eps2A'][bIt].setMax(1.)
                    constrArgs['Parameters']['eps2A'][bIt].setVal(1.)
                    constrArgs['Parameters']['eps2A'][bIt].setConstant()
            if len(constrArgs['Parameters']['eps2B']) == numBins:
                wSq_2A = constrArgs['SumWSq']['1A2A'][bIt] + constrArgs['SumWSq']['1B2A'][bIt]
                if abs(wSq_2A) > 1.e-5:
                    eps2B = (N_1A2AB + N_1B2AB) / (N_1A2A + N_1B2A + N_1A2AB + N_1B2AB)
                    constrArgs['Parameters']['eps2B'][bIt].setVal(eps2B)
                else:
                    assert abs(N_1A2A) <= 1.e-5 and abs(N_1B2A) <= 1.e-5
                    if constrArgs['Parameters']['eps2B'][bIt].getMax() < 1.: constrArgs['Parameters']['eps2B'][bIt].setMax(1.)
                    constrArgs['Parameters']['eps2B'][bIt].setVal(1.)
                    constrArgs['Parameters']['eps2B'][bIt].setConstant()

        # create Poisson constraint
        from P2VV.RooFitWrappers import CombEffConstraint
        constraints.add(CombEffConstraint(**constrArgs))
        return constraints

    def shapes(self):
        return [self._shape]

    def bin_counter(self):
        return self.__bc
    
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
