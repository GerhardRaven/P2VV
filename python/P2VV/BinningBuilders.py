from P2VV.RooFitWrappers import *

def build1DVerticalBinning(Name, Formula, Observable, YStep, MaxStep = 1.):
    from ROOT import Double, TH1F
    from ROOT import RooBrentRootFinder, RooRealBinding

    # Bind the expression to a function and create the root finder object
    binding = RooRealBinding(Formula._target_(), RooArgSet(Observable._target_()))
    root_finder = RooBrentRootFinder(binding)

    # Some parameter values
    obs_min = Observable.getMin()
    obs_max = Observable.getMax()

    # The bin boundaries
    from array import array
    from ROOT import vector
    boundaries = array('d')

    boundaries.append(obs_min)

    obs_orig = Observable.getVal()
    Observable.setVal(obs_max)
    eff_max = Formula.getVal()

    # Loop in the vertical direction and use the root finder to create bin boundaries.
    dt = boundaries[0]
    i = 1
    while dt < obs_max:
        low = boundaries[i - 1]
        high = low + MaxStep
        Observable.setVal(low)
        val = Formula.getVal() + YStep
        result = Double(0.)
        r = root_finder.findRoot(result, low, high, val)
        if r:
            boundaries.append(result)
        else:
            boundaries.append(high)
        dt = boundaries[i]
        i += 1
    boundaries[-1] = obs_max

    # Create a histogram from the bin boundaries and use the function to set its content
    name = Formula['Name'] + '_hist'
    _hist = TH1F(name, name, len(boundaries) - 1, boundaries)
    for i, low in enumerate(boundaries):
        high = boundaries[i + 1]
        mid = (low + high) / 2.
        Observable.setVal(mid)
        val = Formula.getVal()
        _hist.SetBinContent(i + 1, val)
        ## print '{:>5d} {: 7.5f} {: 7.5f} {: 7.5f} {: 7.5f}'.format(i, low, mid, high, val)
        if i == len(boundaries) - 2:
            break
    _hist.SetEntries(len(boundaries) - 1)

    from ROOT import RooBinning
    binning = RooBinning(len(boundaries) - 1, boundaries, Name)
    Observable.setVal(obs_orig)
    eff_func = HistFunc('_'.join((Name, Observable['Name'])), Histogram = _hist,
                                 Observables = [Observable])
    return (binning, eff_func)
