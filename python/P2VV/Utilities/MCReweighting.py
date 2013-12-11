###########################################################################################################################################
## Utilities.MCReweighting: P2VV utilities for reweighting Monte Carlo data to create desired distributions                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   VS, Vasilis Syropoulos, Nikhef, v.syropoulos@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################


# Miscelaneous helping and ploting functions

 # keep reference to canvases, legends, frames
from P2VV.Utilities.Plotting import _P2VVPlotStash

def compareDistributions( **kwargs ):
    physWeightsName = kwargs.pop('physWeightsName', None)
    obsSet          = kwargs.pop('obsSet',          None)
    itNumb          = kwargs.pop('itNumb',          None)
    nullTest        = kwargs.pop('nullTest',        None)
    
    # get datasets
    from ROOT import RooFit, RooDataSet
    data = dict( mcBefore = kwargs.pop('mcData'), Sdata = kwargs.pop('sData') )
    if kwargs.has_key('mcDataPhysRew'): data['mcAfter']    = kwargs.pop('mcDataPhysRew')
    if kwargs.has_key('MomRewData'):    data['MomRewData'] = kwargs.pop('MomRewData')

    # get observables and x ranges
    observables, Kmomenta, muMomenta, trackMomRangeX = [], [], [], {}
    for obs in obsSet:
        obsName = obs.GetName()
        if 'hel' in obsName or 'time' in obsName: observables.append(obs)
        elif obsName.startswith('K') and not obsName.startswith('KK'):
            if   'X' in obsName or 'Y' in obsName: trackMomRangeX[obsName] = (-5e3,5e3)
            elif 'Z' in obsName:                   trackMomRangeX[obsName] = (-5e2,1e5)
            else:                                  trackMomRangeX[obsName] = ( 0., 1e5)
            Kmomenta.append(obs)
        elif obsName.startswith('mu'):
            if   'X' in obsName or 'Y' in obsName: trackMomRangeX[obsName] = (-1e4,1e4)
            elif 'Z' in obsName:                   trackMomRangeX[obsName] = (-1e3,2e5)
            else:                                  trackMomRangeX[obsName] = ( 0., 2e5)
            muMomenta.append(obs)
        elif obsName == 'mdau2': KKMass = obs

    # assymetry plots are compared w.r.t. the sData
    referenceHistName = 'h_' + data['Sdata'].GetName() 

    # start drawing
    from ROOT import TCanvas, RooAbsData, TPaveText, kGreen, kMagenta
    from P2VV.Utilities.Plotting import compareDataSets, makeAssymetryPlot
    from P2VV.Load import LHCbStyle
    from math import pi

    # make canvases
    obsCanv       = TCanvas('anglesTime_%s'%itNumb,'anglesTime_%s'%itNumb)
    KaonCanv      = TCanvas('KaonMomenta_%s'%itNumb,'KaonMomenta_%s'%itNumb)
    muonCanv      = TCanvas('muonMomenta_%s'%itNumb,'muonMomenta_%s'%itNumb)    
    assymKaonCanv = TCanvas('assymKaonMomenta_%s'%itNumb,'assymKaonMomenta_%s'%itNumb)
    assymMuonCanv = TCanvas('assymmuonMomenta_%s'%itNumb,'assymmuonMomenta_%s'%itNumb)
    KKMassCanv    = TCanvas('KKMass_%s'%itNumb,'KKMass_%s'%itNumb)
    obsCanv.Divide(2,2)
    KaonCanv.Divide(4,2)
    muonCanv.Divide(4,2)
    assymKaonCanv.Divide(4,2)
    assymMuonCanv.Divide(4,2)
    
    # set some data drawing options
    colors      = dict( mcBefore = 2, mcAfter = kGreen+3, MomRewData = 4, Sdata = kMagenta+2  )
    stdDrawOpts = dict( DataError = RooAbsData.SumW2, MarkerSize = .6, XErrorSize = 0         )
    dataOpts    = dict()    
    for key in ['mcBefore','Sdata']:
        if data.has_key(key): dataOpts[key] = dict( MarkerColor = colors[key], **stdDrawOpts  )  
    for key in ['mcAfter','MomRewData']:
        if data.has_key(key): dataOpts[key] = dict( MarkerColor = colors[key], **stdDrawOpts  )  
        
    # plot angles and decay time
    print 'P2VV - INFO: compareDistributions: Plotting decay angles and time.'
    for canv, obs, logY, rangeX in zip( [obsCanv.cd(i+1) for i in range(len(observables))], 
                                        observables,
                                        3 * [False] + [True],
                                        [ (-1.,1.), (-1.,1.), (-pi,pi), (0.,14.) ] if not nullTest else \
                                            [ (-pi,pi), (-1.,1.), (-1.,1), (0.,14.) ]
                                        ): 
        compareDataSets( canv, obs, data = data, dataOpts = dataOpts, logy = logY,
                                    frameOpts = dict( Bins = 30, Range=rangeX ),
                         )

    # plot Kaon and muon momenta
    print 'P2VV - INFO: compareDistributions: Plotting track momenta.'
    for canv, assymCanv, obs in zip( 
        [ KaonCanv.cd(k+1) for k in range(len(Kmomenta)) ]      + [ muonCanv.cd(m+1) for m in range(len(muMomenta)) ],
        [ assymKaonCanv.cd(k+1) for k in range(len(Kmomenta)) ] + [ assymMuonCanv.cd(m+1) for m in range(len(muMomenta)) ],
        Kmomenta + muMomenta,
        ): 
        momFrame = compareDataSets( canv, obs, data = data, dataOpts = dataOpts,
                                               frameOpts = dict( Bins = 30, Range=trackMomRangeX[obs.GetName()] )
                                    )
        # make assymetry plots 
        makeAssymetryPlot( assymCanv, momFrame, referenceHistName ) 

    # plot KKMass 
    compareDataSets( KKMassCanv, KKMass, data = data, dataOpts = dataOpts, frameOpts = dict( Bins = 30 ))

    # make a legend and draw it
    legend, assym_legend = TPaveText( .47, .66, .77, .9, 'NDC' ), TPaveText( .269, .247, .569, .489, 'NDC' )
    legend.SetFillColor(0)
    assym_legend.SetFillColor(0)
    entriesNames = dict(mcBefore='SourceBeforePhysRew', mcAfter='SourceAfterPhysRew', MomRewData='SoourceAfterMomRew', Sdata='Target' )
    for key in ['mcBefore', 'mcAfter', 'MomRewData', 'Sdata']:
        if data.has_key(key): 
            legend.AddText('#color[%s]{%s}'%(colors[key],entriesNames[key]))
        if data.has_key(key) and key!='Sdata': 
            assym_legend.AddText('#color[%s]{%s}'%(colors[key],entriesNames[key]))
    for canv, pad in zip([obsCanv, KaonCanv, muonCanv], [4,8,8]):
        canv.cd(pad)
        legend.Draw()
    for canv, pad in zip([assymKaonCanv, assymMuonCanv], [8,8]):
        canv.cd(pad)
        assym_legend.Draw()
    #KKMassCanv.cd()
    #legend.Draw()
    _P2VVPlotStash.append(legend)
    _P2VVPlotStash.append(assym_legend)
    
    # print canvases in file
    for canv in [obsCanv,KaonCanv,muonCanv,assymKaonCanv,assymMuonCanv,KKMassCanv]: canv.Print(canv.GetName()+'.pdf')
    del mcAfterPhysRew
    return [obsCanv,KaonCanv,muonCanv,assymKaonCanv,assymMuonCanv,KKMassCanv]

# clean P2VVPlotStash to save memory
def cleanP2VVPlotStash():
    from P2VV.Utilities.Plotting import _P2VVPlotStash
    while len(_P2VVPlotStash) > 0: 
        for plot in _P2VVPlotStash: _P2VVPlotStash.remove(plot)
    print 'P2VV - INFO: cleanP2VVPlotStash: Emptied P2VVplotStash.'

# combine datasets
def _combineDataSetParts( files, name, weightName='' ):
    print 'P2VV - INFO: combineDataSetParts: Combining the following datasets with common name, %s'%name
    for f in files: print f

    # import args into current workspace
    from P2VV.RooFitWrappers import RooObject
    ws = RooObject().ws()
        
    import ROOT 
    from ROOT import TFile
    dataSets = [ TFile.Open(f,'READ').Get(name) for f in files ]
    for d in dataSets:
        ROOT.SetOwnership(d, True)
    data = dataSets.pop()
    for arg in data.get():
        if not ws[arg.GetName()]: ws.put(arg) 
    for i in range(len(dataSets)):
        d = dataSets.pop()
        data.append(d)
        del d

    print 'P2VV - INFO: combineDataSetParts: Read combined dataset with %d entries.'% data.numEntries()
    
# easily create an observable using information from the PdfConfig class
def _createGetObservable(name):
    from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis
    PdfConfig = Bs2Jpsiphi_2011Analysis()
    obsDict = PdfConfig['obsDict']
    
    if name == 'helcosthetaK': name = 'cpsi'
    if name == 'helcosthetaL': name = 'ctheta'
    if name == 'helphi'      : name = 'phi'
    if name == 'iTag'        : obsDict['iTag']     = ( 'iTag',     'Initial state flavour tag', { 'Untagged' : 0 } ) 
    if name == 'truetime'    : obsDict['truetime'] = ( 'truetime', 'true time', 'ps',            1.5,  0.,  30.   )
    if name == 'KKMassCat'   : obsDict['KKMassCat'] = ( 'KKMassCat', 'KK mass category', dict( [ ('bin%d' % i, i) for i in range(6) ] )                                 )
    
    if type( obsDict[name][2] ) != dict :
        from P2VV.RooFitWrappers import RealVar
        obs = RealVar( obsDict[name][0], Title = obsDict[name][1]
                       , Unit = obsDict[name][2], Value = obsDict[name][3]
                       , MinMax = ( obsDict[name][4], obsDict[name][5] )
                       )
    else :
        from P2VV.RooFitWrappers import Category
        obs = Category( obsDict[name][0], Title = obsDict[name][1], States = obsDict[name][2] )
    obs.setObservable(True)
    del PdfConfig
    return obs 

# Reweighting tools
class UniFunc(object):
    # Diego's class: A function that transform a variable into a flat distribution  
    def __init__(self,var,nbinsmax = None):
	""" Construct a Uniform Function
	A Uniforn function is a function which transformed values 
	have a flat distribution between [0.,1.]
        @var : a list with the initial values
        @nbinsmax: the maximum number of points in the numerical transformation
	"""      
        NMIN = 10
        if nbinsmax==None: nbinsmax= 500
  
	xlist = map(float,var)
        xlist.sort()
	n = len(xlist)
	nbins0 = int((1.*n)/NMIN)
	self.size = min(nbins0,nbinsmax)
	self.xaxis = self.size*[0.]
	self.yaxis = self.size*[0.]
        
        supmin = .5/len(var)
        cte = (1.-2*supmin)/(self.size-1)
	for i in range(self.size):
	    self.yaxis[i]=supmin + i*cte  # Pq no usas append ?
	    jevt = int(i*n*cte)  # Pq no i*(n-1)*cte y te ahorras el min() ?
	    jevt = min(jevt,n-1)
	    self.xaxis[i]=xlist[jevt]

    def value(self,x,xlist,ylist):
        """ returns the y value associated to x between the point in the xlist and y list
        """
	def bin(x,xlist):
	    """ returns the bin index in which boundaries the value of x lies in the xlist
	    """
	    x = float(x)
	    if (x<=xlist[0]): return 0,0
	    if (x>=xlist[-1]): return self.size-1,self.size-1 
	    for i in range(self.size):
		if x < xlist[i]:
		    return max(0,i-1),min(self.size-1,i)
	#print x
	x = float(x)
	#print x
	ww = bin(x,xlist)
	#print ww
	if not "__len__" in dir(ww):
		print "Crazy, " ,x, xlist[0], xlist[-1]

	i,j = ww
	x0 = xlist[i]
	y0 = ylist[i] 
	dx = xlist[j]-x0
	dy = ylist[j]-y0
	dydx = 0.
	if (i != j): dydx = dy/dx # ???????
	y = y0+dydx*(x-x0)
	return y

    def __call__(self,x):
	""" returns the transformed value of x 
	"""
	return self.value(x,self.xaxis,self.yaxis)

    def inverse(self,y):
	""" returns the inverse transformation value. 
	From a unifrom distribution to the original one
	"""
	return self.value(y,self.yaxis,self.xaxis)

def createKaonMomentaBinning(nbins, startingPoint=0, turningPoint=5e4, endpoint=35e4,typeSpec='f'):
    lowBinBounds = [] # 98% of bins are below the turningpoint
    nbins_1, nbins_2 = int(round(nbins * .98)), int(round(nbins * (1-.98)))
    binWidth_1   = ( turningPoint - startingPoint ) / nbins_1
    binWidth_2   = ( endpoint - turningPoint  )     / nbins_2
    assert nbins_1 + nbins_2 == nbins
    # create low edges for histograms bins
    for lowbin in xrange(nbins_1 + 1): lowBinBounds += [lowbin * binWidth_1] 
    for lowbin in xrange(nbins_2):     lowBinBounds += [turningPoint + (lowbin+1) * binWidth_2]
    lowBinBounds.pop(len(lowBinBounds)-1) 
    lowBinBounds.append(endpoint+1)
    assert len(lowBinBounds) == nbins + 1
    from array import array
    return array(typeSpec,lowBinBounds)

# function that reweighits the KK distributions with a 2D histrogram
def TwoDimentionalVerticalReweighting(source, target, nbins, var, weightsName, **kwargs):
    print 'P2VV - INFO: Initialised vertical reweigthing class, TwoDimentionalVerticalReweighting() for variables (%s,%s).'%(var[0],var[1])
    sourWnam = kwargs.pop('SourceWeightName', '')
    targWnam = kwargs.pop('TargetWeightName', '')
    iterIdx  = kwargs.pop('iterationNumber' ,  0)
    
    from ROOT import TH2F, TH1F, TCanvas, gROOT
  
    _valX = lambda ev: ev.find(var[0]).getVal() # value getter of the first variable
    _valY = lambda ev: ev.find(var[1]).getVal() # value getter of the first variable
    
    # get axis ranges
    sourceVar0, sourceVar1, targetVar0, targetVar1 = [],[],[],[]
    for event in source: 
        sourceVar0 += [_valX(event)] 
        sourceVar1 += [_valY(event)] 
    for event in target: 
        targetVar0 += [_valX(event)] 
        targetVar1 += [_valY(event)] 
    xMin, yMin = min(min(sourceVar0),min(targetVar0)), min(min(sourceVar1),min(targetVar1))
    xMax, yMax = max(max(sourceVar0),max(targetVar0)), max(max(sourceVar1),max(targetVar1))
    for l in [sourceVar0, sourceVar1, targetVar0, targetVar1 ]: del l
    
    # import binning
    from P2VV.Utilities.MCReweighting import createKaonMomentaBinning
    binning = createKaonMomentaBinning(nbins, endpoint=max(yMax,xMax),typeSpec='f')
    
    # create 2D histrograms (Kplus_P vs Kminus_P)
    sourceHist  = TH2F('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins, binning, nbins, binning )
    targetHist  = TH2F('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins, binning, nbins, binning )
    
    # weight getter of source distribution
    if sourWnam: source_weight = lambda ev: ev.find(sourWnam).getVal() 
    else:        source_weight = lambda ev: source.weight()
    # fill 2D rewweighting histograms
    for evnt in source: sourceHist.Fill( _valX(evnt), _valY(evnt), source_weight(evnt) )
    for evnt in target: targetHist.Fill( _valX(evnt), _valY(evnt), target.weight() )
       
    # rescale
    assert source.sumEntries()==source.numEntries(), 'P2VV - ERROR: TwoDimentionalVerticalReweighting: Rescale source weights to have sumWeighs = number of entries.'
    if source.numEntries() > target.numEntries(): sourceHist.Scale( target.sumEntries() / source.sumEntries() )
    else: targetHist.Scale( source.sumEntries() / target.sumEntries() )
    
    # calculate weights
    weights = []
    count = 0 # count how many events have a problematic weight
    for event in source:
        bin = sourceHist.FindFixBin( _valX(event),  _valY(event) ) # get the bin with given (Kplus_P,Kminus_P)
        if targetHist.GetBinContent(bin)==0 or sourceHist.GetBinContent(bin)==0:# do not weight the event if not possible with current binning
            weights += [1] 
            count += 1
        else: 
            weights += [targetHist.GetBinContent(bin) / sourceHist.GetBinContent(bin)] # calculate weight
    if count>0: print 'P2VV - INFO: TwoDimentionalVerticalReweighting: %s out of %s events are not weighted.'%(count,source.numEntries())
    
    # fill weights to histogram
    weightsHist = TH1F('Weights', 'Weights',  2*nbins, .9*min(weights), 1.1*min(weights))
    for w in weights: weightsHist.Fill(w)

    # combine weights in case source is already weighted
    if not source.isWeighted() and sourWnam:
        print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Source distribution is already weighted, combining weights.'
        from ROOT import RooArgSet
##>>>>>>> 0a01510b8be7c2f8ef7fe4644810931421e4a46a

        # put all source weights into a list
        sourceWeightList = []
        for event in source: sourceWeightList += [ source_weight(event) ]
        # combine the weights
        combinedWeights  = []
        for sourceWeight, weight in zip(sourceWeightList,weights): combinedWeights += [ sourceWeight * weight ]
        
        # remove initial weights source dataset before writting the combined ones 
        sourceColumns = RooArgSet(source.get())
        sourceColumns.remove( source.get().find(sourWnam) )
        source = source.reduce(sourceColumns)
        weightsName = sourWnam + '_' + weightsName  
    elif source.isWeighted(): print 'P2VV - ERROR: TwoDimentionalVerticalReweighting: Cannot write weights if RooDataSet is already set as weighted.'

    # scale weights to preserve number of events 
    print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Scaling sum of weights to the number of entries.'
    n_events = source.numEntries()
    sumW = sum(combinedWeights)
    scaleWeights = lambda(weight): weight * n_events / sumW 
    combinedWeights = map(scaleWeights, combinedWeights)
    assert len(combinedWeights)==source.numEntries(), 'P2VV - ERROR: TwoDimentionalVerticalReweighting: weights list and source dataset do not have the same length'
 
    # write weights to a new dataset
    source = writeWeights(source, combinedWeights, weightsName, writeDatasetName=source.GetName() + '_' +weightsName)
      
    # plot and print the 2d histograms
    canvSourc, canvTarg, canvWeights = [ TCanvas(n,n) for n in ['source','target','momWeights'] ]
    for hist, canv in zip([sourceHist,targetHist], [canvSourc,canvTarg]): 
        hist.SetStats(False)
        hist.SetAxisRange(0,4.9e4,'Y')
        hist.SetAxisRange(0,4.9e4,'X')
        canv.cd()
        hist.Draw('LEGO')
        canv.Print(canv.GetName() + '_%s.pdf'%iterIdx)
    canvWeights.cd()
    weightsHist.Draw()
    canvWeights.Print(canvWeights.GetName() + '_%s.pdf'%iterIdx )
    
    #for tree in [s,t]: tree.IsA().Destructor(tree)
    del weights, combinedWeights, sourceWeightList
    return source

# function that reweighits a single source distribution to match a given target using a histogram
def OneDimentionalVerticalReweighting(source, target, nbins, var, weightsName, **kwargs):
    print 'P2VV - INFO: Initialised vertical reweigthing class, OneDimentionalVerticalReweighting() for variable %s.'%var
    sourWnam = kwargs.pop('SourceWeightName', '')
    targWnam = kwargs.pop('TargetWeightName', '') 
    iterIdx  = kwargs.pop('iterationNumber' ,  0)

    from ROOT import TH1F, TCanvas, gROOT
    
    # get axis ranges
    gROOT.cd('PyROOT:/')
    s,t = source.buildTree(), target.buildTree(WeightName=targWnam)
    xMin = min( s.GetMinimum(var), t.GetMinimum(var) )
    xMax = max( s.GetMaximum(var), t.GetMaximum(var) )
     
    # create 2D histrograms (Kplus_P vs Kminus_P)
    sourceHist  = TH1F('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins, xMin, xMax )
    targetHist  = TH1F('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins, xMin, xMax )
  
    # value and weight getters of source distribution
    _valX = lambda ev: ev.find(var).getVal()
    if sourWnam: source_weight = lambda ev: ev.find(sourWnam).getVal() 
    else:        source_weight = lambda ev: source.weight()

    # fill rewweighting histograms
    for evnt in source: sourceHist.Fill( _valX(evnt), source_weight(evnt) )
    for evnt in target: targetHist.Fill( _valX(evnt), target.weight()     )
       
    # rescale
    assert source.sumEntries()==source.numEntries(), 'P2VV - ERROR: TwoDimentionalVerticalReweighting: Rescale source weights to have sumWeighs = number of entries.'
    if source.numEntries() > target.numEntries(): sourceHist.Scale( target.sumEntries() / source.sumEntries() )
    else: targetHist.Scale( source.sumEntries() / target.sumEntries() )
    
    # calculate weights
    weights = []
    count = 0 # count how many events have a problematic weight
    for event in source:
        bin = sourceHist.FindFixBin( _valX(event) ) # get the bin with given (Kplus_P,Kminus_P)
        if targetHist.GetBinContent(bin)==0 or sourceHist.GetBinContent(bin)==0:# do not weight the event if not possible with current binning
            weights += [1] 
            count += 1
        else: 
            weights += [targetHist.GetBinContent(bin) / sourceHist.GetBinContent(bin)] # calculate weight
    if count>0: print 'P2VV - INFO: OneDimentionalVerticalReweighting: %s out of %s events are not weighted.'%(count,source.numEntries())
    
    # scale weights to preserve number of events 
    print 'P2VV - INFO: OneDimentionalVerticalReweighting: Scaling sum of weights to the number of entries'
    n_events = source.numEntries()
    sumW = sum(weights)
    scaleWeights = lambda(weight): weight * n_events / sumW 
    weights = map(scaleWeights, weights)
    assert len(weights)==source.numEntries(), 'P2VV - INFO: TwoDimentionalVerticalReweighting: weights list and source dataset do not have the same length.'

    # fill the weights to a histogram
    weightsHist = TH1F('Weights', 'Weights',  2*nbins, .9*min(weights), 1.1*min(weights))
    for w in weights: weightsHist.Fill(w)

    # combine weights in case source is already weighted
    if not source.isWeighted() and sourWnam:
        print 'P2VV - INFO: OneDimentionalVerticalReweighting: Source distribution is already weighted, combining weights.'
        from ROOT import RooArgSet
        # put all source weights into a list
        sourceWeightList = []
        for event in source: sourceWeightList += [ source_weight(event) ]
        
        # combine the weights
        combinedWeights  = []
        for sourceWeight, weight in zip(sourceWeightList,weights): combinedWeights += [ sourceWeight * weight ]
        
        # remove initial weights source dataset before writting the combined ones 
        sourceColumns = RooArgSet(source.get())
        sourceColumns.remove( source.get().find(sourWnam) )
        source = source.reduce(sourceColumns)
        weightsName = sourWnam + '_' + weightsName  
    elif source.isWeighted(): print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Cannot write weights if RooDataSet is already set as weighted.'
    
    # # # test
    # test_s = TH1F('test_s','test_s',3*nbins, xMin,xMax)
    # test_t = TH1F('test_t','test_t',3*nbins, xMin,xMax)
    # sourceEvtList, targetEvtList = [],[]
    # for ev in source: sourceEvtList+=[ _valX(ev) ]
    # for ev in target: targetEvtList+=[ _valX(ev) ]
    # for ev, w in zip(sourceEvtList,combinedWeights): test_s.Fill(ev,w)
    # for evnt in target: test_t.Fill( _valX(evnt), target.weight()     )


    # write weights to a new dataset
    source = writeWeights(source, combinedWeights, weightsName, writeDatasetName=source.GetName() + '_' +weightsName)

   # plot and print the histograms
    # testCanv = TCanvas('test','test')
    # test_t.Draw()
    # test_s.Scale(target.sumEntries() / source.sumEntries())
    # test_s.Draw('same err')
    
    # from P2VV.Utilities.Plotting import _P2VVPlotStash
    # _P2VVPlotStash += [test_s,test_t]
    # return testCanv
    
    # canvSourc, canvTarg, canvWeights = [ TCanvas(n,n) for n in ['source','target','momWeights'] ]
    # for hist, canv in zip([sourceHist,targetHist], [canvSourc,canvTarg]): 
    #     hist.SetStats(False)
    #     #hist.SetAxisRange(0,4.9e4,'Y')
    #     #hist.SetAxisRange(0,4.9e4,'X')
    #     canv.cd()
    #     hist.Draw()
    #     canv.Print(canv.GetName() + '_%s.pdf'%iterIdx)
    # canvWeights.cd()
    # weightsHist.Draw()
    # canvWeights.Print(canvWeights.GetName() + '_%s.pdf'%iterIdx )
    
    # for tree in [s,t]: tree.IsA().Destructor(tree)
    # del weights, combinedWeights, sourceWeightList
    return source

# write weights to a RooDataSet
def writeWeights(dataset, weights, weightsName, writeDatasetName=''):
    print 'P2VV - INFO: writeWeights: Creating dataset with name %s and weight name %s:'%(writeDatasetName,weightsName)
    from ROOT import RooArgSet, RooRealVar, RooDataSet
    from P2VV.RooFitWrappers import RealVar

    weightsVar = RooRealVar( weightsName, weightsName, 1, .9*min(weights), 1.1*max(weights) )
    weightsArgSet  = RooArgSet( weightsVar )
    weightsDataSet = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
    for weight in weights:
        weightsVar.setVal( weight )
        weightsDataSet.add( weightsArgSet )
    dataset.merge( weightsDataSet ) 
    _dataset = RooDataSet(writeDatasetName, writeDatasetName, dataset.get(), 
                          Import=dataset,
                          WeightVar = (weightsName, True)
                          )
    for d in [dataset,weightsDataSet]: d.IsA().Destructor(d)
    del weights
    return _dataset

# Class for multipling a pdf with an angular acceptance and performs an sFit 
class BuildBs2JpsiKKFit():
    def __init__(self,**kwargs):
        print 'P2VV - INFO: Initialised physics reweighting class: BuildBs2JpsiKKFit().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
        pdfConfig = PdfConfig()

        self._dataSetPath = kwargs.pop('dataSetPath', None)
        self._dataSetName = kwargs.pop('dataSetName', None)
        self._weightsName = kwargs.pop('weightsName', 'JpsiKK_sigSWeight')
    
        self._doUntaggedFit = kwargs.pop('doUntaggedFit', '')        
        self._doNullTest    = kwargs.pop('doNullTest',    '')
   
        dataPath    = self._dataSetPath
        dataSetName = self._dataSetName
        dataSetFile =  self._dataSetPath
        # dataPath + 'P2VVDataSets20112012Reco14_I2DiegoMass_6KKMassBins_2TagCats.root'

        # fit options
        fitOpts = dict(  NumCPU    = 8
                         , Optimize  = 2
                         , Minimizer = 'Minuit2'
                         )
        pdfConfig['fitOptions'] = fitOpts
        corrSFitErrCats         = [ 'runPeriod', 'KKMassCat' ]
        randomParVals           = ( ) # ( 1., 12345 )

        # PDF options
        pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
        pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
            = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
        pdfConfig['angEffMomsFiles'] = '' #dataPath + 'Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'
        pdfConfig['anglesEffType']   = ''

        # read data set from file
        from P2VV.Utilities.DataHandling import readData
        dataSet = readData( filePath = dataSetFile, dataSetName = dataSetName,  NTuple = False )
        pdfConfig['signalData'] = dataSet
        pdfConfig['readFromWS'] = True
        
        # build the PDF
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
        pdfBuild = PdfBuilder( **pdfConfig )
        pdf = pdfBuild.pdf()

        # data set with weights corrected for background dilution: for phi_s fit only!
        from P2VV.Utilities.DataHandling import correctWeights
        self._fitData = correctWeights( dataSet, corrSFitErrCats )

        # fix values of some parameters
        for CEvenOdds in pdfBuild['taggingParams']['CEvenOdds'] :
            if not pdfConfig['SSTagging'] :
                CEvenOdds.setConstant('avgCEven.*')
                CEvenOdds.setConstant( 'avgCOdd.*', True )
            else :
                for CEvenOdd in CEvenOdds :
                    CEvenOdd.setConstant('avgCEven.*')
                    CEvenOdd.setConstant( 'avgCOdd.*', True )

        pdfBuild['amplitudes'].setConstant('C_SP')

        # print parameters
        print 120 * '='
        print 'Bs2JpsiKKFit: fit data:'
        self._fitData.Print()
        #print 'Bs2JpsiKKFit: observables in PDF:'
        #pdfObs.Print('v')
        #print 'Bs2JpsiKKFit: parameters in PDF:'
        #pdfPars.Print('v')
        print 'Bs2JpsiKKFit: constraints in PDF:'
        for constr in pdf.ExternalConstraints() : constr.Print()

        self._FitResults  = {} # collect all the fit results
        self._Moments     = {} # collect all moments
       
        self._pdfConfig = pdfConfig
        self._pdfBuild  = pdfBuild
        self._pdf = pdf

    def doFit( self, itNum=0, angAccFile=None ):
        pref = self._pdfConfig['parNamePrefix'] + '_'
        
        # multiply by angular acceptance
        if angAccFile: self._pdf = self._multiplyPdfWithAcc( angAccFile, iterNumb = itNum )

        print 120 * '='
        # fit data
        print 'Bs2JpsiKKFit: fitting %d events (%s)' % (  self._fitData.numEntries(), 'weighted' if  self._fitData.isWeighted() else 'not weighted' )
        fitResult = self._pdf.fitTo( self._fitData, SumW2Error = False, Save = True, Hesse = False, Offset = False, ** self._pdfConfig['fitOptions'] )
        
        # print parameter values
        from P2VV.Imports import parNames, parValues
        print 'Bs2JpsiKK2011Fit: parameters:'
        fitResult.SetName('sFit_%s_Iteration'%itNum)
        fitResult.PrintSpecial( text = False, LaTeX = False, normal = True, ParNames = parNames, ParValues = parValues )
        #fitResult.covarianceMatrix().Print()
        #fitResult.correlationMatrix().Print()
        from ROOT import TFile
        resultFile = TFile.Open('fitResult_iterativeProcedure_%s.root'%itNum,'recreate')
        resultFile.cd()
        fitResult.Write()
        resultFile.Close()
        self._FitResults['iter_%s'%itNum] = fitResult 
    
        print 120 * '=' + '\n'        

    def updateDataParameters(self, oldPars, itNum=0):
        fitResult = self._FitResults['iter_%s'%itNum]
        cloneOldPars = oldPars.copy()
        parkeys = oldPars.keys()
        parkeys.remove('A0Phase')
        for par in parkeys:
            try: oldPars[par] = fitResult.floatParsFinal().find( self._pdfConfig['parNamePrefix'] + '_' + par).getVal()
            except AttributeError: print 'P2VV - WARNING: updateDataParameters: Parameter %s not found in fit result %s using, parameter value will not change.'\
                    %( self._pdfConfig['parNamePrefix'] + '_' + par, fitResult.GetName() )
        print 'P2VV - INFO: updateDataParameters: Data physics parameters updated from sFit.'
        for k in parkeys: print '%20s  %.4f --> %.4f'%(k,cloneOldPars[k],oldPars[k])
        
    def _multiplyPdfWithAcc( self, effFile, iterNumb=None ):
        # read moments file and multiply pure pdf with angular acceptance
        print 'P2VV - INFO:multiplyPdfWithAcc(): multiplying PDF "%s" with angular efficiency moments from file "%s"'\
          % ( self._pdf.GetName(), effFile )
        from P2VV.Utilities.DataMoments import angularMomentIndices, RealMomentsBuilder
        moments = RealMomentsBuilder()
        moments.appendPYList( self._pdfBuild['angleFuncs'].angles, angularMomentIndices( 'weights', self._pdfBuild['angleFuncs'] ) )
        moments.read( effFile )    
        self._Moments['%s_iteration'%iterNumb] = moments # collect all moments  
        if iterNumb: return moments.multiplyPDFWithEff( self._pdf, CoefName = self._pdfConfig['parNamePrefix'] + '_' + 'effC%d' % iterNumb )
        else:        return moments.multiplyPDFWithEff( self._pdf, CoefName = self._pdfConfig['parNamePrefix'] + '_' + 'effnom'            )
        
         
    def getDataSet(self):           return self._pdfConfig['signalData']
    def getPdf(self):               return self._pdf
    def getPdfBuilderObject(self) : return self._pdfBuild
    def getFitResults(self):        return self._FitResults
    def getMomentsDict(self):       return self._Moments['%s_iteration'%itNum]
    def getObservables(self,which=None):
        if   which=='angles': return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName().startswith('hel') ]
        elif which=='time':   return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName()==('time')         ]
        else:                 return             self._pdfBuild['obsSetP2VV']

# Vertical reweighting class to match physics of weighted distribution using a pdf
class MatchPhysics( ):
    def __init__( self, nTupleFile, nTupleName, **kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'
        
        # set global object name prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        setParNamePrefix( 'mc' )

        # transversity amplitudes
        A0Mag2Val    = 0.722**2 / (0.722**2 + 0.480**2 + 0.499**2)
        AperpMag2Val = 0.499**2 / (0.722**2 + 0.480**2 + 0.499**2) 
        AparMag2Val  = 0.480**2 / (0.722**2 + 0.480**2 + 0.499**2)

        A0PhVal    = 0.
        AperpPhVal = 3.07
        AparPhVal  = 3.30

        # build pdf with KK mass bins and category states.
        KKMassStates = dict( [ ('bin%d' % i, i) for i in range(6) ] )
        SWaveAmps = dict(  mc_f_S        = dict()
                         , mc_ASOddPhase = dict()
                           , mc_C_SP       = dict()
                           )
        for k in SWaveAmps.keys():
            for bin in xrange(6): SWaveAmps[k]['bin%s'%bin] = 0

        # CP violation parameters
        phiCPVal    = +0.07 
        lambCPVal = 1.

        # B lifetime parameters
        GammaVal  = 1. / 1.503 
        dGammaVal = 1. / 1.406 - 1. / 1.614
        dMVal     = 17.8
        tResSigma = 0.045

        angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{#mu})', '#phi_{h}'  )

        # load roofit wrappers
        from P2VV.Load import RooFitOutput
        
        # get workspace
        from P2VV.RooFitWrappers import RooObject
        ws = RooObject().ws()
 
        # TODO: This is a hack. Check why you REALLY need to do that
        from math import pi
        ws['helphi'].setRange(-pi,pi)
       
        # angular functions
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )
        self._angleFuncs = angleFuncs

        # get observables
        trueTime  = _createGetObservable('truetime') # get the wrapper instead of the base object !!
        angles    = [_createGetObservable(o) for o in ['helcosthetaK','helcosthetaL','helphi'] ]
        iTag      = _createGetObservable('iTag')
        KKMass = _createGetObservable('KKMass')
        KKMassCat = RooObject._rooobject('KKMassCat')
        runPeriod = RooObject._rooobject('runPeriod')
        self._obsSet = [ trueTime, KKMass ] + angles + [ runPeriod, KKMassCat ]
        
        # add track momenta and B momenta
        self._obsSet += [ RooObject._rooobject('%s_%s'%( part, comp )) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'PX', 'PY', 'PZ', 'P' ) ]
        self._obsSet += [ RooObject._rooobject('B_%s'%comp) for comp in [ 'P', 'Pt' ] ]

        # read ntuple
        cuts = '(runNumber<2524e3) || (runNumber>2546e3)'
        if cuts: print 'P2VV - INFO: MatchPhysics: Accelerate developing precces by cutting on run time ' 

        # read into dataset
        from P2VV.Utilities.DataHandling import readData
        self._data = readData( nTupleFile, dataSetName = nTupleName, NTuple = True, observables = self._obsSet, ntupleCuts = cuts )

        # build pdf
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
        amplitudes = Amplitudes(  A0Mag2 = A0Mag2Val, AperpMag2 = AperpMag2Val
                                  , AparPhase = AparPhVal - A0PhVal, AperpPhase = AperpPhVal - A0PhVal
                                  , f_S = 0., ASOddPhase = 0., C_SP = 1.
                                  )

        # B lifetime
        from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVal, dM = dMVal )
        
        # resolution model
        tResArgs = { }
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = trueTime
        timeResModel = TimeResolution( **tResArgs )

        # CP violation parameters
        from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
        lambdaCP = CPParam( lambdaCP = lambCPVal, phiCP = phiCPVal )
        
        # tagging parameters
        from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
        taggingParams = TaggingParams()

        # coefficients for time functions
        from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

        # build signal PDF
        args = dict(    time                   =  trueTime
                      , iTag                   = iTag
                      , tau                    = lifetimeParams['MeanLifetime']
                      , dGamma                 = lifetimeParams['dGamma']
                      , dm                     = lifetimeParams['dM']
                      , dilution               = taggingParams['dilution']
                      , ADilWTag               = taggingParams['ADilWTag']
                      , avgCEven               = taggingParams['avgCEven']
                      , avgCOdd                = taggingParams['avgCOdd']
                      , coshCoef               = timeBasisCoefs['cosh']
                      , sinhCoef               = timeBasisCoefs['sinh']
                      , cosCoef                = timeBasisCoefs['cos']
                      , sinCoef                = timeBasisCoefs['sin']
                      , resolutionModel        = timeResModel['model']
                      )

        from P2VV.RooFitWrappers import BTagDecay
        protoPdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

        splitCats, splitPars = [ ], [ ]
        splitCats.append( [ KKMassCat ] )
        splitPars.append( [ ws[par] for par in SWaveAmps.iterkeys() ] )
        
        # split PDF for run period and KK-mass bins
        from P2VV.RooFitWrappers import SimultaneousPdf
        self._pdf = SimultaneousPdf( protoPdf.GetName() + '_simul', MasterPdf = protoPdf
                               , SplitCategories = splitCats, SplitParameters = splitPars )
        
        # set values of split parameters
        for par, vals in SWaveAmps.iteritems() :
            for cat, val in vals.iteritems() : ws[ par + '_' + cat ].setVal(val)

        # monte carlo gen paramters
        MCProd = kwargs.pop('MonteCarloProduction')
        if MCProd == 'Sim08':
            from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as mcPars
        self._mcPhysParsSet = mcPars
        
        # create pdf physics parameters dictionary
        self._pdfPhysPars = {}
        for par in self._pdf.Parameters():
            if par.GetName().partition('_')[2].startswith('A') and \
                    ('Mag2' in par.GetName().partition('_')[2] or 'Phase' in par.GetName().partition('_')[2] ): 
                self._pdfPhysPars[par.GetName()] = par
            elif 'f_S' in par.GetName(): 
                self._pdfPhysPars[par.GetName()] = par
        for par in self._pdf.Parameters():
            for physParName in ['dM', 'dGamma', 'Gamma', 'phiCP', 'lambdaCP']:
                 if physParName in par.GetName(): 
                     self._pdfPhysPars[par.GetName()] = par
    
    def setMonteCarloParameters(self, pars=None):
        print 'P2VV - INFO: setMonteCarloParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        for key in self._pdfPhysPars.keys():
            self._pdfPhysPars[key].setVal( self._mcPhysParsSet[ key.partition('_')[2]] ) 
        for k in self._pdfPhysPars.keys(): 
            print '%20s %.4f'%(self._pdfPhysPars[k].GetName(), self._pdfPhysPars[k].getVal())
                    
    def setDataFitParameters(self, dataPars, KKmassCat=None):
        print 'P2VV - INFO: setDataFitParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        for key in self._pdfPhysPars.keys():
            self._pdfPhysPars[key].setVal( dataPars[ key.partition('_')[2]] ) 
        for k in self._pdfPhysPars.keys(): 
            print '%20s %.4f'%(self._pdfPhysPars[k].GetName(), self._pdfPhysPars[k].getVal())
        else: print ' ' # TODO: LoopOver the KKbins and set f_S_i and ASOddPhase_i

    def calculateWeights(self, iterNumb, dataParameters):
        self._iterNumb = iterNumb

        from ROOT import RooArgSet
        normVars =  RooArgSet( self._obsSet )
        
        # Reweights MC verticaly to match the Physics of data.
        nominators, denominators = [], []
        
        def calculatePhysicsWeights(nom=[],den=[]):
            self._physWeights = []
            count = 0
            for idx in xrange(len(nom)): # if pdf.getVal() returns 0 for a given event make it unweighted
                if nom[idx] == 0 or den[idx] == 0: 
                    nom[idx], den[idx] = 1., 1.
                    count += 1 # count how many evetns have a problematic weight assignment 
                self._physWeights += [ nom[idx]/den[idx] ]
            if count>0:print 'P2VV - WARNING: calculateWeights: For %s events out of %s pdf value is zero.'%(count,len(nom))
            
            # scale weights to preserve sample size.
            n_events = self._data.numEntries()
            sumW = sum(self._physWeights)
            scaleWeights = lambda(weight): weight * n_events / sumW 
            self._physWeights = map(scaleWeights, self._physWeights)

        # loop over events and calculate physics weights
        self._pdf.attachDataSet( self._data ) # make the pdf directly dependant on data
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'    
        self.setMonteCarloParameters()
        for nev in xrange(self._data.numEntries()):
            self._data.get(nev)
            denominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating nominators for phyisics matching weights'
        self.setDataFitParameters(dataParameters)
        for nev in xrange(self._data.numEntries()):
            self._data.get(nev)
            nominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating phyisics matching weights'
        calculatePhysicsWeights(nom=nominators, den=denominators)
        
    def writeWeights(self, weightsName='weightPhys'):
        from ROOT import RooArgSet, RooRealVar, RooDataSet
        from P2VV.RooFitWrappers import RealVar

        self._weightsName = weightsName
        weightsVar        = RooRealVar( self._weightsName, self._weightsName, 1, .9*min(self._physWeights), 1.1*max(self._physWeights) )
        weightsArgSet     = RooArgSet( weightsVar )
        weightsDataSet    = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        
        for weight in self._physWeights:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        
        del self._physWeights

        self._data.merge( weightsDataSet )
        del weightsDataSet

        self._data.SetName('MC_AfterPhysRew_%s_iteration'%self._iterNumb )
        print 'P2VV - INFO: Phyisics matching weights added to dataset: MC_AfterPhysRew_%s_iteration'%self._iterNumb
        self._weightedData = RooDataSet( self._data.GetName(),self._data.GetTitle(),
                                         self._data.get(), 
                                         Import    = self._data,
                                         WeightVar = (self._weightsName, True) )
               
    def getPdf(self):               return self._pdf
    def getAngleFunctions(self):    return angleFuncs
    def getMcTime(self):            return  [ o for o in self._pdf.Observables() if 'time' in o.GetName() ]
    def getParNamePrefix(self):     return 'mc'
    def getDataSet(self, weighted=False):          
        if weighted: return self._weightedData
        else:        return self._data
    def getProjDataset(self): return getDataset.reduce(getPdf.ConditionalObservables() + self.pdf.indexCat())
    def plotWeights(self, Range=() ):
        from ROOT import TCanvas
        from P2VV.Utilities.Plotting import _P2VVPlotStash
        c = TCanvas('PhysWeights_%s'%self._iterNumb, 'PhysWeights_%s'%self._iterNumb)
        frame = self.getDataSet().get().find(self._weightsName).frame(Bins=150,Range=(.6,1.6))
        self.getDataSet().plotOn(frame, MarkerSize=.5, XErrorSize=0.)
        frame.Draw()
        c.Print('physicsWeights_%s.pdf'%self._iterNumb)
        _P2VVPlotStash +=[c]
     

# Class that matches weighted distributions with horizontal reweighting and then recalculates the decay angles.
class MatchWeightedDistributions():
    def __init__( self,  **kwargs ):
        print 'P2VV - INFO: Initialised kinematic reweighting class: matchWeightedDistributions().'
        self._inTree          = kwargs.pop('inTree', '')  # distribution to be modified.
        self._outTree         = kwargs.pop('outTree')     # distribution to be matched with.
        self._inWeightName    = kwargs.pop('inWeightName')
        self._outWeightName   = kwargs.pop('outWeightName')
        self._mcObsSet        = kwargs.pop('observables')
        self._copyVars        = kwargs.pop('nonObsVars', None)
        self._mimicedVars     = dict( inDistr={} )
        
        self._nBins = kwargs.pop('nBins',          '1000'  )
        self._vars  = kwargs.pop('reweightVars', 'Kminus_P')
        self._itNum = kwargs.pop('itNum',            0     )
       
        from ROOT import RooDataSet, gROOT
        if self._inTree and type(self._inTree) == RooDataSet:
            self._inSumW = self._inTree.sumEntries()
            gROOT.cd('PyROOT:/')
            self._inTree  = self._inTree.buildTree( WeightName=self._inWeightName)
        if type(self._outTree) == RooDataSet:
            self._outSumW = self._outTree.sumEntries()
            gROOT.cd('PyROOT:/')
            self._outTree = self._outTree.buildTree( WeightName=self._outWeightName )
                 
        # make some usefull sets.
        if self._copyVars:
            self._muonSet, self._KaonSet, self._BmomSet = [], [],[]
            for var in self._copyVars:
                if var.GetName().startswith('mu'):      self._muonSet.append(var)
                elif var.GetName().startswith('K'):     self._KaonSet.append(var)
                elif var.GetName().startswith('B'):     self._BmomSet.append(var)
                elif var.GetName().startswith('mdau2'): self._KKMass = [var]
            
    def _mimicWeights(self):
        for var in self._vars:
            print 'P2VV - INFO: _mimicWeights: Mimicing weighteed distribution of source variable:', var
            self._mimicedVars['inDistr'][var] = self._MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins, 'in')
            
            if self._mimicedVars.has_key('outDistr'):
                if not self._mimicedVars['outDistr'].has_key(var):
                    print 'P2VV - INFO: _mimicWeights: Mimicing weighteed distribution of target variable:', var
                    self._mimicedVars['outDistr'][var] = self._MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out') 
            else:
                print 'P2VV - INFO: _mimicWeights: Mimicing weighteed distribution of target variable:', var
                self._mimicedVars['outDistr'] = {}
                self._mimicedVars['outDistr'][var] = self._MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out')            

    def _MimicWeightedDistribution( self, t, var, wPref, Nbins, varSpec ):
        from ROOT import gDirectory, TH1F, TRandom, TCanvas
        gDirectory.cd()

        # Draw histogram to be mimiced
        nam = 'h_' + var + str(varSpec) + str(self._itNum)
        c = TCanvas(nam,nam)
        c.cd()
        t.Draw( '%s>>%s(%s,%s,%s)'%(var,nam,Nbins,t.GetMinimum(var),t.GetMaximum(var)), wPref )
        hist = gDirectory.Get(nam)

        # construct mimic histogram
        mimicHist = TH1F('mimic_' + hist.GetName(), 'mimic_' + hist.GetTitle(), Nbins, t.GetMinimum(var),t.GetMaximum(var)) 
        newDistribution = []  # List of new mimiced distribution
        
        # generate n random numbers within each bin range untill n = n_th bin content.   
        rndm = TRandom()
        for b in xrange( 1, hist.GetNbinsX() ):
            count = 1
            while count <= hist.GetBinContent(b):
                number    = rndm.Uniform( hist.GetBinLowEdge(b), hist.GetBinLowEdge(b) + hist.GetBinWidth(b) )                 
                mimicHist.Fill(number)
                newDistribution += [number]
                count += 1
            if hist.GetBinContent(b) != mimicHist.GetBinContent(b) and hist.GetBinContent(b) - mimicHist.GetBinContent(b) > 1 :
                print 'P2VV - WARNING: MimicWeightedDistribution: Mimicinc of bin %s is of by %s events. Original bin content is %s'\
                    %(b,abs(hist.GetBinContent(b)- mimicHist.GetBinContent(b)), hist.GetBinContent(b))
        
        del c
        return newDistribution

    def reweight( self, itNum, data, reweightingType='horizontal' ):
        # get rid of the previous tree if any
        self._inTree = data
        self._itNum = itNum
        
        # get physics weight variable
        if not self._inTree.get().find( self._inWeightName ):
            from ROOT import RooRealVar, RooNumber
            RooInf = RooNumber.infinity()
            self._physWeightsVar = RooRealVar( self._inWeightName, self._inWeightName, -RooInf, RooInf )
        else: self._physWeightsVar = self._inTree.get().find( self._inWeightName )

        # get dataset scales
        self._inSumW = self._inTree.sumEntries()
        self._dataSetsScale = self._outSumW / self._inSumW
        self._physRewScale  = self._inSumW / self._inTree.numEntries()
        assert self._physRewScale < 1.01, \
            'P2VV - WARNING: Scale between original MC and physics reweighted MC is %s, remove weights and then perform momentum reweighting.' %self._physRewScale 
        
        # convert roodataset to tree
        from ROOT import RooDataSet
        from ROOT import gROOT
        gROOT.cd('PyROOT:/')
        self._inTree = self._inTree.buildTree(WeightName=self._inWeightName)

        # mimic the weights, transform Kaon momenta and recalculate angles 
        if reweightingType == 'horizontal':
            self._mimicWeights()
            self._TransformKaonMomentaAndRecalculateAngles()
        else: print 'do verical reweighting'

    def _constructTransformation( self ):
        source = self._mimicedVars['inDistr'][self._vars[0]]
        target = self._mimicedVars['outDistr'][self._vars[0]]
        
        # Umc, Udat transform the distributions into flat ones
        self._Tsource = UniFunc(source, nbinsmax = self._nBins)
        self._Ttarget = UniFunc(target, nbinsmax = self._nBins)
    
    def _TransformKaonMomenta( self, k1, k2 ):
        # apply the inverse of the target transformation to the source  
        trans = lambda x: self._Ttarget.inverse( self._Tsource(x) )
        return ( trans(k1), trans(k2) )

    def _TransformKaonMomentaAndRecalculateAngles( self ):
        # put the recalculated angles plus time in a RooDataSet
        from ROOT import RooDataSet, RooArgSet, RooFit
        for obs in self._mcObsSet:
            if   obs.GetName()=='helcosthetaK': helcosthetaK = obs
            elif obs.GetName()=='helcosthetaL': helcosthetaL = obs
            elif obs.GetName()=='helphi':          helphi    = obs
            elif obs.GetName()=='time':             time     = obs
            elif obs.GetName()=='truetime':       truetime   = obs

        # get references to Kaon variables 
        for var in self._KaonSet: 
            if 'plus' in var.GetName():
                if   'X' in var.GetName(): KPX = var
                elif 'Y' in var.GetName(): KPY = var
                elif 'Z' in var.GetName(): KPZ = var
                elif var.GetName().endswith('P'):     KP  = var
            elif 'minus' in var.GetName():
                if   'X' in var.GetName(): KmX = var
                elif 'Y' in var.GetName(): KmY = var
                elif 'Z' in var.GetName(): KmZ = var
                elif var.GetName().endswith('P'):     Km  = var

        # important bookeeping check, make sure you do not flip B_P and B_Pt
        assert(self._BmomSet[0].GetName()=='B_P' and self._BmomSet[1].GetName()=='B_Pt')
        BP, BPt = self._BmomSet[0], self._BmomSet[1]
        
        # construct datasets 
        recalculatedVars = RooArgSet( [helcosthetaK,helcosthetaL,helphi] + self._KaonSet + self._BmomSet ) 
        copiedVars       = RooArgSet( self._KKMass + self._muonSet + [truetime,time,self._physWeightsVar] )
        recalculatedData = RooDataSet( 'MomRewMC_%s_Iter'%self._itNum, 'MomRewMC_%s_Iter'%self._itNum, recalculatedVars         )        
        copiedData       = RooDataSet( 'copiedData',                   'copiedData',                   self._inTree, copiedVars )

        # get particle masses
        from ROOT import TDatabasePDG
	MeV = 1000 # TDatabasePDG is in GeV, this is the factor needed to go to MeV.
        PDG = TDatabasePDG()
        Mmu = PDG.GetParticle('mu-').Mass()*MeV
        Mk  = PDG.GetParticle('K-').Mass()*MeV
        
        # define useful functions
        from ROOT import TVector3, TLorentzVector, HelicityAngles 
        from math import sqrt
        _VM2LV  = lambda v,m : TLorentzVector( v, sqrt( m*m + v.Mag2() ) ) # TVector3,mass to TLorentzVector 
        _E2V    = lambda entry, label : TVector3( getattr(entry,label+'_PX'),getattr(entry,label+'_PY'),getattr(entry,label+'_PZ')) # entry to TVenctor3
        _B3PMAG = lambda K1,K2,mu1,mu2: (K1+K2+mu1+mu2).Mag() # B momentum.
        _BPT    = lambda K1,K2,mu1,mu2: sqrt( (K1+K2+mu1+mu2).x()**2 + (K1+K2+mu1+mu2).y()**2 ) # B transverse momentum.
               
        # build the Kmomenta transformation
        print 'P2VV - INFO: _ReweightAndTransformAngles: Building Kaon momenta transformation.'
        self._constructTransformation()
        print 'P2VV - INFO: _ReweightAndTransformAngles: Transforming K+, K- momenta and recalculating decay angles.'
        for entry in self._inTree:
            k1_3P = _E2V(entry,'Kplus')
            k2_3P = _E2V(entry,'Kminus')
            mag1,mag2 = self._TransformKaonMomenta( k1_3P.Mag(), k2_3P.Mag() )
            k1_3P.SetMag( mag1 )
            k2_3P.SetMag( mag2 )
            
            mu1_3P = _E2V( entry , 'muplus' )
            mu2_3P = _E2V( entry , 'muminus')
           
            helangles = HelicityAngles( _VM2LV( k1_3P,  Mk  ), _VM2LV( k2_3P,  Mk  ),
                                        _VM2LV( mu1_3P, Mmu ), _VM2LV( mu2_3P, Mmu )
                                        )
           
            # fill the new dataset.
            helcosthetaK.setVal( helangles[0] )
            helcosthetaL.setVal( helangles[1] )
            helphi.setVal(       helangles[2] )
            
            KPX.setVal( k1_3P.x() )
            KPY.setVal( k1_3P.y() )
            KPZ.setVal( k1_3P.z() )
            KP.setVal( k1_3P.Mag() )
            
            KmX.setVal( k2_3P.x() )
            KmY.setVal( k2_3P.y() )
            KmZ.setVal( k2_3P.z() )
            Km.setVal( k2_3P.Mag() )

            BP.setVal( _B3PMAG(k1_3P,k2_3P,mu1_3P,mu2_3P) )
            BPt.setVal( _BPT(k1_3P,k2_3P,mu1_3P,mu2_3P)   )

            recalculatedData.addFast( recalculatedVars )
        
        recalculatedData.merge( copiedData )
        self._recalculatedData = RooDataSet(recalculatedData.GetName(), recalculatedData.GetTitle(), 
                                            RooArgSet(recalculatedVars, copiedVars), 
                                            Import = recalculatedData, 
                                            WeightVar = (self._physWeightsVar, True)
                                            )
        self._recalculatedDataNoPhysWeights = recalculatedData
  
    def getDataSet(self, tree=False, noPhysWeights=False): 
        if tree:
            from ROOT import gROOT
            gROOT.cd('PyROOT:/')
            return self._recalculatedData.buildTree(WeightName = self._inWeightName)
        else: 
            if noPhysWeights: return self._recalculatedDataNoPhysWeights
            else:             return self._recalculatedData

# physics parameter imports
parValues6KKmassBins20112012 = dict( # sFit in 2011 + 2012 data  
                                      A0Mag2            = 5.2085e-01
                                     ,A0Phase           = 0. 
                                     ,ASOddPhase_bin0   = 8.0380e-01 
                                     ,ASOddPhase_bin1   = 2.3096e+00
                                     ,ASOddPhase_bin2   = 4.5944e-01
                                     ,ASOddPhase_bin3   =-3.6016e-01
                                     ,ASOddPhase_bin4   =-6.7912e-01
                                     ,ASOddPhase_bin5   =-8.2363e-01
                                     ,AparPhase         = 3.2547e+00 
                                     ,AperpMag2         = 2.5354e-01
                                     ,AperpPhase        = 3.1874e+00
                                     ,Gamma             = 6.6062e-01
                                     ,dGamma            = 8.7053e-02
                                     ,phiCP             = 8.2805e-02
                                     ,betaTimeEff_p2011 =-8.2950e-03
                                     ,betaTimeEff_p2012 =-1.3746e-02
                                     ,dM                = 1.7765e+01 
                                     ,f_S_bin0          = 4.4442e-01 
                                     ,f_S_bin1          = 6.3233e-02
                                     ,f_S_bin2          = 8.6777e-03
                                     ,f_S_bin3          = 8.8876e-03
                                     ,f_S_bin4          = 4.4337e-02
                                     ,f_S_bin5          = 2.0997e-01
                                     ,lambdaCP          = 9.6384e-01
                                     ,wTagDelP0OS       = 1.3725e-02
                                     ,wTagDelP0SS       =-1.5955e-02
                                     ,wTagDelP1OS       = 6.9765e-02
                                     ,wTagDelP1SS       = 1.5007e-02
                                     ,wTagP0OS          = 3.9049e-01
                                     ,wTagP0SS          = 4.4015e-01
                                     ,wTagP1OS          = 1.0363e+00
                                     ,wTagP1SS          = 9.4426e-01
                                     )

parValuesNoKKBinsWideKKWindow = dict(  # sFit in 2011 data
                                       AperpMag2        = 0.24663
                                      ,AperpPhase       = 3.0232
                                      ,A0Mag2           = 0.52352
                                      ,A0Phase          = 0 # cosntrained
                                      ,AparPhase        = 3.2144
                                      ,f_S              = 0.0463
                                      ,ASOddPhase       = -0.0666
                                      ,dM               = 17.6739
                                      ,dGamma           = 0.1021
                                      ,Gamma            = 0.6728
                                      ,phiCP            = 0.0847
                                      ,lambdaCP         = 0.9275
                                       )

parValuesMcSim08_6KKmassBins = dict(  # Sim08 generating conditions  
                                      A0Mag2            = 0.722**2 / (0.722**2 + 0.480**2 + 0.499**2)
                                     ,A0Phase           = 0. 
                                     ,ASOddPhase_bin0   = 0.
                                     ,ASOddPhase_bin1   = 0.
                                     ,ASOddPhase_bin2   = 0.
                                     ,ASOddPhase_bin3   = 0.
                                     ,ASOddPhase_bin4   = 0.
                                     ,ASOddPhase_bin5   = 0.
                                     ,AparPhase         = 3.30
                                     ,AperpMag2         = 0.499**2 / (0.722**2 + 0.480**2 + 0.499**2) 
                                     ,AperpPhase        = 3.07
                                     ,Gamma             = 1. / 1.503 
                                     ,dGamma            = 1. / 1.406 - 1. / 1.614
                                     ,phiCP             = +0.07 
                                     ,dM                = 17.8
                                     ,f_S_bin0          = 0.
                                     ,f_S_bin1          = 0.
                                     ,f_S_bin2          = 0.
                                     ,f_S_bin3          = 0.
                                     ,f_S_bin4          = 0.
                                     ,f_S_bin5          = 0.
                                     ,lambdaCP          = 1.
                                      )

parValuesMc2011Gen =  dict( AperpMag2        = 0.16
                           ,AperpPhase       = -0.17
                           ,A0Mag2           =  0.6
                           ,A0Phase          = 0 # cosntrained
                           ,AparPhase        = 2.50
                           ,f_S         = 0
                           ,ASOddPhase  = 0
                           ,dM               = 17.8
                           ,dGamma           = 0.06
                           ,Gamma            = 0.679
                           ,phiCP            = -0.04
                           ,lambdaCP         = 1.
                            )
parValuesMc2012Gen =  dict( AperpMag2        = 0.2488
                           ,AperpPhase       = 3.07
                           ,A0Mag2           = 0.5209
                           ,A0Phase          = 0. # cosntrained
                           ,AparPhase        = 3.30
                           ,f_S         = 0.
                           ,ASOddPhase  = 0.
                           ,dM               = 17.8
                           ,dGamma           = 0.0917
                           ,Gamma            = 0.6653
                           ,phiCP            = 0.07
                           ,lambdaCP         = 1.
                            )
parValuesMc2012Fit =  dict(# '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting_part1.root'
                           #  AperpMag2        =  2.5628e-01 
                           # ,AperpPhase       = -2.8377e+00 
                           # ,A0Mag2           =  5.1761e-01
                           # ,A0Phase          =  0. # cosntrained
                           # ,AparPhase        =  3.2728e+00 
                           # ,f_S_bin0         =  0.
                           # ,ASOddPhase_bin0  =  0.
                           # ,dM               =  1.7630e+01 
                           # ,dGamma           =  8.0522e-02 
                           # ,Gamma            =  6.8023e-01
                           # ,phiCP            = -1.2965e-01
                           # ,lambdaCP         =  1.0684e+00
                           
                           # For null test
                           #'/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting.root'
                           # using 2012 angular acceptance weights
                            AperpMag2        = 2.4841e-01
                           ,AperpPhase       = 3.4126e+00
                           ,A0Mag2           = 5.2148e-01
                           ,A0Phase          = 0. # fixed
                           ,AparPhase        = 3.2927e+00
                           ,f_S         =  0.
                           ,ASOddPhase  =  0.
                           ,dM               = 1.7629e+01
                           ,dGamma           = 9.3441e-02
                           ,Gamma            = 6.7722e-01
                           ,phiCP            = 0.07 # fixed to MC12Gen, not enough sensitivity
                           ,lambdaCP         = 1.   # fixed to Mc12Gen, not enough sensitivity
                            )

