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
    if kwargs.has_key('mkkRewData'):    data['mkkRewData'] = kwargs.pop('mkkRewData')

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
    obsCanv         = TCanvas('anglesTime_%s'%itNumb,'anglesTime_%s'%itNumb)
    assymObsCanv    = TCanvas('assymAnglesTime_%s'%itNumb,'assymAnglesTime_%s'%itNumb)
    KaonCanv        = TCanvas('KaonMomenta_%s'%itNumb,'KaonMomenta_%s'%itNumb)
    muonCanv        = TCanvas('muonMomenta_%s'%itNumb,'muonMomenta_%s'%itNumb)    
    assymKaonCanv   = TCanvas('assymKaonMomenta_%s'%itNumb,'assymKaonMomenta_%s'%itNumb)
    assymMuonCanv   = TCanvas('assymmuonMomenta_%s'%itNumb,'assymmuonMomenta_%s'%itNumb)
    KKMassCanv      = TCanvas('KKMass_%s'%itNumb,'KKMass_%s'%itNumb)
    assymKKMassCanv = TCanvas('assymKKMass_%s'%itNumb,'assymKKMass_%s'%itNumb)
    obsCanv.Divide(2,2)
    assymObsCanv.Divide(2,2)
    KaonCanv.Divide(4,2)
    muonCanv.Divide(4,2)
    assymKaonCanv.Divide(4,2)
    assymMuonCanv.Divide(4,2)
    
    # set some data drawing options
    colors      = dict( mcBefore = 2, mcAfter = kGreen+3, MomRewData = 4, Sdata = kMagenta+2, mkkRewData = 1 )
    stdDrawOpts = dict( DataError = RooAbsData.SumW2, MarkerSize = .6, XErrorSize = 0 )
    dataOpts    = dict()    
    for key in ['mcBefore','Sdata', 'mcAfter','MomRewData', 'mkkRewData']:
        if data.has_key(key): dataOpts[key] = dict( MarkerColor = colors[key], **stdDrawOpts  )  
            
    # plot angles and decay time
    print 'P2VV - INFO: compareDistributions: Plotting decay angles and time.'
    for canv, assymCanv, obs, logY in zip( 
        [obsCanv.cd(i+1) for i in range(len(observables))], 
        [assymObsCanv.cd(i+1) for i in range(len(observables))], 
        observables,
        3 * [False] + [True],
#        [ (-1.,1.), (-1.,1.), (-pi,pi), (0.,14.) ]
        ): 
        anglesFrames= compareDataSets( canv, obs, data = data, dataOpts = dataOpts, logy = logY,
                                       frameOpts = dict( Bins = 30 ), #, Range=rangeX ),
                                       )
        # make assymetry plots 
        makeAssymetryPlot(assymCanv, anglesFrames, referenceHistName ) 

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
    KKMassFrame = compareDataSets( KKMassCanv, KKMass, data = data, dataOpts = dataOpts, frameOpts = dict( Bins = 70 ))
    makeAssymetryPlot( assymKKMassCanv, KKMassFrame, referenceHistName )

    # make a legend and draw it
    legend, assym_legend = TPaveText( .47, .66, .77, .9, 'NDC' ), TPaveText( .269, .247, .569, .489, 'NDC' )
    for l in [legend,assym_legend]: l.SetFillColor(0)
    entriesNames = dict(mcBefore='mcBeforePhysRew', mcAfter='mcfterPhysRew',\
                            MomRewData='mcAfterMomRew', Sdata='data', mkkRewData = 'mcAfterMkkRew' )
    for key in ['mcBefore', 'mcAfter', 'MomRewData', 'Sdata', 'mkkRewData']:
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
    assymKKMassCanv.cd()
    assym_legend.Draw()
    _P2VVPlotStash.append(legend)
    _P2VVPlotStash.append(assym_legend)
    
    # print canvases in file
    for canv in [obsCanv,KaonCanv,muonCanv,assymKaonCanv,assymMuonCanv,KKMassCanv,assymObsCanv,assymKKMassCanv]: canv.Print(canv.GetName()+'.pdf')
    return [obsCanv,KaonCanv,muonCanv,assymKaonCanv,assymMuonCanv,KKMassCanv,assymObsCanv,assymKKMassCanv]

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
def TwoDimentionalVerticalReweighting(source, target, nbins, var, **kwargs):
    print 'P2VV - INFO: Initialised vertical reweigthing class, TwoDimentionalVerticalReweighting() for variables (%s,%s).'%(var[0],var[1])
    iterIdx  = kwargs.pop('iterationNumber' ,  0)
    plot     = kwargs.pop('xCheckPlots', False  )

    from ROOT import TH2F, TH1F, TCanvas
  
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
    
    # fill 2D rewweighting histograms
    for evnt in source: sourceHist.Fill( _valX(evnt), _valY(evnt), source.weight() )
    for evnt in target: targetHist.Fill( _valX(evnt), _valY(evnt), target.weight() )
       
    # rescale
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
    if plot: 
        weightsHist = TH1F('Weights', 'Weights',  2*nbins, .9*min(weights), 1.1*min(weights))
        for w in weights: weightsHist.Fill(w)
      
    # plot and print the 2d histograms
    if plot: 
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
    
        del source, target
        return weights, source
    else: 
        del source, target
        return weights

# function that reweighits a single source distribution to match a given target using a histogram
def OneDimentionalVerticalReweighting(source, target, nbins, var, **kwargs):
    print 'P2VV - INFO: Initialised vertical reweigthing class, OneDimentionalVerticalReweighting() for variable %s.'%var
    iterIdx  = kwargs.pop('iterationNumber' ,  0)
    plot     = kwargs.pop('xCheckPlots', False  )

    from ROOT import TH1F, TCanvas
    
    # dataset value getter 
    _valX = lambda ev: ev.find(var).getVal() 

    # get axis ranges
    sourceVar, targetVar = [],[]
    for event in source: sourceVar += [_valX(event)] 
    for event in target: targetVar += [_valX(event)]
    xMin, xMax = min(min(sourceVar),min(targetVar)), max(max(sourceVar),max(targetVar))
    for l in [ sourceVar, targetVar ]: del l
     
    # create 2D histrograms (Kplus_P vs Kminus_P)
    sourceHist  = TH1F('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins, xMin, xMax )
    targetHist  = TH1F('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins, xMin, xMax )
  
    # fill reweighting histograms
    for evnt in source: sourceHist.Fill( _valX(evnt), source.weight() )
    for evnt in target: targetHist.Fill( _valX(evnt), target.weight() )
       
    # rescale
    if source.numEntries() > target.numEntries(): sourceHist.Scale( target.sumEntries() / source.sumEntries() )
    else: targetHist.Scale( source.sumEntries() / target.sumEntries() )
    
    # calculate weights
    weights = []
    count = 0 # count how many events have a problematic weight
    for event in source:
        bin = sourceHist.FindFixBin( _valX(event) ) # get the bin with given var value
        if targetHist.GetBinContent(bin)==0 or sourceHist.GetBinContent(bin)==0:# do not weight the event if not possible with current binning
            weights += [1] 
            count += 1
        else: 
            weights += [targetHist.GetBinContent(bin) / sourceHist.GetBinContent(bin)] # calculate weight
    if count>0: print 'P2VV - INFO: OneDimentionalVerticalReweighting: %s out of %s events are not weighted.'%(count,source.numEntries())

    if plot: # check the result of the reweighting 
        weightsHist = TH1F('Weights', 'Weights',  3*nbins, .9*min(weights), 1.1*min(weights))
        for w in weights: weightsHist.Fill(w)

        test_s = TH1F('test_s','test_s',3*nbins, xMin,xMax)
        test_t = TH1F('test_t','test_t',3*nbins, xMin,xMax)
        sourceEvtList, targetEvtList = [],[]
        for ev in source: sourceEvtList+=[ _valX(ev) ]
        for ev in target: targetEvtList+=[ _valX(ev) ]
        for ev, w in zip(sourceEvtList,weights): test_s.Fill(ev,w)
        for evnt in target: test_t.Fill(_valX(evnt), target.weight())
        # plot and print the histograms
        testCanv = TCanvas('test','test')
        test_t.Draw()
        test_s.Scale(target.sumEntries() / source.sumEntries())
        test_s.Draw('same err')
       
        from P2VV.Utilities.Plotting import _P2VVPlotStash
        _P2VVPlotStash += [test_s,test_t]
        testCanv.Print( testCanv.GetName() + '_%s.pdf'%iterIdx )
        
        del source, target
        return weights, testCanv
    else: 
        del source, target
        return weights
 
class WeightedDataSetsManager(dict):
    def __init__( self, **kwargs ):       
        print 'P2VV - INFO: Initialised datasets manager class WeightedDataSetsManager().'
        self['initSource'] = kwargs.pop('source', '')
        self['dataSets']   = dict( initSource = self['initSource'] )
        
        self['WeightsLists']    = {}
        self['combinedWeights'] = []   
        self['latestDataSetPointer'] = 'initSource'
        
        self['iterationNumber'] = 0
        self['saveIntermediateDatasets'] = False
        
    def appendWeights( self, weightsName, weightsList ):
        from numpy import array
        self['WeightsLists'][weightsName] = array( weightsList )
        
        # combine weights
        if not len( self['WeightsLists'] ) <= 1:
            print 'P2VV - INFO: WeightedDataSetsManager: Combining weights, named %s, with existing ones.'%weightsName
            self['combinedWeights'] = 1
            for wList in self['WeightsLists'].keys(): self['combinedWeights'] *= array( self['WeightsLists'][wList] )
        else: self['combinedWeights'] = self['WeightsLists'][weightsName]

        # scale weights to preserve number of events 
        print 'P2VV - INFO: WeightedDataSetsManager: Scaling sources sum of weights to the number of entries.'
        n_events = self['initSource'].numEntries()
        sumW = sum( self['combinedWeights'] )
        self['combinedWeights'] =  ( n_events / sumW ) * self['combinedWeights'] 
        
        # write weights
        wName = ''
        for name in self['WeightsLists'].keys(): wName += name + '_'
        wName += str(self['iterationNumber'])
        self['dataSets'][weightsName] = self.writeWeights(self['dataSets'][self['latestDataSetPointer']], 
                                                          'weight_' + wName, wName, 
                                                          )

        # delete intermediate dataset, except if it is the initial source and or you want to keep them for plotting 
        if not self['saveIntermediateDatasets'] and not self['latestDataSetPointer'] == 'initSource':
            print 'P2VV - INFO: appendWeights: Deleting dataset named ' + self['dataSets'][self['latestDataSetPointer']].GetName()
            del self['dataSets'][self['latestDataSetPointer']]

        # bookkeeping flag 
        self['latestDataSetPointer'] = weightsName
        self['combinedWeightsName']  = 'weight_' + wName
        
    def writeWeights( self, dataset, weightsName, writeDatasetName ):
        print 'P2VV - INFO: writeWeights: Creating dataset with name %s and weight name %s:'%(writeDatasetName,weightsName)
        from ROOT import RooArgSet, RooRealVar, RooDataSet
        weightsVar = RooRealVar( weightsName, weightsName, 1, .9*min(self['combinedWeights']), 1.1*max(self['combinedWeights']) )
        weightsArgSet  = RooArgSet( weightsVar )
        weightsDataSet = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        for weight in self['combinedWeights']:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        _dataset  = RooDataSet( writeDatasetName, writeDatasetName, dataset.get(), Import = dataset )
        _dataset.merge( weightsDataSet )
        _Wdataset = RooDataSet( writeDatasetName, writeDatasetName, _dataset.get(), Import = _dataset, WeightVar = (weightsName,True) )

        del weightsDataSet, dataset, _dataset
        return _Wdataset

    def getDataSet( self, which='' ): 
        dataSetKey = which if which else self['latestDataSetPointer']
        return self['dataSets'][dataSetKey]

    def plotWeights(self, which = '', Range=() ):
        from ROOT import TCanvas, TH1F
        from P2VV.Utilities.Plotting import _P2VVPlotStash

        which = which if which else self['latestDataSetPointer']
        c = TCanvas( 'weights_' + which + str(self['iterationNumber']), 'weights_' + which + str( self['iterationNumber']) )        
        plotRange = min(self['WeightsLists'][which]), max(self['WeightsLists'][which]) if not Range else Range    
        weghtsHist = TH1F('weights_'+which, 'weights_'+which, 200, plotRange[0], plotRange[1])
        for weight in self['WeightsLists'][which]: weghtsHist.Fill(weight)
        c.cd()
        weghtsHist.SetMarkerSize(.5)
        weghtsHist.Draw('err')
        c.Print( 'weights_%s_.pdf'%which )
        _P2VVPlotStash +=[c,weghtsHist]

    def clear( self ): 
        del self['dataSets'], self['WeightsLists'],  self['combinedWeights']
        self['dataSets'] =  dict( initSource = self['initSource'] )
        self['WeightsLists'] = dict()
        self['latestDataSetPointer'] = 'initSource'

# Class for multipling a pdf with an angular acceptance and performs an sFit 
class BuildBs2JpsiKKFit():
    def __init__( self,**kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: BuildBs2JpsiKKFit().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
        
        # specify running period
        MCProd = kwargs.pop('MonteCarloProduction', 'Sim08')
        if   '2011' in MCProd: runPeriod = '2011'
        elif '2012' in MCProd: runPeriod = '2012'
        else: runPeriod = '3fb'
        self._pdfConfig = PdfConfig( RunPeriods = runPeriod )
        
        # blind / unblind (phi_s,dGama)
        #self._pdfConfig['blind'] = {}
        
        # give parameters of the pdf a prefix
        self._pdfConfig['parNamePrefix'] = 'data' 
        namePF = self._pdfConfig['parNamePrefix']

        self._dataSetPath = kwargs.pop('dataSetPath', None)
        self._dataSetName = kwargs.pop('dataSetName', None)
        self._weightsName = kwargs.pop('weightsName', 'JpsiKK_sigSWeight')
    
        self._doUntaggedFit = kwargs.pop('doUntaggedFit', '')        
        self._doNullTest    = kwargs.pop('doNullTest',    '')
     
        # fit options
        fitOpts = dict( NumCPU = 2, Optimize  = 2, Minimizer = 'Minuit2' )
        self._pdfConfig['fitOptions'] = fitOpts
        corrSFitErrCats         = [ 'runPeriod', 'KKMassCat' ] if ('2011' not in MCProd or '2012' not in MCProd) else [ 'KKMassCat' ]
        randomParVals           = ( ) # ( 1., 12345 )

        # PDF options
        # time acceptance
        if runPeriod=='3fb':
            self._pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2011' ) ] )['file']\
                = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
            self._pdfConfig['timeEffHistFiles'].getSettings( [ ( 'runPeriod', 'p2012' ) ] )['file']\
                = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
        elif runPeriod=='2011':
            self._pdfConfig['timeEffHistFiles'] = dict(  
                file  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
                , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                )
        elif runPeriod=='2012':
            self._pdfConfig['timeEffHistFiles'] = dict(  
                file = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
                , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
                , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                )
            
        # angular acceptance
        self._pdfConfig['angEffMomsFiles'] = ''
        self._pdfConfig['anglesEffType']   = ''
        
        # read data set from file
        from P2VV.Utilities.DataHandling import readData
        dataSet = readData( filePath = self._dataSetPath, dataSetName = self._dataSetName,  NTuple = False )
        self._pdfConfig['signalData'] = dataSet
        self._pdfConfig['readFromWS'] = True
        
        # build the PDF
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
        self._pdfBuild = PdfBuilder( **self._pdfConfig )
        self._pdf = self._pdfBuild.pdf()
        
        # data set with weights corrected for background dilution: for phi_s fit only!
        from P2VV.Utilities.DataHandling import correctWeights
        self._fitData = correctWeights( dataSet, corrSFitErrCats )

        # fix values of some parameters
        for CEvenOdds in self._pdfBuild['taggingParams']['CEvenOdds'] :
            if not self._pdfConfig['SSTagging'] :
                if namePF:
                    CEvenOdds.setConstant( self._pdfConfig['parNamePrefix'] + '_' + 'avgCEven.*')
                    CEvenOdds.setConstant( self._pdfConfig['parNamePrefix'] + '_' + 'avgCOdd.*', True )
                else: 
                    CEvenOdds.setConstant( 'avgCEven.*')
                    CEvenOdds.setConstant( 'avgCOdd.*', True )
            else :
                for CEvenOdd in CEvenOdds :
                    if namePF:
                        CEvenOdd.setConstant( self._pdfConfig['parNamePrefix'] + '_' + 'avgCEven.*')
                        CEvenOdd.setConstant( self._pdfConfig['parNamePrefix'] + '_' + 'avgCOdd.*', True )
                    else: 
                        CEvenOdd.setConstant('avgCEven.*')
                        CEvenOdd.setConstant( 'avgCOdd.*', True )
        for par in self._pdf.getParameters(self._fitData):
            if 'C_SP' in par: par.setConstant()
            #self._pdfBuild['amplitudes'].setConstant('C_SP')

        # print parameters
        print 120 * '='
        print 'Bs2JpsiKKFit: fit data:'
        self._fitData.Print()
        print 'Bs2JpsiKKFit: observables in PDF:'
        self._pdf.getObservables(self._fitData).Print('v')
        print 'Bs2JpsiKKFit: parameters in PDF:'
        self._pdf.getParameters(self._fitData).Print('v')
        print 'Bs2JpsiKKFit: constraints in PDF:'
        for constr in self._pdf.ExternalConstraints() : constr.Print()

        self._FitResults  = {} # collect all the fit results
        self._Moments     = {} # collect all moments
       
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
            if par.startswith('__'):
                fitResultParKey = '__%s_'%self._pdfConfig['parNamePrefix'] + par.partition('__')[2]
            else:
                fitResultParKey =  self._pdfConfig['parNamePrefix'] + '_' + par
            try: oldPars[par] = fitResult.floatParsFinal().find(fitResultParKey).getVal()
            except AttributeError: 
                print 'P2VV - WARNING: updateDataParameters: Parameter %s not found in fit result %s, parameter value will not change.'\
                    %( fitResultParKey, fitResult.GetName() )
        print 'P2VV - INFO: updateDataParameters: Updating physics parameters from sFit.'
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
    def getBlindString(self):       return self._pdfConfig['blind']
    def getObservables(self,which=None):
        if   which=='angles': return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName().startswith('hel') ]
        elif which=='time':   return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName()==('time')         ]
        else:                 return             self._pdfBuild['obsSetP2VV']

# Vertical reweighting class to match physics of weighted distribution using a pdf
class MatchPhysics( ):
    def __init__( self, nTupleFile, nTupleName, **kwargs ):      
        # monte carlo gen conditions specifier
        MCProd = kwargs.pop('MonteCarloProduction', 'Sim08')        
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'        
        print 'P2VV - INFO: Matching physics on %s mc sample.'%MCProd

        # set global object name prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        self._namePF = 'mc'
        setParNamePrefix( self._namePF )
        
        # blind parameters in MC pdf
        blindStr = kwargs.pop( 'BlindPdf', {} )

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
        time      = _createGetObservable('time')
        angles    = [_createGetObservable(o) for o in ['helcosthetaK','helcosthetaL','helphi'] ]
        iTag      = _createGetObservable('iTag')
        KKMass = _createGetObservable('KKMass')
        KKMassCat = RooObject._rooobject('KKMassCat')
        runPeriod = RooObject._rooobject('runPeriod')
        self._obsSet = [ trueTime, time, KKMass ] + angles + [ runPeriod, KKMassCat ]
        
        # add track momenta and B momenta
        self._obsSet += [ RooObject._rooobject('%s_%s'%( part, comp )) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'PX', 'PY', 'PZ', 'P' ) ]
        self._obsSet += [ RooObject._rooobject('B_%s'%comp) for comp in [ 'P', 'Pt' ] ]

        # read ntuple
        from P2VV.Utilities.DataHandling import readData
        if   MCProd == 'Sim08_2011':         cuts = 'runPeriod==2011'
        elif MCProd == 'Sim08_2012':         cuts = 'runPeriod==2012'
        elif MCProd == 'Sim08_2011_reduced': cuts = 'runPeriod==2011 && runNumber>2543.93e3 && runNumber<2544e3' # 2543.87 # 2542e3
        elif MCProd == 'Sim08_2012_reduced': cuts = 'runPeriod==2012 && runNumber>2523e3 && runNumber<2525.35e3'
        elif MCProd == 'Sim08_reduced':      cuts = '(runNumber<2524e3) || (runNumber>2546e3)'
        else: cuts = ''
        self._data = readData( nTupleFile, dataSetName=nTupleName, NTuple=True, observables=self._obsSet, ntupleCuts=cuts)

        # build pdf
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
        amplitudes = Amplitudes(  A0Mag2 = A0Mag2Val, AperpMag2 = AperpMag2Val
                                  , AparPhase = AparPhVal - A0PhVal, AperpPhase = AperpPhVal - A0PhVal
                                  , f_S = 0., ASOddPhase = 0., C_SP = 1.
                                  )

        # B lifetime
        if blindStr and 'dGamma' in blindStr: dGammaVar = dict( Value = dGammaVal, Blind = blindStr['dGamma']  )
        else: dGammaVar = dict( Value = dGammaVal )
        from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVar, dM = dMVal )
        
        # resolution model
        tResArgs = { }
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = trueTime
        timeResModel = TimeResolution( **tResArgs )

        # CP violation parameters
        if blindStr and 'phiCP' in blindStr: phiCPVar = dict( Value = phiCPVal, Blind = blindStr['phiCP']  )
        else: phiCPVar = dict( Value = phiCPVal )
        from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
        lambdaCP = CPParam( lambdaCP = lambCPVal, phiCP = phiCPVar )
        
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
        if 'Sim08' in MCProd: from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as mcPars
        if 'Sim06' in MCProd: from P2VV.Utilities.MCReweighting import parValuesMcSim06_6KKmassBins as mcPars
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

        #del self._data
    
    def setMonteCarloParameters(self, pars=None):
        print 'P2VV - INFO: setMonteCarloParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        for key in self._pdfPhysPars.keys():
            if key.startswith('__'): self._pdfPhysPars[key].setVal( self._mcPhysParsSet['_' + key.partition('__mc')[2]]  )
            else:                    self._pdfPhysPars[key].setVal( self._mcPhysParsSet[      key.partition('_')[2]]     ) 
        for k in self._pdfPhysPars.keys(): 
            print '%20s %.4f'%(self._pdfPhysPars[k].GetName(), self._pdfPhysPars[k].getVal())
                    
    def setDataFitParameters(self, dataPars, KKmassCat=None):
        print 'P2VV - INFO: setDataFitParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        for key in self._pdfPhysPars.keys():
            if key.startswith('__'): self._pdfPhysPars[key].setVal( dataPars[ '_' + key.partition('__mc')[2]]  ) 
            else:                    self._pdfPhysPars[key].setVal( dataPars[ key.partition('_')[2]] ) 
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
            
        # loop over events and calculate physics weights
        self._pdf.attachDataSet( self._data ) # make the pdf directly dependant on data
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'    
        self.setMonteCarloParameters()
        for nev in xrange(self._data.numEntries()):
            self._data.get(nev)
            denominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: calculatePhysicsWeights: Calculating nominators for phyisics matching weights'
        self.setDataFitParameters(dataParameters)
        for nev in xrange(self._data.numEntries()):
            self._data.get(nev)
            nominators.append( self._pdf.getVal(normVars) )
        print 'P2VV - INFO: Calculating phyisics matching weights'
        calculatePhysicsWeights(nom=nominators, den=denominators)
        
        return self._physWeights
 
    def getDataSet(self):           return self._data
    def getPdf(self):               return self._pdf
    def getAngleFunctions(self):    return self._angleFuncs
    def getMcTime(self):            return  [ o for o in self._pdf.Observables() if 'time' in o.GetName() ]
    def getParNamePrefix(self):     return self._namePF   

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

# sFit in 2011 + 2012 data  
parValues6KKmassBins20112012 = dict(
    __dGamma__	       = 0.086993
    ,__phiCP__	       = 0.0830558
    ,A0Mag2	       = 0.5208
    ,ASOddPhase_bin0   = 0.803468
    ,ASOddPhase_bin1   = 2.3092
    ,ASOddPhase_bin2   = 0.460888
    ,ASOddPhase_bin3   = -0.360896
    ,ASOddPhase_bin4   = -0.679218
    ,ASOddPhase_bin5   = -0.824088
    ,AparPhase	       = 3.25566
    ,AperpMag2	       = 0.253554
    ,AperpPhase	       = 3.18827
    ,Gamma	       = 0.660572
    ,betaTimeEff_p2011 = -0.00832011
    ,betaTimeEff_p2012 = -0.0137996
    ,dM	               = 17.7651
    ,f_S_bin0	       = 0.444554
    ,f_S_bin1	       = 0.0632622
    ,f_S_bin2	       = 0.00864546
    ,f_S_bin3	       = 0.00887104
    ,f_S_bin4	       = 0.0443292
    ,f_S_bin5	       = 0.210105
    ,lambdaCP	       = 0.963988
    ,wTagDelP0OS       = 0.0137251
    ,wTagDelP0SS       = -0.0159546
    ,wTagDelP1OS       = 0.0697653
    ,wTagDelP1SS       = 0.0150089
    ,wTagP0OS	       = 0.39053
    ,wTagP0SS	       = 0.440126
    ,wTagP1OS	       = 1.03623
    ,wTagP1SS	       = 0.944176
    )

 # sFit in 2011 data  
parValues6KKmassBins2011 = dict(
     A0Mag2          = 0.530647
    ,A0Phase         = 0.
    ,ASOddPhase_bin0 = 0.646488	
    ,ASOddPhase_bin1 = 0.970226	
    ,ASOddPhase_bin2 = 0.424538	
    ,ASOddPhase_bin3 = -0.394219	
    ,ASOddPhase_bin4 = -0.854113	
    ,ASOddPhase_bin5 = -0.607166	
    ,AparPhase	     = 3.32684	
    ,AperpMag2	     = 0.246338	
    ,AperpPhase	     = 3.39293	
    ,Gamma	     = 0.657435	
    ,__dGamma__	     = 0.122035	
    ,__phiCP__	     = -0.133016	
    ,betaTimeEff     = -0.00830863	
    ,dM	             = 17.7666	 
    ,f_S_bin0	     = 0.39986	
    ,f_S_bin1	     = 0.0538139	
    ,f_S_bin2	     = 0.0112028	
    ,f_S_bin3	     = 0.0149245	
    ,f_S_bin4	     = 0.0247391	
    ,f_S_bin5	     = 0.203116	
    ,lambdaCP	     = 0.951561	
    ,wTagDelP0OS     = 0.0137092	
    ,wTagDelP0SS     = -0.0159789	
    ,wTagDelP1OS     = 0.0698641	
    ,wTagDelP1SS     = 0.0149033	
    ,wTagP0OS	     = 0.379543	
    ,wTagP0SS	     = 0.435735	
    ,wTagP1OS	     = 1.02522	
    ,wTagP1SS	     = 0.979372	
     )

 # sFit in 2011 reduced data sample  
parValues6KKmassBins2011_reduced = dict(
    A0Mag2	     = 0.530141
    ,A0Phase         = 0.
    ,ASOddPhase_bin0 = 0.64447
    ,ASOddPhase_bin1 = 0.958849
    ,ASOddPhase_bin2 = 0.411171
    ,ASOddPhase_bin3 = -0.395949
    ,ASOddPhase_bin4 = -0.859725
    ,ASOddPhase_bin5 = -0.599493
    ,AparPhase	     = 3.31748
    ,AperpMag2	     = 0.245888
    ,AperpPhase	     = 3.39694
    ,Gamma	     = 0.657411
    ,__dGamma__	     = 0.12224
    ,__phiCP__	     = -0.13062
    ,betaTimeEff     = -0.00830066
    ,dM	             = 17.7666
    ,f_S_bin0	     = 0.400578
    ,f_S_bin1	     = 0.0549151
    ,f_S_bin2	     = 0.0122189
    ,f_S_bin3	     = 0.014571
    ,f_S_bin4	     = 0.0244505
    ,f_S_bin5	     = 0.204648
    ,lambdaCP	     = 0.962228
    ,wTagDelP0OS     = 0.0137092
    ,wTagDelP0SS     = -0.0159797
    ,wTagDelP1OS     = 0.0698618
    ,wTagDelP1SS     = 0.0149051
    ,wTagP0OS	     = 0.379598
    ,wTagP0SS	     = 0.435766
    ,wTagP1OS	     = 1.02458
    ,wTagP1SS	     = 0.978629
    )

 # sFit in 2012 data  
parValues6KKmassBins2012 = dict(
     A0Mag2	      = 0.516663
    ,A0Phase         = 0.
    ,ASOddPhase_bin0 = 0.880057
    ,ASOddPhase_bin1 = 2.40856
    ,ASOddPhase_bin2 = 0.691195
    ,ASOddPhase_bin3 = -0.447733
    ,ASOddPhase_bin4 = -0.615622
    ,ASOddPhase_bin5 = -0.902559
    ,AparPhase	     = 3.15983
    ,AperpMag2	     = 0.257204
    ,AperpPhase	     = 3.02308
    ,Gamma	     = 0.675951
    ,__dGamma__       = 0.0600561
    ,__phiCP__       = -0.121879
    ,dM	             = 17.7657
    ,f_S_bin0	     = 0.459275
    ,f_S_bin1        = 0.0712949
    ,f_S_bin2	     = 0.00404748
    ,f_S_bin3	     = 0.00334672
    ,f_S_bin4	     = 0.0579512
    ,f_S_bin5	     = 0.214111
    ,lambdaCP	     = 0.981036
    ,wTagDelP0OS     = 0.0137176
    ,wTagDelP0SS     = -0.0159754
    ,wTagDelP1OS     = 0.0698766
    ,wTagDelP1SS     = 0.0151004
    ,wTagP0OS        = 0.388797
    ,wTagP0SS	     = 0.440893
    ,wTagP1OS	     = 1.01434
    ,wTagP1SS	     = 0.964957
)

# sFit in 2011 data no KK mass bins
parValuesNoKKBinsWideKKWindow = dict(  
    AperpMag2        = 0.24663
    ,AperpPhase       = 3.0232
    ,A0Mag2           = 0.52352
    ,A0Phase          = 0 # cosntrained
    ,AparPhase        = 3.2144
    ,f_S              = 0.0463
    ,ASOddPhase       = -0.0666
    ,dM               = 17.6739
    ,__dGamma__       = 0.1021
    ,Gamma            = 0.6728
    ,__phiCP__        = 0.0847
    ,lambdaCP         = 0.9275
    )

# MC generating conditions  
parValuesMcSim08_6KKmassBins = dict(  
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
    ,__dGamma__        = 1. / 1.406 - 1. / 1.614
    ,__phiCP__         = +0.07 
    ,dM                = 17.8
    ,f_S_bin0          = 0.
    ,f_S_bin1          = 0.
    ,f_S_bin2          = 0.
    ,f_S_bin3          = 0.
    ,f_S_bin4          = 0.
    ,f_S_bin5          = 0.
    ,lambdaCP          = 1.
    )

# Sim06 generating conditions  
parValuesMcSim06_6KKmassBins = dict(
     A0Mag2            = 0.6
    ,A0Phase           = 0. 
    ,ASOddPhase_bin0   = 0.
    ,ASOddPhase_bin1   = 0.
    ,ASOddPhase_bin2   = 0.
    ,ASOddPhase_bin3   = 0.
    ,ASOddPhase_bin4   = 0.
    ,ASOddPhase_bin5   = 0.
    ,AparPhase         = 2.5
    ,AperpMag2         = 0.16
    ,AperpPhase        = -0.17
    ,Gamma             = 0.679
    ,__dGamma__        = 0.06
    ,__phiCP__         = -0.04
    ,dM                = 17.8
    ,f_S_bin0          = 0.
    ,f_S_bin1          = 0.
    ,f_S_bin2          = 0.
    ,f_S_bin3          = 0.
    ,f_S_bin4          = 0.
    ,f_S_bin5          = 0.
    ,lambdaCP          = 1.
    )

# For null test
parValuesMc2012Fit =  dict(
    # sFit on '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting_part1.root'
    #  AperpMag2        =  2.5628e-01 
    # ,AperpPhase       = -2.8377e+00 
    # ,A0Mag2           =  5.1761e-01
    # ,A0Phase          =  0. # cosntrained
    # ,AparPhase        =  3.2728e+00 
    # ,f_S_bin0         =  0.
    # ,ASOddPhase_bin0  =  0.
    # ,dM               =  1.7630e+01 
    # ,__dGamma__           =  8.0522e-02 
    # ,Gamma            =  6.8023e-01
    # ,__phiCP__            = -1.2965e-01
    # ,lambdaCP         =  1.0684e+00
    
    # sFit on '/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting.root'
    # using 2012 angular acceptance weights
     AperpMag2        = 2.4841e-01
    ,AperpPhase       = 3.4126e+00
    ,A0Mag2           = 5.2148e-01
    ,A0Phase          = 0. # fixed
    ,AparPhase        = 3.2927e+00
    ,f_S         =  0.
    ,ASOddPhase  =  0.
    ,dM               = 1.7629e+01
    ,__dGamma__           = 9.3441e-02
    ,Gamma            = 6.7722e-01
    ,__phiCP__            = 0.07 # fixed to MC12Gen, not enough sensitivity
    ,lambdaCP         = 1.   # fixed to Mc12Gen, not enough sensitivity
    )
