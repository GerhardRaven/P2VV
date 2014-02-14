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
    Legend          = kwargs.pop('legend',         False)
    
    # get datasets
    from ROOT import RooFit, RooDataSet
    data = dict( Sdata = kwargs.pop('sData') )
    for key in [ 'mcData', 'mkkRewData','BmomRewData', 'mcDataPhysRew', 'MomRewData', 'BmomMkkRewData' ]:
        if kwargs.has_key(key):
            if kwargs[key]: data[key] = kwargs.pop(key)
    names = ''
    for d in data.itervalues(): names +=  d.GetName() + '\n'
    print 'P2VV - INFO: compareDistributions: Making comparition plots for datasets:', names
    
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
        elif obsName == 'B_P':   B_P = obs
        elif obsName == 'B_Pt':   B_Pt = obs
        
    # assymetry plots are compared w.r.t. the sData
    referenceHistName = 'h_' + data['Sdata'].GetName() 

    # start drawing
    from ROOT import TCanvas, RooAbsData, TPaveText, kGreen, kMagenta, kOrange
    from P2VV.Utilities.Plotting import compareDataSets, makeAssymetryPlot
    from P2VV.Load import LHCbStyle
    from math import pi

    # make canvases
    namePF          = '' 
    for name in data.keys(): namePF += name + '_' 
    namePF += kwargs.pop('prodData')
    obsCanv         = TCanvas('anglesTime_%s_'%itNumb + namePF       ,'anglesTime_%s'%itNumb)
    assymObsCanv    = TCanvas('assymAnglesTime_%s_'%itNumb + namePF  ,'assymAnglesTime_%s'%itNumb)
    KaonCanv        = TCanvas('KaonMomenta_%s_'%itNumb + namePF      ,'KaonMomenta_%s'%itNumb)
    muonCanv        = TCanvas('muonMomenta_%s_'%itNumb + namePF      ,'muonMomenta_%s'%itNumb)    
    assymKaonCanv   = TCanvas('assymKaonMomenta_%s_'%itNumb + namePF ,'assymKaonMomenta_%s'%itNumb)
    assymMuonCanv   = TCanvas('assymmuonMomenta_%s_'%itNumb + namePF ,'assymmuonMomenta_%s'%itNumb)
    KKMassCanv      = TCanvas('KKMass_%s_'%itNumb + namePF           ,'KKMass_%s'%itNumb)
    BmomCanv        = TCanvas('B_P_%s_'%itNumb + namePF              ,'B_P_%s'%itNumb)
    obsCanv.Divide(2,2)
    assymObsCanv.Divide(2,2)
    muonCanv.Divide(4,2)
    assymKaonCanv.Divide(4,2)
    KaonCanv.Divide(4,2)
    assymMuonCanv.Divide(4,2)
    BmomCanv.Divide(2,2)
    KKMassCanv.Divide(2,2)

    # set some data drawing options
    colors      = dict(mcData=2, mcDataPhysRew=kGreen+3, MomRewData=4, Sdata=kMagenta+2, 
                       mkkRewData=1, BmomRewData=5, BmomMkkRewData = kOrange+8  )
    stdDrawOpts = dict( DataError = RooAbsData.SumW2, MarkerSize = .6, XErrorSize = 0. )
    dataOpts    = dict()    
    for key in colors.keys():
        if data.has_key(key): dataOpts[key] = dict( MarkerColor = colors[key], **stdDrawOpts  )
        
    # plot angles and decay time
    for canv, assymCanv, obs, logY, assymYrange in zip( 
        [obsCanv.cd(i+1) for i in xrange(len(observables))], 
        [assymObsCanv.cd(i+1) for i in xrange(len(observables))], 
        observables,
        3 * [False] + [True],
        3 * [(-.1,.1),] + [(-3,3)]
        ):
        print 'P2VV - INFO: compareDistributions: Plotting %s.'%obs.GetName()
        anglesFrames= compareDataSets( canv, obs, data = data, dataOpts = dataOpts, logy = logY,
                                       frameOpts = dict( Bins = 30 ), 
                                       )
        # make assymetry plots
        # print 'P2VV - INFO: compareDistributions: Creating assymentry plot for %s'%obs.GetName()
        makeAssymetryPlot(assymCanv, anglesFrames, referenceHistName, len(data.keys()), yRange=assymYrange ) 

    # plot Kaon and muon momenta
    print 'P2VV - INFO: compareDistributions: Plotting track momenta.'
    for canv, assymCanv, obs, assymYrange in zip( 
        [ KaonCanv.cd(k+1) for k in xrange(len(Kmomenta)) ]      + [ muonCanv.cd(m+1) for m in xrange(len(muMomenta)) ],
        [ assymKaonCanv.cd(k+1) for k in xrange(len(Kmomenta)) ] + [ assymMuonCanv.cd(m+1) for m in xrange(len(muMomenta)) ],
        Kmomenta + muMomenta,
        16 * [(-.4,.4),] 
        ):
        print 'P2VV - INFO: compareDistributions: Plotting %s.'%obs.GetName()
        momFrame = compareDataSets( canv, obs, data = data, dataOpts = dataOpts,
                                    frameOpts = dict( Bins = 30, Range=trackMomRangeX[obs.GetName()] )
                                    )
        # make assymetry plots
        # print 'P2VV - INFO: compareDistributions: Creating assymentry plot for %s'%obs.GetName()
        makeAssymetryPlot( assymCanv, momFrame, referenceHistName, len(data.keys()), yRange=assymYrange ) 

    # plot KKMass
    print 'P2VV - INFO: compareDistributions: Plotting KKmass.'
    KKMassFrame = compareDataSets( KKMassCanv.cd(1), KKMass, data = data, dataOpts = dataOpts, frameOpts = dict( Bins = 40 ))
    # print 'P2VV - INFO: compareDistributions: Creating assymentry plot for KKmass'
    makeAssymetryPlot( KKMassCanv.cd(2), KKMassFrame, referenceHistName, len(data.keys()) )

    # plot B_P and B_Pt
    for pad, obs, rangeY in zip( [1,2], [ B_P,B_Pt], [(0.,5e5),(0.,3e4)] ):
        print 'P2VV - INFO: compareDistributions: Plotting %s.'%obs.GetName()
        BmomFrame = compareDataSets( BmomCanv.cd(pad), obs, data = data, dataOpts = dataOpts, frameOpts = dict( Bins=70, Range=rangeY ))
        # print 'P2VV - INFO: compareDistributions: Creating assymentry plot for %s'%obs.GetName()
        makeAssymetryPlot( BmomCanv.cd(pad+2), BmomFrame, referenceHistName, len(data.keys()) )

    # make a legend and draw it
    if Legend:
        legend, assym_legend = TPaveText(.587, .66, .893, .902, 'NDC' ), TPaveText( .269, .247, .569, .489, 'NDC' )
        for l in [legend,assym_legend]: l.SetFillColor(0)
        entriesNames = dict(mcData='mcNominal', mcDataPhysRew='mcPhysRew', mkkRewData = 'mcMkkRew', \
                                MomRewData='mKKMomRew', Sdata='data', BmomRewData = 'mcBmomRew' )
        for key in ['mcData', 'BmomRewData', 'mcDataPhysRew', 'mkkRewData', 'MomRewData', 'Sdata']:
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
        KKMassCanv.cd(1)
        legend.Draw()
        BmomCanv.cd(1)
        legend.Draw()
        _P2VVPlotStash.append(legend)
        _P2VVPlotStash.append(assym_legend)
        
    # print canvases in file
    for canv in [obsCanv,KaonCanv,muonCanv,assymObsCanv,assymMuonCanv,assymKaonCanv,KKMassCanv,BmomCanv]: canv.Print(canv.GetName() +  '.pdf')
    return [obsCanv,KaonCanv,muonCanv,assymObsCanv,assymMuonCanv,assymKaonCanv,KKMassCanv,BmomCanv]

# clean P2VVPlotStash to save memory
def cleanP2VVPlotStash():
    from P2VV.Utilities.Plotting import _P2VVPlotStash
    while len(_P2VVPlotStash) > 0: 
        for plot in _P2VVPlotStash: _P2VVPlotStash.remove(plot)
    print 'P2VV - INFO: cleanP2VVPlotStash: Emptied P2VVplotStash.'
    
# easily create an observable using information from the PdfConfig class
def _createGetObservable(name):
    from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis 
    PdfConfig = Bs2Jpsiphi_RunIAnalysis( RunPeriods = '3fb' )
    obsDict = PdfConfig['obsDict']
    
    if name == 'helcosthetaK': name = 'cpsi'
    if name == 'helcosthetaL': name = 'ctheta'
    if name == 'helphi'      : name = 'phi'
    if name == 'iTag'        : obsDict['iTag']       = ( 'iTag',     'Initial state flavour tag', { 'Untagged' : 0 } ) 
    if name == 'truetime'    : obsDict['truetime']   = ( 'truetime', 'true time', 'ps',            1.5,  0.,  30.   )
    if name == 'KKMassCat'   : obsDict['KKMassCat']  = ( 'KKMassCat', 'KK mass category', dict( [ ('bin%d' % i, i) for i in range(6) ] )                                 )

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

# combine moments function with waiting feature.
def combineMoments( accFile1, accFile2, outName, Prefix='', delete=True):
    from time import sleep
    import os
    from P2VV.Utilities.DataMoments import combineMoments
    while True:
        try:
            open(accFile1)
            open(accFile2)
            combineMoments( [ accFile1, accFile2 ], outName, prefix = Prefix, printMoms = False )
            if delete: 
                for angAcc in [accFile1,accFile2]: os.remove(angAcc)
            break
        except IOError:
            print 'P2VV - INFO: combineMoments: Waiting for the followig flies to combine efficiency moments:\n  %s\n  %s'%(accFile1,accFile2)
            sleep(60)

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

# function that reweighits the KK distributions with a 2D histrogram
def TwoDimentionalVerticalReweighting(source, target, nbins, var, **kwargs):
    print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Reweighting variables (%s, %s) in sample %s.'%(var[0],var[1],source.GetName())
    iterIdx        = kwargs.pop('iterationNumber',   0   )
    plot           = kwargs.pop('xCheckPlots',     False )
    equalStatsBins = kwargs.pop('equalStatBins',   False )
    combineWeights = kwargs.pop('combWeights',     False )

    from ROOT import TH2D, TH1D, TCanvas
  
    _valX = lambda ev: ev.find(var[0]).getVal() # value getter of the first variable
    _valY = lambda ev: ev.find(var[1]).getVal() # value getter of the first variable
    
    # get axis ranges
    from P2VV.RooFitWrappers import RooObject
    xMin, yMin = RooObject._rooobject(var[0]).getMin(), RooObject._rooobject(var[1]).getMin()
    xMax, yMax = RooObject._rooobject(var[0]).getMax(), RooObject._rooobject(var[1]).getMax()
        
    # import / create binning
    print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Using %s binnning with %sx%s bins.'%('almost equal statistics' if equalStatsBins else 'uniform',nbins[0],nbins[1])
    if equalStatsBins: # create equal statistics binning
        from array import array
        
        fixBinSumW = True
        if fixBinSumW: print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Binning target distribution. Bins per dimentsion will have fixed number of weights!'
        
        targetListX,  targetListY = [],[]
        for ev in target: targetListX += [ (ev.find(var[0]).getVal(), target.weight()) ] if fixBinSumW else [ (ev.find(var[0]).getVal(),1) ]
        for ev in target: targetListY += [ ev.find(var[1]).getVal() ]

# BELOW THIS IS TESTING REGION, DO NOT DELETE BEFORE OR AFTER THIS BOX
        targetListX.sort()
        targetListY.sort()


        from math import sqrt, cos, sin ,pi
        rhoDistr = []
        rho = lambda x,y: sqrt(x**2 + y**2)  
        theta = pi / 4
        x = lambda r: r * cos(theta)
        y = lambda r: r * sin(theta)
    
    #    for x_, y_ in zip(targetListX,targetListY): rhoDistr += [ rho(x_[0],y_) ] 
    #    rhoDistr.sort()

        eventCount = 0
        binstat = target.numEntries() / nbins if not fixBinSumW else target.sumEntries() / nbins 
        lowBounds = [ rho(xMin,yMin) ]
        while True:
            if eventCount + 1 >= target.numEntries(): break
            binContentCount = 0
            while binContentCount < binstat and eventCount < target.numEntries():
#                lowBound_i = rhoDistr[eventCount]
                lowBound_i = rho( targetListX[eventCount][0], targetListX[eventCount][0] )
                binContentCount += 1 if not fixBinSumW else targetListX[eventCount][1]  
                eventCount += 1
            lowBounds += [ lowBound_i ]

        lowboundsX = map( x, lowBounds )
        lowboundsY = map( y, lowBounds )
        
        lowboundsX, lowboundsY = array('f',lowboundsX), array('f',lowboundsY)

        print 'bincontent:',   binstat
# END OF TESTING BOX
        
        sourceHist = TH2D('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins[0], lowboundsX, nbins[1], lowboundsY )
        targetHist = TH2D('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins[0], lowboundsX, nbins[1], lowboundsY )

    else: # import binning
        sourceHist = TH2D('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins[0], xMin, xMax, nbins[1], yMin, yMax )
        targetHist = TH2D('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins[0], xMin, xMax, nbins[1], yMin, yMax )

    # create 2D histrograms and fill
    #  fill source 2D histogram
    if combineWeights and source.isWeighted(): # alow reweighting independant from the previous ones
        sourcePreviousWeights = []
        for evnt in source:
            sourceHist.Fill( _valX(evnt), _valY(evnt), source.weight() )
            sourcePreviousWeights += [ source.weight() ]
    else: 
        for evnt in source: sourceHist.Fill( _valX(evnt), _valY(evnt), source.weight() )

    #  fill target 2D histogram
    for evnt in target: targetHist.Fill( _valX(evnt), _valY(evnt), target.weight() )

    # rescale
    if source.numEntries() > target.numEntries(): sourceHist.Scale( target.sumEntries() / source.sumEntries() )
    else: targetHist.Scale( source.sumEntries() / target.sumEntries() )

    # calculate weights
    weights = []
    prblEvens, t_zero, s_zero, both_zero, t_negative, s_negative = 0, 0, 0, 0, 0, 0 # count how many events have a problematic weight
    for event in source:
        bin = sourceHist.FindFixBin( _valX(event),  _valY(event) ) # get bin with the given varible values
        # make sure the event is not problematic
        if not targetHist.GetBinContent(bin) <= 0 and not sourceHist.GetBinContent(bin) <= 0 :
            weights += [targetHist.GetBinContent(bin) / sourceHist.GetBinContent(bin)]
        else:
            if targetHist.GetBinContent(bin)==0 and sourceHist.GetBinContent(bin)==0: both_zero += 1
            if targetHist.GetBinContent(bin)<=0:
                if    targetHist.GetBinContent(bin)==0: t_zero += 1
                else: t_negative += 1
            elif sourceHist.GetBinContent(bin)<=0:
                if    sourceHist.GetBinContent(bin)==0: s_zero += 1
                else: s_negative += 1
            weights += [0.]
            prblEvens += 1

    print 'P2VV - INFO: TwoDimentionalVerticalReweighting of variables (%s): %s problematic events%s'\
                        %(var[0]+' '+var[1],prblEvens,' due to:' if prblEvens else '.' )
    summary =  '  Source binning                      : %s\n'%s_zero  
    summary += '  Target binning                      : %s\n'%t_zero
    summary += '  Both                                : %s\n'%both_zero
    summary += '  Source events with negative weights : %s\n'%s_negative
    summary += '  Target events with negative weights : %s'%t_negative    
    if prblEvens:  print summary

    # combine previous weights with the latest ones
    if combineWeights and source.isWeighted():
        print 'P2VV - INFO: TwoDimentionalVerticalReweighting: Combining previous source weights with the latest ones.'
        from numpy import array
        weights = array(weights) * array(sourcePreviousWeights)

    # plot and print the 2d histograms
    if plot: 
        canvSourc, canvTarg = [ TCanvas(n,n) for n in ['source','target'] ]
        for hist, canv in zip([sourceHist,targetHist], [canvSourc,canvTarg]): 
            hist.SetStats(False)
            canv.cd()
            hist.Draw('LEGO')
            canv.Print(canv.GetName() + '_%s.pdf'%iterIdx)

        testS0 = TH1D('%s_Source'%var[0],'test%s_Source'%var[0], nbins[0], xMin, xMax )
        testT0 = TH1D('%s_Target'%var[0],'test%s_Target'%var[0], nbins[0], xMin, xMax )
        testS1 = TH1D('%s_Source'%var[1],'test%s_Source'%var[1], nbins[1], yMin, yMax )
        testT1 = TH1D('%s_Target'%var[1],'test%s_Target'%var[1], nbins[1], yMin, yMax )
        source0, source1, = [], []
        for ev in source: 
            source0 += [ _valX(ev) ]
            source1 += [ _valY(ev) ]
        for var0, var1, weight in zip( source0, source1, weights ):
            testS0.Fill(var0,weight)
            testS1.Fill(var1,weight)
        for ev in target:
            testT0.Fill(_valX(ev),target.weight())
            testT1.Fill(_valY(ev),target.weight())
        testS0.Scale( target.sumEntries() / source.numEntries() )
        testS1.Scale( target.sumEntries() / source.numEntries() )

        can = TCanvas('test','test')
        can.Divide(2,2)
        can.cd(1)
        testT0.Draw()
        testS0.Draw('same err')
        can.cd(2)
        testT1.Draw()
        testS1.Draw('same err')

        from P2VV.Utilities.Plotting import _P2VVPlotStash
        _P2VVPlotStash += [testS0,testT0,testS1,testT1]
        
        # from ROOT import TFile
        # f = TFile.Open('histograms.root','recreate')
        # f.cd()
        # for hist in [testS0,testT0,testS1,testT1,can, sourceHist, targetHist]: hist.Write() 
        # f.Close()
        # assert False
        
        del source, target
        return weights, [can, testS0,testT0,testS1,testT1,sourceHist,targetHist]
    else: 
        del source, target
        return weights


# function that reweighits a single source distribution to match a given target using a histogram
def OneDimentionalVerticalReweighting(source, target, nbins, var, **kwargs):
    print 'P2VV - INFO: OneimentionalVerticalReweighting: Reweighting variable %s in sample %s.'%(var,source.GetName())
    iterIdx        = kwargs.pop('iterationNumber',     0 )
    plot           = kwargs.pop('xCheckPlots',     False )
    equalStatsBins = kwargs.pop('equalStatBins',   False )
    combineWeights = kwargs.pop('combWeights',     False )

    from ROOT import TH1D, TCanvas

    # dataset value getter 
    _valX = lambda ev: ev.find(var).getVal() 

    # get axis ranges
    from P2VV.RooFitWrappers import RooObject
    xMin, xMax = RooObject._rooobject(var).getMin(), RooObject._rooobject(var).getMax(var)

    # create histrograms
    print 'P2VV - INFO: OneDimentionalVerticalReweighting: Using %s binnning with #bins = %s.'%('equal statistics' if equalStatsBins else 'uniform',nbins)
    if equalStatsBins: # create equal statistics binning
        from array import array
        
        fixBinSumW = False
        if fixBinSumW: print 'P2VV - INFO: OneDimentionalVerticalReweighting: Binning target distribution. Bins will have fixed number of weights!'
        
        targetList = []
        for ev in target: targetList += [ ( ev.find(var).getVal(), target.weight() )]if fixBinSumW else [ (ev.find(var).getVal(),1) ]
        targetList.sort()

        binstat = target.sumEntries() / nbins if fixBinSumW else target.numEntries() / nbins
        lowbounds = [ xMin ]
        eventCount = 0
        while True:
            if eventCount + 1 >= target.numEntries(): break
            binContentCount = 0
            while binContentCount < binstat and eventCount < target.numEntries():
                lowbound_i = targetList[eventCount][0]        # variable value
                binContentCount += targetList[eventCount][1]  # event weight
                eventCount += 1
            lowbounds += [lowbound_i]

        sourceHist  = TH1D('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins, array('f',lowbounds) )
        targetHist  = TH1D('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins, array('f',lowbounds) )
    else:
        sourceHist  = TH1D('h_'+source.GetName(), 'h_'+source.GetTitle(), nbins, xMin, xMax )
        targetHist  = TH1D('h_'+target.GetName(), 'h_'+target.GetTitle(), nbins, xMin, xMax )

    # fill reweighting histograms
    if combineWeights and source.isWeighted(): # alow reweighting independant from the previous ones
        sourcePreviousWeights = []
        for evnt in source:
            sourceHist.Fill( _valX(evnt), source.weight() )
            sourcePreviousWeights += [ source.weight() ]
    else: 
        for evnt in source: sourceHist.Fill( _valX(evnt), source.weight() )
    # fill target 2D histogram
    for evnt in target: targetHist.Fill( _valX(evnt), target.weight() )
       
    # rescale
    if source.numEntries() > target.numEntries(): sourceHist.Scale( target.sumEntries() / source.sumEntries() )
    else: targetHist.Scale( source.sumEntries() / target.sumEntries() )

    # calculate weights
    weights = []
    prblEvens, t_zero, s_zero, both_zero, s_negative, t_negative = 0, 0, 0, 0, 0, 0 # count how many events have a problematic weight
    for event in source:
        bin = sourceHist.FindFixBin( _valX(event) ) # get the bin with given var value
        # make sure the event is not problematic
        if not targetHist.GetBinContent(bin) <= 0 and not sourceHist.GetBinContent(bin) <= 0 :
            weights += [targetHist.GetBinContent(bin) / sourceHist.GetBinContent(bin)] # calculate weight
        else:
            if targetHist.GetBinContent(bin)==0 and sourceHist.GetBinContent(bin)==0: both_zero += 1
            if targetHist.GetBinContent(bin)<=0:
                if    targetHist.GetBinContent(bin)==0: t_zero += 1
                else: t_negative += 1
            elif sourceHist.GetBinContent(bin)<=0:
                if    sourceHist.GetBinContent(bin)==0: s_zero += 1
                else: s_negative += 1
            weights += [0.]
            prblEvens += 1

    print 'P2VV - INFO: OneDimentionalVerticalReweighting of variable %s: %s problematic events%s'%(var,prblEvens,' due to:' if prblEvens else '.' )
    summary =  '  Source binning                      : %s\n'%s_zero  
    summary += '  Target binning                      : %s\n'%t_zero
    summary += '  Both                                : %s\n'%both_zero
    summary += '  Source events with negative weights : %s\n'%s_negative
    summary += '  Target events with negative weights : %s'%t_negative    
    if prblEvens: print summary

    # combine previous weights with the latest ones
    if combineWeights and source.isWeighted():
        print 'P2VV - INFO: OneDimentionalVerticalReweighting: Combining previous source weights with the latest ones.'
        from numpy import array
        weights = array(weights) * array(sourcePreviousWeights)

    # check the result of the reweighting 
    if plot: 
        test_s = TH1D('test_s','test_s', nbins, xMin,xMax)
        test_t = TH1D('test_t','test_t', nbins, xMin,xMax)
        sourceEvtList, targetEvtList = [],[]
        for ev in source: sourceEvtList+=[ _valX(ev) ]
        for ev in target: targetEvtList+=[ _valX(ev) ]
        for ev, w in zip(sourceEvtList,weights): test_s.Fill(ev,w)
        for evnt in target: test_t.Fill(_valX(evnt), target.weight())
        # plot and print the histograms
        testCanv = TCanvas('test','test')
        test_t.Draw()
        test_s.Scale(target.sumEntries() / source.numEntries())
        test_s.Draw('same err')
       
        from P2VV.Utilities.Plotting import _P2VVPlotStash
        _P2VVPlotStash += [test_s,test_t]
        del source, target
        return weights, testCanv, test_t
    else: 
        del source, target
        return weights
 
class WeightedDataSetsManager(dict):
    def __init__( self, **kwargs ):       
        self['initSource']        = kwargs.pop('source', '')
        self['dataSets']          = dict( initSource = self['initSource'] )
        self['permanentDataSets'] = dict()
        
        self['WeightsLists']         = {}
        self['permanetnWeigtsLists'] = {}
        self['combinedWeights']      = []   
        
        self['latestDataSetPointer'] = 'initSource'
        self['iterationNumber'] = 0
        self['saveIntermediateDatasets'] = False

        print 'P2VV - INFO: WeightedDataSetsManager: Initialsed for sample %s.'%self['initSource'].GetName()

    def appendWeights( self, weightsName, weightsList, combWithPrevious = False, scale = True ):
        # check if the weights is of numpy type array
        import numpy
        if not type(weightsList)==numpy.ndarray: weightsList = numpy.array(weightsList)

        # optionally cmbine with the previous weights list
        if combWithPrevious and len(self['WeightsLists'].keys()) != 0: 
            print 'P2VV - INFO: appendWeights: Combining weights "%s" with "%s" (dot product).'%(weightsName,self._wName)
            prev = self['WeightsLists'][ self['latestDataSetPointer'] ]
            weightsList = prev * weightsList
        self['WeightsLists'][weightsName] = weightsList

        # scale weights to preserve number of events 
        if scale:
            print 'P2VV - INFO: appendWeights: Scaling sources sum of weights to the number of entries.'
            n_events = self['initSource'].numEntries()
            sumW = sum( self['WeightsLists'][weightsName] )
            self['WeightsLists'][weightsName] =  ( n_events / sumW ) * self['WeightsLists'][weightsName] 

        # bookkeeping
        self._wName = ''
        for name in self['WeightsLists'].keys(): self._wName += name + '_'
        self._wName += str(self['iterationNumber'])
        self['latestDataSetPointer'] = weightsName
        self['combinedWeightsName']  = 'weight_' + self._wName 
        
        # write weights
        self['dataSets'][weightsName] = self.writeWeights(self['dataSets']['initSource'], \
                                        'weight_' + self._wName, self['initSource'].GetName() + '_' + self._wName )

        # delete unused datasets
        if not self['saveIntermediateDatasets']:
            removeKeys = self['dataSets'].keys()
            for k in ['initSource', self['latestDataSetPointer']]: removeKeys.remove(k)
            for key in removeKeys: 
                print 'P2VV - INFO: appendWeights: Deleting dataset named ' + self['dataSets'][key].GetName()
                del self['dataSets'][key]
        
    def writeWeights( self, dataset, weightsName, writeDatasetName ):
        print 'P2VV - INFO: writeWeights: Creating dataset with name %s and weight name %s:'%(writeDatasetName,weightsName)
        from ROOT import RooArgSet, RooRealVar, RooDataSet
        Weights = self['WeightsLists'][self['latestDataSetPointer']]
        weightsVar = RooRealVar( weightsName, weightsName, 1, .9*min(Weights), 1.1*max(Weights) )
        weightsArgSet  = RooArgSet( weightsVar )
        weightsDataSet = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        for weight in Weights:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        _dataset  = RooDataSet( writeDatasetName, writeDatasetName, dataset.get(), Import = dataset )
        _dataset.merge( weightsDataSet )
        _Wdataset = RooDataSet( writeDatasetName, writeDatasetName, _dataset.get(), Import = _dataset, WeightVar = (weightsName,True) )

        del weightsDataSet, dataset, _dataset, Weights
        return _Wdataset

    def getDataSet( self, which='' ): 
        dataSetKey = which if which else self['latestDataSetPointer']
        if which in self['permanentDataSets'].keys(): return self['permanentDataSets'][which]
        else: return self['dataSets'][dataSetKey]

    def setDataSet ( self, data, weightsName ):
        # set weights name
        self['combinedWeightsName'] = data.GetName()[:-1] + weightsName + data.GetName()[-2:]
        
        # delete previous data and point to the new one
        if self['saveIntermediateDatasets']: self['dataSets'][weightsName] = data
        else:
            for key in self['dataSets'].keys(): 
                if key == self['latestDataSetPointer']: del self['dataSets'][key]
            self['dataSets'][weightsName] = data
        
        #set dataset pointer
        self['latestDataSetPointer'] = weightsName

    def plotWeights( self, which = '', Range=() ):
        from ROOT import TCanvas, TH1F, gStyle
        from P2VV.Utilities.Plotting import _P2VVPlotStash

        which = which if which else self['latestDataSetPointer']
        c = TCanvas( 'weights_' + which + str(self['iterationNumber']), 'weights_' + which + str( self['iterationNumber']) )        
    
        # plot range
        weightsList = self['WeightsLists'][which]        
        plotRange = ( -3, 3 ) if not Range else Range    
                    
        # create histogram fill and draw
        weghtsHist = TH1F('weights_'+which, 'weights_'+which, 200, plotRange[0], plotRange[1])
        for weight in weightsList: weghtsHist.Fill(weight)
        c.cd()
        gStyle.SetOptStat("nmruo")
        weghtsHist.SetMarkerSize(.5)
        weghtsHist.Draw('err')
        c.Print( 'weights_%s_.pdf'%which )
        _P2VVPlotStash +=[c,weghtsHist]

    def getWeightName( self ): return self._wName    
        
    def clear( self ):
        del self['dataSets'], self['WeightsLists'],  self['combinedWeights']
        
        # restore initial and permanent datasets 
        self['dataSets'] =  dict( initSource = self['initSource'] )
                
        # restore weights container and flag
        self['WeightsLists'] = dict()
        self['latestDataSetPointer'] = 'initSource'


# Vertical reweighting class to match physics of weighted distribution using a pdf
class MatchPhysics( ):
    def __init__( self, nTupleFile, nTupleName, **kwargs ):      
        # monte carlo gen conditions specifier
        MCProd = kwargs.pop('MonteCarloProduction', '2011')        
        print 'P2VV - INFO: Initialised physics reweighting class: MatchPhysics().'        
    
        # set global object name prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        self._namePF = 'mc'
        setParNamePrefix( self._namePF )
        
        # use pdf config to grab default configuration
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_RunIAnalysis as PdfConfig
        pdfConfig = PdfConfig( RunPeriods = '3fb' )
        
        # blind parameters in MC pdf
        blind = kwargs.pop('blind', False)
        blindStr = pdfConfig['blind'] if blind else False

        # transversity amplitudes
        A0Mag2Val    = 0.722**2 / (0.722**2 + 0.480**2 + 0.499**2)
        AperpMag2Val = 0.499**2 / (0.722**2 + 0.480**2 + 0.499**2) 
        AparMag2Val  = 0.480**2 / (0.722**2 + 0.480**2 + 0.499**2)

        A0PhVal    = 0.
        AperpPhVal = 3.07
        AparPhVal  = 3.30

        # build pdf with KK mass bins and category states.
        KKMassStates = dict( [ ('bin%d' % i, i) for i in range(6) ] )
        SWaveAmps    = dict( mc_f_S        = dict()
                            ,mc_ASOddPhase = dict()
                            ,mc_C_SP       = dict()
                             )
        for k in SWaveAmps.keys():
            for bin in xrange(6): SWaveAmps[k]['bin%s'%bin] = 0

        # CP violation parameters
        phiCPVal  = +0.07 
        lambCPVal = 1.

        # B lifetime parameters
        GammaVal  = 1. / 1.503 
        dGammaVal = 1. / 1.406 - 1. / 1.614
        dMVal     = 17.8
        tResSigma = 0.045

        angleNames = ( 'cos(#theta_{K})', 'cos(#theta_{#mu})', '#phi_{h}' )

        # load roofit wrappers
        from P2VV.Load import RooFitOutput
        
        # get workspace
        from P2VV.RooFitWrappers import RooObject
        ws = RooObject().ws()
      
        # angular functions
        if ws['helphi']: # there is a funny think going on with the range of helphi not sth to worry.
            from math import pi
            ws['helphi'].setRange(-pi,pi)
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )
        self._angleFuncs = angleFuncs

        # get observables (get the wrapper instead of the base objects!!)
        trueTime  = _createGetObservable('truetime')
        time      = RooObject._rooobject('time') if ws['time'] else _createGetObservable('time')
        angles    = [ RooObject._rooobject(o) for o in ['helcosthetaK','helcosthetaL','helphi'] ]
        iTag      = _createGetObservable('iTag')
        KKMass    = RooObject._rooobject('KKMass') if ws['KKMass'] else _createGetObservable('KKMass')
        KKMassCat = RooObject._rooobject('KKMassCat') if ws['KKMassCat'] else _createGetObservable('KKMassCat')
        self._obsSet = [ trueTime, time, KKMass ] + angles + [ KKMassCat ]
        self._normSet = angles

        # set momenta range and put them in obsSet
        print 'P2VV - INFO: Setting track and B momenta ranges.'
        from P2VV.Utilities.MCReweighting import trackMomentaRanges
        from P2VV.RooFitWrappers import RealVar
        for var in [ '%s_%s'%( part, comp ) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'PX', 'PY', 'PZ', 'P' ) ] + ['B_P','B_Pt']:
            try: var = RooObject._rooobject(var)
            except KeyError: var = RealVar( var, Title = var, Unit = 'MeV/c^2', MinMax = ( trackMomentaRanges[var][0], trackMomentaRanges[var][1]) )
            self._obsSet += [ var ]
        
        # read ntuple
        from P2VV.Utilities.DataHandling import readData
        readOpts = { 'ntupleCuts' : 'mass>5350 && mass<5355' } if kwargs.pop('Reduced', False) else  { }
        self._data = readData( nTupleFile, dataSetName=nTupleName, NTuple=True, observables=self._obsSet, **readOpts)
        self._data.SetName( 'mcData_' + MCProd )
           
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

        for p in self._pdf.Parameters():
            if 'C_SP' in p.GetName():
                try: binIndex = int(p.GetName().partition('bin')[2])
                except ValueError: assert False, 'P2VV - ERROR: MatchPhysics: Cannot set C_SP factors'
                p.setVal( pdfConfig['CSPValues'][binIndex] )
                p.setConstant()
    
    def setMonteCarloParameters(self, pars=None, mcParPrefix=''):
        print 'P2VV - INFO: setMonteCarloParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        from P2VV.Utilities.MCReweighting import parValuesMcSim08_6KKmassBins as mcPars
                
        for par in self._pdf.Parameters():
           if par.GetName() == 'dummyBlindState' or par.isConstant(): continue
           key = par.GetName().replace('mc_', mcParPrefix)
           if mcPars.has_key(key): 
               par.setVal( mcPars[key] )
               print '%20s %.4f'%(par.GetName(), par.getVal())
           else: print 'P2VV - ERROR: setMonteCarloParameters: Cannot find parameter %s in %s dictionary.'(par.GetName(),mcPars['name'])
                            
    def setDataFitParameters(self, dataPars, dataParPrefix=''):
        print 'P2VV - INFO: setDataFitParameters: Setting the following parameters to the monte carlo pdf, named %s.'%self._pdf.GetName()
        from P2VV.Parameterizations.FullPDFs import PdfConfiguration
        pdfConfig = PdfConfiguration() 
        pdfConfig.readParametersFromFile( filePath = dataPars )
                
        # print currect pdf parameter values
        for par in self._pdf.Parameters():
            if par.GetName() == 'dummyBlindState' or par.isConstant(): continue
            key = par.GetName().replace('mc_', dataParPrefix)
            if pdfConfig.parameters().has_key(key):
                par.setVal( pdfConfig.parameters()[key][0] )
                if 'phiCP' in key or 'dGamma' in key: print '%20s %s'%(par.GetName(), '(blinded)' )
                else: print '%20s %.4f'%(par.GetName(), par.getVal())
            else: assert False, 'P2VV - ERROR:setDataFitParameters: Cannot find parameter %s in data physics parameters file %s'%(key,dataPars)
        
    def calculateWeights(self, iterNumb, dataParameters):
        print 'P2VV - INFO: Matching physics on mc sample.'
        self._iterNumb = iterNumb

        from ROOT import RooArgSet
        normVars =  RooArgSet( self._normSet )
        
        # Reweights MC verticaly to match the Physics of data.
        nominators, denominators = [], []
        
        def calculatePhysicsWeights(nom=[],den=[]):
            self._physWeights = []
            count = 0
            for idx in xrange(len(nom)): # if pdf.getVal() returns 0 for a given event make it unweighted
                if nom[idx] == 0 or den[idx] == 0: 
                    nom[idx], den[idx] = 0., 0.
                    count += 1 # count how many evetns have a problematic weight
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
    def setDataSet(self, data):     self._data = data 
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
            print 'P2VV - INFO: _mimicWeights: Mimicing weighted distribution of source variable:', var
            self._mimicedVars['inDistr'][var] = self._MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins, 'in')
            
            if self._mimicedVars.has_key('outDistr'):
                if not self._mimicedVars['outDistr'].has_key(var):
                    print 'P2VV - INFO: _mimicWeights: Mimicing weighted distribution of target variable:', var
                    self._mimicedVars['outDistr'][var] = self._MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out') 
            else:
                print 'P2VV - INFO: _mimicWeights: Mimicing weighted distribution of target variable:', var
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

    def reweight( self, itNum, data ):
        # get rid of the previous tree if any
        self._inTree = data
        self._itNum = itNum
        
        # get physics weight variable
        if not self._inTree.get().find( self._inWeightName ):
            from ROOT import RooRealVar, RooNumber
            RooInf = RooNumber.infinity()
            self._physWeightsVar = RooRealVar( self._inWeightName, self._inWeightName, -RooInf, RooInf )
        else: self._physWeightsVar = self._inTree.get().find( self._inWeightName )
        
        # convert roodataset to tree
        from ROOT import gROOT
        gROOT.cd('PyROOT:/')
        self._inTree = self._inTree.buildTree(WeightName=self._inWeightName)
        self._outTree = self._outTree.buildTree(WeightName=self._outWeightName)

        # mimic the weights, transform Kaon momenta and recalculate angles 
        self._mimicWeights()
        self._TransformKaonMomentaAndRecalculateAngles()
        
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
        name = self._inTree.GetName()[:-1] + 'horKKmom' + self._inTree.GetName()[-2:]
        self._recalculatedData = RooDataSet( name, name, recalculatedVars  )        
        copiedData             = RooDataSet('copiedData', 'copiedData', self._inTree, copiedVars )

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
        print 'P2VV - INFO: _ReweightAndTransformAngles: Transforming p(K+),p(K-) and recalculating decay angles.'
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

            self._recalculatedData.add( recalculatedVars )
        
        self._recalculatedData.merge( copiedData )
        self._recalculatedData = RooDataSet(self._recalculatedData.GetName(), self._recalculatedData.GetTitle(), 
                                            RooArgSet(recalculatedVars, copiedVars), 
                                            Import = self._recalculatedData, 
                                            WeightVar = (self._physWeightsVar, True)
                                            )
        del self._inTree, self._outTree
  
    def getDataSet( self ): return self._recalculatedData


## micelanous impots / helping stuff

# MC generating conditions  
parValuesMcSim08_6KKmassBins = dict(
     name              = 'Sim08_Conditions'
    ,A0Mag2            = 0.722**2 / (0.722**2 + 0.480**2 + 0.499**2)
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
    ,lambdaCP          = 1.    
    ,dM                = 17.8
    ,f_S_bin0          = 0.
    ,f_S_bin1          = 0.
    ,f_S_bin2          = 0.
    ,f_S_bin3          = 0.
    ,f_S_bin4          = 0.
    ,f_S_bin5          = 0.
    )

# trackand B momenta ranges 
trackMomentaRanges = dict(
    Kminus_PX  = [ -5e5 , 5e5 ],
    Kminus_PY  = [ -5e5 , 5e5 ],
    Kminus_PZ  = [    0 , 5e5 ],
    Kminus_P   = [    0 , 5e5 ],
    Kplus_PX   = [ -5e5 , 5e5 ],
    Kplus_PY   = [ -5e5 , 5e5 ],
    Kplus_PZ   = [    0 , 5e5 ],
    Kplus_P    = [    0 , 5e5 ],
    muminus_PX = [ -9e5 , 6e5 ],
    muminus_PY = [ -9e5 , 6e5 ],
    muminus_PZ = [    0 , 1e6 ],
    muminus_P  = [    0 , 1e6 ],
    muplus_PX  = [ -5e5 , 6e5 ],
    muplus_PY  = [ -5e5 , 6e5 ],
    muplus_PZ  = [ -5e5 , 1e6 ],
    muplus_P   = [    0 , 1e6 ],
    B_P        = [    0 , 1e6 ],
    B_Pt       = [    0 , 2e5 ]    
    )

# bokkeeping dictionaries for plotting
# plotingScenarios = dict( BmommkkphysKKmom = [ ('mcData','BmomRewData'), ('BmomRewData','mkkRewData'), ('mkkRewData','mcDataPhysRew'), ('mcDataPhysRew','MomRewData') ],
#                          BmomphysKKmom    = [ ('mcData','BmomRewData'), ('BmomRewData','mcDataPhysRew'), ('mcDataPhysRew','MomRewData') ],
#                          mkkphysKKmom     = [ ('mcData','mkkRewData'), ('mkkRewData','mcDataPhysRew'), ('mcDataPhysRew','MomRewData') ],
#                          physKKmom        = [ ('mcData','mcDataPhysRew' ), ('mcDataPhysRew','MomRewData') ]
#                          )

plotingScenarios = dict( BmommkkphysKKmom     = [ ('mcData','BmomRewData'), 
                                                ('mcData','BmomRewData','mkkRewData'), 
                                                ('mcData','BmomRewData','mkkRewData','mcDataPhysRew'), 
                                                ('mcData','BmomRewData','mkkRewData','mcDataPhysRew','MomRewData')
                                                  ],
                         TwoDBmommkkphysKKmom = [ ('mcData','BmomMkkRewData'), 
                                                  ('mcData','BmomMkkRewData','mcDataPhysRew'), 
                                                  ('mcData','BmomMkkRewData','mcDataPhysRew','MomRewData')
                                                  ],
                         mkkBmomphysKKmom     = [ ('mcData','mkkRewData'), 
                                                ('mcData','mkkRewData','BmomRewData'),                
                                                ('mcData','mkkRewData','BmomRewData','mcDataPhysRew') , 
                                                ('mcData','mkkRewData','BmomRewData','mcDataPhysRew','MomRewData') 
                                                  ],
                         BmomphysKKmom        = [ ('mcData','BmomRewData'), 
                                                ('mcData','BmomRewData','mcDataPhysRew'), 
                                                ('mcData','BmomRewData','mcDataPhysRew','MomRewData') 
                                                  ],
                         mkkphysKKmom         = [ ('mcData','mkkRewData'), 
                                                ('mcData','mkkRewData','mcDataPhysRew'), 
                                                ('mcData','mkkRewData','mcDataPhysRew','MomRewData')
                                                  ],
                         physKKmom            = [ ('mcData','mcDataPhysRew'), 
                                                ('mcData','mcDataPhysRew','MomRewData') 
                                                  ]
                         )

weightNamesDataKeysMap = dict( mcData         = 'initSource', 
                               mcDataPhysRew  = 'phys',
                               mkkRewData     = 'mKK',
                               MomRewData     = 'KKmom',
                               BmomRewData    = 'Bmom',
                               BmomMkkRewData = 'Bmom_mKK'
                               )
