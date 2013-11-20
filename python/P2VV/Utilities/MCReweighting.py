###########################################################################################################################################
## Utilities.MCReweighting: P2VV utilities for reweighting Monte Carlo data to create desired distributions                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   VS, Vasilis Syropoulos, Nikhef, v.syropoulos@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################

 # keep reference to canvases, legends, frames
from P2VV.Utilities.Plotting import _P2VVPlotStash

def compareDistributions( **kwargs ):
    # TODO: Provide observables set as an argument, do not get the pbsset from the dataset
    mcData          = kwargs.pop('mcData',          None)
    momRewData      = kwargs.pop('momRewData',      None)
    sData           = kwargs.pop('sData',           None)
    obsSet          = kwargs.pop('obsSet',          None)
    itNumb          = kwargs.pop('itNumb',          None)
    physWeightsName = kwargs.pop('physWeightsName', None)
    nullTest        = kwargs.pop('nullTest',        None)

    # get mcData before and after physics reweighting
    from ROOT import RooFit, RooDataSet
    mcBeforePhysRew = mcData
    mcAfterPhysRew  = RooDataSet( 'MC_BeforePhysRew_%s'%itNumb, 'MC_BeforePhysRew_%s'%itNumb, 
                                  mcBeforePhysRew.get(), 
                                  Import    = mcData,
                                  WeightVar = (physWeightsName, True)
                                  )
   
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

    # assymetry plots are compared w.r.t. the sData
    referenceHistName = 'h_' + sData.GetName() 

    # start drawing
    from ROOT import TCanvas, RooAbsData, TPaveText, kGreen, kMagenta
    from P2VV.Utilities.Plotting import compareDataSets, makeAssymetryPlot
    from P2VV.Load import LHCbStyle
    from math import pi

    obsCanv = TCanvas('observables_%s'%itNumb,'observables_%s'%itNumb)
    obsCanv.Divide(2,2)
    KaonCanv = TCanvas('KaonMomenta_%s'%itNumb,'KaonMomenta_%s'%itNumb)
    KaonCanv.Divide(4,2)
    muonCanv = TCanvas('muonMomenta_%s'%itNumb,'muonMomenta_%s'%itNumb)
    muonCanv.Divide(4,2)
    assymKaonCanv = TCanvas('assymKaonMomenta_%s'%itNumb,'assymKaonMomenta_%s'%itNumb)
    assymKaonCanv.Divide(4,2)
    assymMuonCanv = TCanvas('assymmuonMomenta_%s'%itNumb,'assymmuonMomenta_%s'%itNumb)
    assymMuonCanv.Divide(4,2)

    # set some data drawing options
    stdDrawOpts = dict( DataError=RooAbsData.SumW2, MarkerSize = .6, XErrorSize = 0 )
    colors      = dict( mcBefore=2, mcAfter=kGreen+3, MomRewData=4, Sdata=kMagenta+2)
 
    # plot angles and decay time
    print 'P2VV - INFO: compareDistributions: Plotting decay angles and time.'
    for canv, obs, logY, rangeX, rangeY in zip( [obsCanv.cd(i+1) for i in range(len(observables))], 
                                                observables,
                                                3 * [False] + [True],
                                                [ (-1.,1.), (-1.,1.), (-pi,pi), (0.,14.) ] if not nullTest else \
                                                [ (-pi,pi), (-1.,1.), (-1.,1),   (0.,14.) ],
                                                [ [400.,1650.], [500.,1200.], [600.,1300.], [.1,5e4] ] if not nullTest else\
                                                [ [1e3,2*5e3], [2e3,2*5e3], [15e2,2*65e2], [.1,2*5e5] ]
                                         # else [ [0,26500], [0,19000], [0,17000],   [.1,5e5] ]
                                                ): 
        compareDataSets( canv, obs, 
                         data      = dict( mcBefore   = mcBeforePhysRew, 
                                           mcAfter    = mcAfterPhysRew, 
                                           MomRewData = momRewData,      
                                           Sdata      = sData                                                        
                                           ),
                         dataOpts  = dict( mcBefore   = dict( MarkerColor = colors['mcBefore'],   **stdDrawOpts ),
                                           mcAfter    = dict( MarkerColor = colors['mcAfter'],    **stdDrawOpts ),
                                           MomRewData = dict( MarkerColor = colors['MomRewData'], **stdDrawOpts ),
                                           Sdata      = dict( MarkerColor = colors['Sdata'],      **stdDrawOpts ) 
                                           ),
                         frameOpts = dict( Bins = 30, Range=rangeX ),
                         logy      = logY,
                         RangeY    = rangeY
                         )
        
    # plot Kaon and muon momenta
    print 'P2VV - INFO: compareDistributions: Plotting track momenta.'
    for canv, assymCanv, obs, rangeY in zip( 
        [ KaonCanv.cd(k+1) for k in range(len(Kmomenta)) ]      + [ muonCanv.cd(m+1) for m in range(len(muMomenta)) ],
        [ assymKaonCanv.cd(k+1) for k in range(len(Kmomenta)) ] + [ assymMuonCanv.cd(m+1) for m in range(len(muMomenta)) ],
        Kmomenta + muMomenta,
        2*[ (0.,3000), (0.,2500), (0.,2500), (0.,3000)] + 2*[ (0.,4000), (0.,3000), (0.,3000), (0.,4000)] if not nullTest else \
        2*[ (0.,3000), (0.,4500), (0.,4500), (0.,3000)] + 2*[ (0.,2000), (0.,3000), (0.,3000), (0.,2000) ]
        ): 
        momFrame = compareDataSets( canv, 
                                    obs, 
                                    data      = dict( mcBefore   = mcBeforePhysRew, 
                                                      mcAfter    = mcAfterPhysRew, 
                                                      MomRewData = momRewData,      
                                                      Sdata      = sData
                                                      ),
                                    dataOpts  = dict( mcBefore   = dict( MarkerColor = colors['mcBefore'],    **stdDrawOpts ),
                                                      mcAfter    = dict( MarkerColor = colors['mcAfter'],     **stdDrawOpts ),
                                                      MomRewData = dict( MarkerColor = colors['MomRewData'],  **stdDrawOpts ),
                                                      Sdata      = dict( MarkerColor = colors['Sdata'],       **stdDrawOpts )
                                                      ),
                                    frameOpts = dict( Bins = 50, Range=trackMomRangeX[obs.GetName()] ),
                                    RangeY    = rangeY if not nullTest else None
                                    )
    
        makeAssymetryPlot( assymCanv, momFrame, referenceHistName ) 
 
    # make a legend and draw it
    legend = TPaveText( .47, .66, .77, .9, 'NDC' )
    legend.SetFillColor(0)
    legend.AddText('#color[%s]{%s}'%( colors['mcBefore'],   'SourceBeforePhysRew' ) )
    legend.AddText('#color[%s]{%s}'%( colors['mcAfter'],    'SourceAfterPhysRew'  ) )
    legend.AddText('#color[%s]{%s}'%( colors['MomRewData'], 'SoourceAfterMomRew'  ) )
    legend.AddText('#color[%s]{%s}'%( colors['Sdata'],           'Target'         ) )

    for canv, pad in zip( [obsCanv, KaonCanv, muonCanv, assymKaonCanv, assymMuonCanv], [4,8,8,8,8] ): 
        canv.cd(pad)
        legend.Draw()
    _P2VVPlotStash.append(legend) # keep reference to the legend
    
    mcAfterPhysRew.IsA().Destructor(mcAfterPhysRew)
    return [obsCanv,KaonCanv,muonCanv,assymKaonCanv,assymMuonCanv,momFrame]



class UniFunc(object):
    """ Diego's Uniform Function: A function that transform a variable into a flat distribution
    """

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


# combine datasets
def _combineDataSetParts( files, name, weightName='' ):
    print 'P2VV - INFO: combineDataSetParts: Combining the following datasets with common name, %s'%name
    for f in files: print f
 
    from ROOT import TFile
    dataSets = [ TFile.Open(f,'READ').Get(name) for f in files ]
    data = dataSets.pop()
    for d in dataSets: data.append(d)

    # import args into current workspace
    from P2VV.RooFitWrappers import RooObject
    ws = RooObject().ws()
    for arg in data.get():
        if not ws[arg.GetName()]: ws.put(arg) 

    # destroy and delete unnecessary stuff
    for d in dataSets: d.IsA().Destructor(d)
    del dataSets
    return data


# easily create an observable using information from the PdfConfig class
def _createObservable(name):
    from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis
    PdfConfig = Bs2Jpsiphi_2011Analysis()
    obsDict = PdfConfig['obsDict']
    
    if name == 'helcosthetaK': name = 'cpsi'
    if name == 'helcosthetaL': name = 'ctheta'
    if name == 'helphi'      : name = 'phi'
    if name == 'iTag'        : obsDict['iTag']    = ( 'iTag', 'Initial state flavour tag', { 'Untagged' : 0 } ) 
    #if name == 'truetime'    : obsDict['truetime'] = ( 'iTag', 'Initial state flavour tag', { 'Untagged' : 0 } ) 

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
    #return obs 



# Vertical reweighting class of MC to match the physics of sWeighted data.
class MatchPhysics( ):
    def __init__( self, nTupleFile, **kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
        self._pdfConfig = PdfConfig()

        # efficiencies files
        self._pdfConfig['timeEffType']      = kwargs.pop('timeEffType', '')
        self._pdfConfig['timeEffHistFiles'] = kwargs.pop('timeEffHistFiles', {} )
        self._pdfConfig['anglesEffType']    = kwargs.pop('anglesEffType', '')
        self._pdfConfig['angEffMomsFiles']  = kwargs.pop('angEffMomsFiles', '')
        
        # decay time resolution configuration
        self._pdfConfig['timeResType']           = 'eventNoMean'
        self._pdfConfig['numTimeResBins']        = 40
        
        # cpv paramaetrization
        self._pdfConfig['lambdaCPParam'] = 'lambPhi'

        # KK mass parameters
        self._pdfConfig['paramKKMass']       = 'simultaneous'
        self._pdfConfig['KKMassBinBounds']   = [ 990., 1050. ]
        self._pdfConfig['CSPValues']         = [ 0.498]

        # tagging configuration
        self._pdfConfig['constrainTagging'] = ''
       
        # external constrints
        from P2VV.Imports import extConstraintValues
        self._pdfConfig['externalConstr']   = {}
        extConstraintValues.setVal( 'DM',      ( 17.63, 0.11 ) )
        extConstraintValues.setVal( 'P0OS',    (  0.392, 0.008, 0.392 ) )
        extConstraintValues.setVal( 'DelP0OS', (  0.0110, 0.0034 ) )
        extConstraintValues.setVal( 'P1OS',    (  1.000,  0.023  ) )
        extConstraintValues.setVal( 'DelP1OS', (  0.000,  0.001  ) )
        extConstraintValues.setVal( 'P0SS',    (  0.350, 0.017, 0.350 ) )
        extConstraintValues.setVal( 'DelP0SS', ( -0.019, 0.005   ) )
        extConstraintValues.setVal( 'P1SS',    (  1.00,  0.16    ) )
        extConstraintValues.setVal( 'DelP1SS', (  0.00,  0.01    ) )

        # read data
        if type(nTupleFile) == list: 
            from P2VV.Utilities.MCReweighting import _combineDataSetParts
            self._pdfConfig['signalData'] = _combineDataSetParts(nTupleFile, kwargs.pop('nTupleName', 'JpsiKK') )
        else:
            from P2VV.Utilities.DataHandling import readData
            self._pdfConfig['signalData'] = readData( filePath    = nTupleFile, 
                                                      dataSetName = kwargs.pop('nTupleName', 'JpsiKK'),
                                                      NTuple      = False 
                                                        )
        self._pdfConfig['runPeriods'] = []
        self._pdfConfig['readFromWS'] = True

        from math import pi
        from P2VV.RooFitWrappers import RooObject
        ws = RooObject().ws()['helphi'].setRange(-pi,pi)

        # monte carlo paramters 
        self._mcPars     = kwargs.pop('mcParameters', {} ) 
                
        # give a prefix to all objects created by the pdf builder
        self._pdfConfig['parNamePrefix'] = kwargs.pop('parNamePrefix','mc') 
        
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
        self._pdfBuilder = PdfBuilder( **self._pdfConfig )
        self._pdf = self._pdfBuilder.pdf()

        # physics parameters names
        self._physicsParameters = {}
        for par in self._pdf.Parameters():
            if par.GetName().partition('_')[2].startswith('A') and \
                    ('Mag2' in par.GetName().partition('_')[2] or 'Phase' in par.GetName().partition('_')[2] ): 
                self._physicsParameters[par.GetName()] = par
            elif 'f_S' in par.GetName() and 'bin' in par.GetName(): 
                self._physicsParameters[par.GetName()] = par
            elif par.GetName().partition('_')[2] in ['dM', 'dGamma', 'Gamma', 'phiCP', 'lambdaCP']: 
                self._physicsParameters[par.GetName()] = par

        # import mc physics parameters
        self._AngAmpsParsVals = kwargs.pop('monteCarloParams')
        

    def setMonteCarloParameters(self, pars=None):
        print 'P2VV - INFO: setMonteCarloParameters: Setting the following parameters as mc parameters.'
        for key in self._physicsParameters.keys():
            self._physicsParameters[key].setVal( self._AngAmpsParsVals[ key.partition('_')[2]] ) 
        for k in self._physicsParameters.keys(): 
            print '%20s %.4f'%(self._physicsParameters[k].GetName(), self._physicsParameters[k].getVal())
                    
    def setDataFitParameters(self, dataPars, KKmassCat=None):
        print 'P2VV - INFO: setDataFitParameters: Setting the following parameters as data parameters.'
        for key in self._physicsParameters.keys():
            self._physicsParameters[key].setVal( dataPars[ key.partition('_')[2]] ) 
        for k in self._physicsParameters.keys(): 
            print '%20s %.4f'%(self._physicsParameters[k].GetName(), self._physicsParameters[k].getVal())
        else: print ' ' # TODO: LoopOver the KKbins and set f_S_i and ASOddPhase_i


    def calculateWeights(self, iterNumb, dataParameters):
        self._iterNumb = iterNumb

        from ROOT import RooArgSet
        normVars =  RooArgSet( self._pdfBuilder['obsSetP2VV'] )
        
        # Reweights MC verticaly to match the Physics of data.
        nominators, denominators = [], []
        
        def _calculatePhysicsWeights(nom=[],den=[]):
            self._physWeights = []
            count = 0
            for idx in xrange(len(nom)): # if pdf.getVal() returns 0 for a given event make it unweighted
                if nom[idx] == 0 or den[idx] == 0: 
                    nom[idx], den[idx] = 1., 1.
                    count += 1 # count how many evetns have a problematic weight assignment 
                self._physWeights += [ nom[idx]/den[idx] ]
            print 'P2VV - WARNING: calculateWeights: For %s events out of %s pdf value is zero.'%(count,len(nom))
            return self._physWeights

        # loop over events and calculate physics weights
        self._pdf.attachDataSet( self._pdfConfig['signalData'] ) # make the pdf directly dependant on data
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'    
        self.setMonteCarloParameters()
        for nev in xrange(self._pdfConfig['signalData'].numEntries()):
            self._pdfConfig['signalData'].get(nev)
            denominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating nominators for phyisics matching weights'
        self.setDataFitParameters(dataParameters)
        for nev in xrange(self._pdfConfig['signalData'].numEntries()):
            self._pdfConfig['signalData'].get(nev)
            nominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating phyisics matching weights'
        _calculatePhysicsWeights(nom=nominators, den=denominators)


    def writeWeights(self, weightsName='weightPhys'):
        from ROOT import RooArgSet, RooRealVar, RooDataSet, RooNumber
        from P2VV.RooFitWrappers import RealVar
        RooInf  = RooNumber.infinity()

        self._weightsName = weightsName
        weightsVar        = RooRealVar( self._weightsName, self._weightsName, -RooInf, RooInf   )
        weightsArgSet     = RooArgSet( weightsVar )
        weightsDataSet    = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        
        for weight in self._physWeights:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        
        self._pdfConfig['signalData'].merge( weightsDataSet )
        self._pdfConfig['signalData'].SetName('MC_physicsReweighted_%s_iteration'%self._iterNumb )
        print 'P2VV - INFO: Phyisics matching weights added to dataset: '+'MC_physicsReweighted_%s_iteration'%self._iterNumb

        del self._physWeights, weightsVar, weightsArgSet, 
        weightsDataSet.IsA().Destructor(weightsDataSet)
        
    

    def getPdf(self):               return self._pdf
    def getAngleFunctions(self):    return self._pdfBuilder['angleFuncs']
    def getDataSet(self):           return self._pdfConfig['signalData']
    def getMcTime(self)  :          return  [ o for o in self._pdf.Observables() if 'time' in o.GetName() ]
    def getProjDataset(self):
        return getDataset.reduce(getPdf.ConditionalObservables() + self.pdf.indexCat())
    def getNtupleVars(self):
        varSet = []
        for name in [ '%s_%s' % ( part, comp ) for part in ['Kplus','Kminus','muplus','muminus'] for comp in ('P','PX','PY','PZ') ] + ['B_P', 'B_Pt']:
            varSet.append( self._pdf.ws()[name] )
        return varSet


# Match MC to sWeighted data with horizontal reweighting of B_P and recalculate angles.
class MatchWeightedDistributions():
    def __init__( self,  **kwargs ):
        print 'P2VV - INFO: Initialised kinematic reweighting class: matchWeightedDistributions().'
        self._inTree          = kwargs.pop('inTree', '')  # distribution to be modified.
        self._outTree         = kwargs.pop('outTree')     # distribution to be matched with.
        self._inWeightName    = kwargs.pop('inWeightName')
        self._outWeightName   = kwargs.pop('outWeightName')
        self._mcObsSet        = kwargs.pop('observables')
        self._copyVars        = kwargs.pop('spectatorVars', None)
               
        self._nBins = kwargs.pop('nBins',    '1000'     )
        self._vars  = kwargs.pop('whichVars', 'Kminus_P')
        self._itNum = kwargs.pop('itNum',         0     )
       
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
                if var.GetName().startswith('mu'):  self._muonSet.append(var)
                elif var.GetName().startswith('K'): self._KaonSet.append(var)
                elif var.GetName().startswith('B'): self._BmomSet.append(var)
                else: continue
            

    def _mimicWeights(self):
        # TODO: Mimic output distribution only once
        self._mimicedVars = dict( inDistr={} )
        for var in self._vars:
            print 'P2VV - INFO: Input Distribution: Mimicing weights of variable: ', var
            self._mimicedVars['inDistr'][var]  = self._MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins, 'in')

            if self._mimicedVars.has_key('outDistr'):
                print 'Entrer first if statement'
                print self._mimicedVars.kyes()
                if not self._mimicedVars['outDistr'].has_key(var):
                    print 'Enter second if statement'
                    print self._mimicedVars.keys()
                    print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
                    self._mimicedVars['outDistr'][var] = self._MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out') 

            else:
                print 'Enter else statement'
                print self._mimicedVars.keys()
                print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
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
        self._physWeightsVar = self._inTree.get().find( self._inWeightName ) # get physics weights var         
       
        # get dataset scales
        self._inSumW = self._inTree.sumEntries()
        self._dataSetsScale = self._outSumW / self._inSumW
        self._physRewScale  = self._inSumW / self._inTree.numEntries()
        assert self._physRewScale <1.01, \
            'P2VV - WARNING: Scale between original MC and physics reweighted MC is %s, remove weights and then perform momentum reweighting.' %self._physRewScale 
        
        # convert roodataset to tree
        from ROOT import RooDataSet
        if type(self._inTree)==RooDataSet: 
            from ROOT import gROOT
            gROOT.cd('PyROOT:/')
            self._inTree = self._inTree.buildTree(WeightName=self._inWeightName) 

        # mimic the weights 
        self._mimicWeights()

        # reweight and recalculate angles
        self._ReweightAndTransformAngles(self._inTree, self._mimicedVars['inDistr'][self._vars[0]],
                                                       self._mimicedVars['outDistr'][self._vars[0]],)


    def _ReweightAndTransformAngles( self,t, pin, pout, Nbins= None ):
        """ t: TTree, pin: original momentum distribution (python list), pout : the momentum distribution you want (python list)
        Nbins controls the number of points for the transformation functions
        """
        if Nbins==None: Nbins=self._nBins

        # transformation of input and output distributions to uniform.
        Udat = UniFunc(pout, nbinsmax = Nbins)
        Umc = UniFunc(pin, nbinsmax = Nbins)

        trans = lambda x : Udat.inverse( Umc(x) )
        trans2 = lambda x,y : ( trans(x) , trans(y)  )

        # put the newly recalculated angles plus time in a RooDataSet.
        from ROOT import RooDataSet, RooArgSet, RooFit
        for obs in self._mcObsSet:
            if   obs.GetName()=='helcosthetaK': helcosthetaK = obs
            elif obs.GetName()=='helcosthetaL': helcosthetaL = obs
            elif obs.GetName()=='helphi':          helphi    = obs
            elif obs.GetName()=='time':             time     = obs
            elif obs.GetName()=='truetime':         time     = obs

        recalculatedVars = RooArgSet( [helcosthetaK,helcosthetaL,helphi] + self._KaonSet + self._BmomSet ) 
        copiedVars       = RooArgSet( self._muonSet + [time, self._physWeightsVar] )
        recalculatedData = RooDataSet( 'MomRewMC_%s_Iter'%self._itNum, 'MomRewMC_%s_Iter'%self._itNum, recalculatedVars         )
        copiedData       = RooDataSet( 'copiedData',                   'copiedData',                   self._inTree, copiedVars )

        from ROOT import TDatabasePDG
	MeV = 1000 # TDatabasePDG is in GeV, this is the factor needed to go to MeV.
        PDG = TDatabasePDG()
        Mmu = PDG.GetParticle('mu-').Mass()*MeV
        Mk  = PDG.GetParticle('K-').Mass()*MeV
        
        from ROOT import TVector3, TLorentzVector, HelicityAngles 
        from math import sqrt
        _VM2LV  = lambda v,m : TLorentzVector( v, sqrt( m*m + v.Mag2() ) ) # TVector3,mass to TLorentzVector 
        _E2V    = lambda entry, label : TVector3( getattr(entry,label+'_PX'),getattr(entry,label+'_PY'),getattr(entry,label+'_PZ')) # entry to TVenctor3
        _B3PMAG = lambda K1,K2,mu1,mu2: (K1+K2+mu1+mu2).Mag() # B momentum.
        _BPT    = lambda K1,K2,mu1,mu2: sqrt( (K1+K2+mu1+mu2).x()**2 + (K1+K2+mu1+mu2).y()**2 ) # B transverse momentum.
       
        # get Kaon variables references. 
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

        assert(self._BmomSet[0].GetName()=='B_P' and self._BmomSet[1].GetName()=='B_Pt')# important check
        BP, BPt = self._BmomSet[0], self._BmomSet[1]
        
        # from datetime import datetime
        # startTime = datetime.now()
        
        print 'P2VV - INFO: Recalculating decay angles after kinematic distributions matching.'
        for entry in t:
            k1_3P = _E2V(entry,'Kplus')
            k2_3P = _E2V(entry,'Kminus')
            mag1,mag2 = trans2( k1_3P.Mag(), k2_3P.Mag() )
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
     
        #print 'Recalculating loop took that much to complete  ', (datetime.now()-startTime)
        #assert False
        
        recalculatedData.merge( copiedData )
        self._recalculatedData = RooDataSet(recalculatedData.GetName(), recalculatedData.GetTitle(), 
                                            RooArgSet(recalculatedVars, copiedVars), 
                                            Import = recalculatedData, 
                                            WeightVar = (self._physWeightsVar, True)
                                            ) 

   
    def getDataSet(self, tree=False): 
        if tree:
            from ROOT import gROOT
            gROOT.cd('PyROOT:/')
            return self._recalculatedData.buildTree()
        else:    return self._recalculatedData


# class for building the pdf for the efficiency moments calculation  
class EfficiencyMomentsPdfBuilder():
    def __init__(self, **kwargs):
        # set global object name prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        setParNamePrefix( 'moms' )
        
        # get workspace
        from P2VV.RooFitWrappers import RooObject
        ws = RooObject().ws()
        
        # angular functions
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

        # get obeservables
        from P2VV.Utilities.MCReweighting import _createObservable
        for o in  ['helcosthetaK','helcosthetaL','helphi','time','truetime','iTag']:
            if not ws[o]: _createObservable(o)
        angles = [ ws[o] for o in ['helcosthetaK','helcosthetaL','helphi'] ]
        iTag = ws['iTag']
        obsSet = angles + [ws['truetime']]

        # pdf apramter valiues
        from math import sqrt, cos, sin
        pars = kwargs.pop('pdfParVals')
        AparMag2 = sqrt( 1 - pars['AperpMag2'] - pars['A0Mag2'] )

        # angular amplitude function
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
        amplitudes = Amplitudes(  ReApar  = sqrt( AparMag2 / pars['A0Mag2'])            * cos(pars['AparPhase']),
                                  ImApar  = sqrt(AparMag2  / pars['A0Mag2'])            * sin(pars['AparPhase']), 
                                  ReAperp = sqrt(pars['AperpMag2'] / pars['A0Mag2'])    * cos(pars['AperpPhase']),
                                  ImAperp = sqrt(pars['AperpMag2'] / pars['A0Mag2']) * sin(pars['AperpPhase']),
                                  ReAS    = 0.,
                                  ImAS    = 0.
                                  )
                                  
        from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        lifetimeParams = LifetimeParams( Gamma = pars['Gamma'], dGamma = pars['dGamma'], dM = pars['dM'] )
        
        tResArgs = { }
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs['time'] = ws['truetime']
        timeResModel = TimeResolution(  )

        # CP violation parameters
        from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
        lambdaCP = CPParam( lambdaCPSq = 1., phiCP = pars['phiCP'] )

        # tagging parameters
        from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
        taggingParams = TaggingParams()

        # coefficients for time functions
        from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

        # build underlying physics PDF
        args = dict(    time            = ws['truetime']
                      , iTag            = iTag
                      , tau             = lifetimeParams['MeanLifetime']
                      , dGamma          = lifetimeParams['dGamma']
                      , dm              = lifetimeParams['dM']
                      , dilution        = taggingParams['dilution']
                      , ADilWTag        = taggingParams['ADilWTag']
                      , avgCEven        = taggingParams['avgCEven']
                      , avgCOdd         = taggingParams['avgCOdd']
                      , coshCoef        = timeBasisCoefs['cosh']
                      , sinhCoef        = timeBasisCoefs['sinh']
                      , cosCoef         = timeBasisCoefs['cos']
                      , sinCoef         = timeBasisCoefs['sin']
                      , resolutionModel = timeResModel['model']
                      )

        from P2VV.RooFitWrappers import BTagDecay
        self._pdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )
        setParNamePrefix( '' )
        

    def getPdf(self): return self._pdf




class BuildBs2JpsiKK2011sFit():
    def __init__(self,**kwargs):
        print 'P2VV - INFO: Initialised physics reweighting class: buildBs2JpsiKK2011sFit().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
        pdfConfig = PdfConfig()

        self._dataSetPath = kwargs.pop('dataSetPath', None)
        self._dataSetName = kwargs.pop('dataSetName', None)
        self._weightsName = kwargs.pop('weightsName', 'JpsiKK_sigSWeight')
    
        self._doUntaggedFit = kwargs.pop('doUntaggedFit', '')        

        pdfConfig['timeEffHistFiles'] = kwargs.pop('timeEffHistFile', {})
    
        parFileIn  = kwargs.pop( 'parFileIn',  '' )
        parFileOut = kwargs.pop( 'parFileOut', '' )

        if self._doUntaggedFit: print 'P2VV - INFO: BuildBs2JpsiKK2011sFit: Building untagged fit.' 
        
        # fit options
        pdfConfig['fitOptions'] = kwargs.pop( 'fitOpts', None)
        if not pdfConfig['fitOptions']:
            pdfConfig['fitOptions'] = dict(    NumCPU    = 8
                                             , Optimize  = 2
                                             , Minimizer = 'Minuit2'
                                             , Offset    = True
                                             #               , Hesse     = False
                                             , Timer     = True
                                             #               , Verbose   = True
                                             )
        self._FitResults  = {} # collect all the fit results
        self._Moments     = {} # collect all moments
        self._PDFS        = {} # collect allpdfs

        # PDF options 
        KKmassParam = kwargs.pop( 'KKmassBins', None )
        pdfConfig['timeEffType']          = 'paper2012' if not self._doUntaggedFit else ''
        pdfConfig['anglesEffType']        = ''
        pdfConfig['KKMassBinBounds']      = [ 990., 1050. ] if not KKmassParam else [990., 1020. - 12., 1020., 1020. + 12., 1050.]
        pdfConfig['CSPValues']            = [ 0.498]        if not KKmassParam else [ 0.959, 0.770, 0.824, 0.968 ] 

        KKMassPars = pdfConfig['obsDict']['KKMass']
        pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                         , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )
        
        pdfConfig['constrainTagging']   = '' if not  self._doUntaggedFit else 'fixed'
        pdfConfig['timeResType']           = 'eventNoMean'
        pdfConfig['numTimeResBins']        = 40
        pdfConfig['lambdaCPParam'] = 'lambPhi'
        
        if self._doUntaggedFit:
            pdfConfig['externalConstr'] = {}
            pdfConfig['externalConstr']['dM'] = ( 17.63, 0.11 )
            pdfConfig['externalConstr']['timeResSigmaSF'] = (1.45,   0.06 )
        else:
            pdfConfig['externalConstr']['dM'] = ( 17.63, 0.11 )
            pdfConfig['externalConstr'].pop('betaTimeEff')

        from P2VV.Imports import extConstraintValues
        extConstraintValues.setVal( 'DM',      ( 17.63, 0.11 ) )
        extConstraintValues.setVal( 'P0OS',    (  0.392, 0.008, 0.392 ) )
        extConstraintValues.setVal( 'DelP0OS', (  0.0110, 0.0034 ) )
        extConstraintValues.setVal( 'P1OS',    (  1.000,  0.023  ) )
        extConstraintValues.setVal( 'DelP1OS', (  0.000,  0.001  ) )
        extConstraintValues.setVal( 'P0SS',    (  0.350, 0.017, 0.350 ) )
        extConstraintValues.setVal( 'DelP0SS', ( -0.019, 0.005   ) )
        extConstraintValues.setVal( 'P1SS',    (  1.00,  0.16    ) )
        extConstraintValues.setVal( 'DelP1SS', (  0.00,  0.01    ) )

        # read Data.
        if type(self._dataSetPath)==list: 
            dataSet = _combineDataSetParts(self._dataSetPath, self._dataSetName, weightName=self._weightsName)
        else: 
            from P2VV.Utilities.DataHandling import readData
            dataSet = readData( filePath = self._dataSetPath, dataSetName = self._dataSetName,  NTuple = False )
        pdfConfig['signalData'] = dataSet
        pdfConfig['runPeriods'] = []
        pdfConfig['readFromWS'] = True
        self._fitData = dataSet
        
        # build pdf.
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
        pdfConfig['parNamePrefix'] =  'data' # give the dataPdf parameters a prefix.
        self._pdfBuild = PdfBuilder( **pdfConfig )
        self._pdf = self._pdfBuild.pdf()
  
        print pdfConfig['fitOptions']
        if not 'Optimize' in pdfConfig['fitOptions'] or pdfConfig['fitOptions']['Optimize'] < 2 :
            # unset cache-and-track
            for par in self._pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

        if parFileIn :
            # read parameters from file
            pdfConfig.readParametersFromFile( filePath = parFileIn )
            pdfConfig.setParametersInPdf(self._pdf)
        
        from P2VV import RooFitDecorators
        # Taggging calibration for untagged pdf
        if  self._doUntaggedFit:
            self._pdf.ws()['data_wTagP1OS'].setMin(0)
            self._pdf.ws()['data_wTagP1SS'].setMin(0)

            self._pdf.ws()['data_wTagP0OS'].setVal(.5)
            self._pdf.ws()['data_wTagP0SS'].setVal(.5)
            self._pdf.ws()['data_wTagP1OS'].setVal(0)
            self._pdf.ws()['data_wTagP1SS'].setVal(0)
            self._pdf.ws()['data_wTagDelP0OS'].setVal(0)
            self._pdf.ws()['data_wTagDelP0SS'].setVal(0)
            
            self._pdf.ws()['data_wTagP1OS'].setConstant(True)
            self._pdf.ws()['data_wTagP1SS'].setConstant(True)
            self._pdf.ws()['data_wTagP0OS'].setConstant(True)
            self._pdf.ws()['data_wTagP0SS'].setConstant(True)
            self._pdf.ws()['data_wTagP1OS'].setConstant(True)
            self._pdf.ws()['data_wTagP1SS'].setConstant(True)
            self._pdf.ws()['data_wTagDelP0OS'].setConstant(True)
            self._pdf.ws()['data_wTagDelP0SS'].setConstant(True)    

            # These pars drop out in the mc pdf
            self._pdf.ws()['data_ASOddPhase_bin0'].setVal(0.)
            self._pdf.ws()['data_f_S_bin0'].setVal(0.)
            self._pdf.ws()['data_ASOddPhase_bin0'].setConstant(True)
            self._pdf.ws()['data_f_S_bin0'].setConstant(True)

            # These pars have little sensitivity 
            from P2VV.Utilities.MCReweighting import parValuesMc2012Fit as pars
            self._pdf.ws()['data_phiCP'].setVal(pars['phiCP'])
            self._pdf.ws()['data_phiCP'].setConstant()
            self._pdf.ws()['data_lambdaCP'].setVal(pars['lambdaCP'])
            self._pdf.ws()['data_lambdaCP'].setConstant()

        self._pdfConfig = pdfConfig


    def doFit( self, itNum=0, angAccFile=None ):
        pref = self._pdfConfig['parNamePrefix'] + '_'
        
        # multiply by angular acceptance
        if angAccFile: self._pdf = self._multiplyPdfWithAcc( angAccFile, iterNumb = itNum )

        # get observables and parameters in PDF
        pdfObs  = self._pdf.getObservables(self._pdfConfig['signalData'])
        pdfPars = self._pdf.getParameters(self._pdfConfig['signalData'])

        # create fit data
        self._fitData = self._pdfConfig['signalData'].reduce(pdfObs)
        self._pdfConfig['signalData'].IsA().Destructor(self._pdfConfig['signalData'])
        
        # float/fix values of some parameters
        for CEvenOdds in self._pdfBuild['taggingParams']['CEvenOdds'] :
            for CEvenOdd in CEvenOdds :
                CEvenOdd.setConstant( pref + 'avgCEven.*')
                CEvenOdd.setConstant( pref +  'avgCOdd.*', True )
        
        self._pdfBuild['tagCatsOS'].parameter(pref + 'wTagDelP1OS').setVal(0.)
        self._pdfBuild['tagCatsSS'].parameter(pref + 'wTagDelP1SS').setVal(0.)
        self._pdfBuild['tagCatsOS'].setConstant(pref + 'wTagDelP1')
        self._pdfBuild['tagCatsSS'].setConstant(pref + 'wTagDelP1')
        
        self._pdfBuild['amplitudes'].setConstant(pref + 'C_SP')


        print 120 * '='
        print 'Bs2JpsiKK2011Fit: fitting %d events (%s)' % ( self._fitData.numEntries(), 'weighted' if self._fitData.isWeighted() else 'not weighted' )

        fitResult = self._pdf.fitTo( self._fitData, SumW2Error = False, Save = True,  **self._pdfConfig['fitOptions'] )

        # print parameter values
        from P2VV.Imports import parNames, parValues2011 as parValues
        print 'Bs2JpsiKK2011Fit: parameters:'
        fitResult.SetName('sFit_%s_Iteration'%itNum)
        fitResult.PrintSpecial( text = False, LaTeX = False, normal = True, ParNames = parNames, ParValues = parValues )
        fitResult.covarianceMatrix().Print()
        fitResult.correlationMatrix().Print()
        from ROOT import TFile
        resultFile = TFile.Open('fitResult_iterativeProcedure_%s.root'%itNum,'recreate')
        resultFile.cd()
        fitResult.Write()
        resultFile.Close()
        self._FitResults['iter_%s'%itNum] = fitResult 
    
        print 120 * '=' + '\n'        

        if type(self._dataSetPath)==list: 
            self._pdfConfig['signalData'] = _combineDataSetParts(self._dataSetPath, self._dataSetName, weightName=self._weightsName)
        else: 
            from P2VV.Utilities.DataHandling import readData
            self._pdfConfig['signalData'] = readData(filePath=self._dataSetPath, dataSetName=self._dataSetName,  NTuple=False)



    def updateDataParameters(self, oldPars, itNum=0):
        fitResult = self._FitResults['iter_%s'%itNum]
        parkeys = oldPars.keys()
        parkeys.remove('A0Phase')
        for par in parkeys:
            try: oldPars[par] = fitResult.floatParsFinal().find( self._pdfConfig['parNamePrefix'] + '_' + par).getVal()
            except AttributeError: print 'P2VV - WARNING: updateDataParameters: Parameter %s not found in fit result %s using, parameter value will not change.'\
                    %( self._pdfConfig['parNamePrefix'] + '_' + par, fitResult.GetName() ) 
        
    def _multiplyPdfWithAcc( self, effFile, iterNumb=None ):
        # read moments file and multiply pure pdf with angular acceptance
        print 'P2VV - INFO:multiplyPdfWithAcc(): multiplying PDF "%s" with angular efficiency moments from file "%s"'\
          % ( self._pdf.GetName(), effFile )
        from P2VV.Utilities.DataMoments import angularMomentIndices, RealMomentsBuilder
        moments = RealMomentsBuilder()
        moments.appendPYList( self._pdfBuild['angleFuncs'].angles, angularMomentIndices( 'weights', self._pdfBuild['angleFuncs'] ) )
        moments.read( effFile )    
        self._Moments['%s_iteration'%iterNumb] = moments # collect all moments  
        if iterNumb: return moments.multiplyPDFWithEff( self._pdf, CoefName = 'effC%d' % iterNumb )
        else:        return moments.multiplyPDFWithEff( self._pdf, CoefName = 'effnom'            )
        
         
    def getDataSet(self):           return self._fitData
    def getPdf(self):               return self._pdf
    def getPdfBuilderObject(self) : return self._pdfBuild
    def getFitResults(self):        return self._FitResults
    def getMomentsDict(self):       return self._Moments['%s_iteration'%itNum]
    def getObservables(self,which=None):
        if   which=='angles': return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName().startswith('hel') ]
        elif which=='time':   return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName()==('time')         ]
        else:                 return             self._pdfBuild['obsSetP2VV']


         
parValuesNoKKBinsWideKKWindow = dict( AperpMag2        = 0.24663
                                      ,AperpPhase       = 3.032
                                      ,A0Mag2           = 0.52352
                                      ,A0Phase          = 0 # cosntrained
                                      ,AparPhase        = 3.2176
                                      ,f_S_bin0         = 0.046222
                                      ,ASOddPhase_bin0  = -0.066683
                                      ,dM               = 17.677
                                      ,dGamma           = 0.1023
                                      ,Gamma            = 0.67283
                                      ,phiCP            = 0.085636
                                      ,lambdaCP         = 0.92737
                                       )
parValuesMc2011Gen =  dict( AperpMag2        = 0.16
                           ,AperpPhase       = -0.17
                           ,A0Mag2           =  0.6
                           ,A0Phase          = 0 # cosntrained
                           ,AparPhase        = 2.50
                           ,f_S_bin0         = 0
                           ,ASOddPhase_bin0  = 0
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
                           ,f_S_bin0         = 0.
                           ,ASOddPhase_bin0  = 0.
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
                           
                           #'/project/bfys/vsyropou/data/P2VVDataSetsMC2012_wideKKMass_noKKMassBins_2TagCats_TrivialWeights_forReweighting.root'
                            AperpMag2        =  2.4592e-01
                           ,AperpPhase       =  2.9835e+00
                           ,A0Mag2           =  5.2006e-01
                           ,A0Phase          =  0. # cosntrained
                           ,AparPhase        =  2.9608e+00
                           ,f_S_bin0         =  0.
                           ,ASOddPhase_bin0  =  0.
                           ,dM               =  1.7804e+01
                           ,dGamma           =  9.2625e-02
                           ,Gamma            =  6.7736e-01
                           ,phiCP            = -1.1723e-02
                           ,lambdaCP         =  1.0354e+00 
                            )

