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
    mcData          = kwargs.pop('mcData',          None)
    momRewData      = kwargs.pop('momRewData',      None)
    sData           = kwargs.pop('sData',           None)
    itNumb          = kwargs.pop('itNumb',          None)
    physWeightsName = kwargs.pop('physWeightsName', None)
    nullTest        = kwargs.pop('nullTest',        None)

    # set truetime to time
    if mcData.get().find('truetime'): 
        truetimeFlag = True
        mcData.get().find('truetime').SetName('time')
        momRewData.get().find('truetime').SetName('time')

    # get mcData before and after physics reweighting
    from ROOT import RooFit, RooDataSet
    mcBeforePhysRew = mcData
    mcAfterPhysRew  = RooDataSet( 'MC_BeforePhysRew_%s'%itNumb, 'MC_BeforePhysRew_%s'%itNumb, mcData.get(), 
                                  RooFit.Import(mcData),
                                  RooFit.WeightVar(physWeightsName)
                                  )
   
    # get observables and x ranges
    observables, Kmomenta, muMomenta = [], [], []
    trackMomRangeX = {}
    for obs in mcData.get():
        obsName = obs.GetName()
        if 'hel' in obsName or 'time' in obsName: observables.append(obs)
        elif obsName.startswith('K') and not obsName.startswith('KK'): 
            if   'X' in obsName or 'Y' in obsName: trackMomRangeX[obsName] = (-5e3,5e3)
            elif 'Z' in obsName:                   trackMomRangeX[obsName] = (-5e2,1e5)
            else:                                  trackMomRangeX[obsName] = (0,1e5)
            Kmomenta.append(obs)
        elif obsName.startswith('mu'):  
            if   'X' in obsName or 'Y' in obsName: trackMomRangeX[obsName] = (-1e4,1e4)
            elif 'Z' in obsName:                   trackMomRangeX[obsName] = (-1e3,2e5)
            else:                                  trackMomRangeX[obsName] = (0,2e5)
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
    # stdDrawOpts = dict( DataError=RooAbsData.SumW2, MarkerSize = .5, XErrorSize = 0 )
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
    legend.AddText('#color[%s]{%s}'%( colors['mcBefore'], 'SourceBeforePhysRew' ) )
    legend.AddText('#color[%s]{%s}'%( colors['mcAfter'], 'SourceAfterPhysRew' ) )
    legend.AddText('#color[%s]{%s}'%( colors['MomRewData'], 'SoourceAfterMomRew' ) )
    legend.AddText('#color[%s]{%s}'%( colors['Sdata'], 'Target' ) )

    for canv, pad in zip( [obsCanv, KaonCanv, muonCanv, assymKaonCanv, assymMuonCanv], [4,8,8,8,8] ): 
        canv.cd(pad)
        legend.Draw()
    _P2VVPlotStash.append(legend) # keep reference to the legend
    
    if truetimeFlag: # restore the truetime in mcData 
        mcData.get().find('time').SetName('truetime')
        momRewData.get().find('time').SetName('truetime')
    
    del mcAfterPhysRew
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
def combineDataSetParts( files, name, weightName='' ):
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


# Vertical reweighting class of MC to match the physics of sWeighted data.
class MatchPhysics( dict ):
    def __init__( self,nTupleFile, **kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'
        self._nTupleFile = nTupleFile
        self._nTupleName = kwargs.pop('nTupleName', 'DecayTree')
        self._mcPars     = kwargs.pop('mcParameters',    {}    ) 
        
        self._trueTime    = kwargs.pop('trueTime',  True  ) 
        self._modelSwave  = kwargs.pop('modelSwave', True ) 
        self._selectData  = kwargs.pop('selectData', False)
        
        self._multiplyWithTimeAcc = kwargs.pop('multiplyWithTimeAcc', False)
        self._timeEffFiles        = kwargs.pop('timeEffFiles','')
        self._timeEffHistNames    = kwargs.pop('timeEffHistNames','')

        self._allWeights = {} # collect all the weights throughout the loop
        self._pref       = 'mc' # prefix for all objects created by this class
        
        # set the prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        setParNamePrefix( self._pref )

    def buildMonteCarloPdf( self, dataPdfBuilder=None ):
        if not dataPdfBuilder: assert False, 'P2VV - ERORR: Build data pdf first and provide the pdf builder object.'
        from math import pi, sin, cos, sqrt

        # job parameters
        tResModel   = ''
        trigger     = ''
        timeInt     = False

        # transversity amplitudes
        A0Mag2Val    = self._mcPars['A0Mag2']    # 0.60
        AperpMag2Val = self._mcPars['AperpMag2'] #  0.16
        AparMag2Val  = 1. - A0Mag2Val - AperpMag2Val
        A0PhVal    = self._mcPars['A0Phase'] # 0.
        AperpPhVal = self._mcPars['AperpPhase'] #-0.17
        AparPhVal  = self._mcPars['AparPhase'] #2.50
    
        # CP violation parameters
        phiCPVal      = self._mcPars['phiCP'] #-0.04

        # B lifetime parameters
        GammaVal  = self._mcPars['Gamma'] #0.679
        dGammaVal = self._mcPars['dGamma'] #0.060
        dMVal     = self._mcPars['dM'] #17.8
        tResSigma = 0.045

        ####################################
        ## create variables and read data ##
        ####################################

        # import RooFit wrappers
        from P2VV.Load import RooFitOutput

        # get / create observables.
        self._mcObsSet = [ obs for obs in dataPdfBuilder['obsSetP2VV'] ]        
        for obs in self._mcObsSet: 
            if obs.GetName()   == 'time':        time = obs
            elif obs.GetName() =='helcosthetaK': helcosthetaK = obs
            elif obs.GetName() =='helcosthetaL': helcosthetaL = obs
            elif obs.GetName() =='helphi':       helphi = obs 

        from P2VV.RooFitWrappers import RealVar, Category  
        trueTime = RealVar( 'truetime', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0., 20. ))
        iTag = Category('iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )

        self._mcObsSet.remove( time )
        self._mcObsSet.append( trueTime )
        self._mcObsSet.append( iTag )

        # angular functions
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        self._angleFuncs = AngleFuncs( cpsi = helcosthetaK, ctheta = helcosthetaL, phi = helphi )
       
        # Read data        
        self._createNtupleVars()
        if self._selectData: self._selectMonteCarloData()
        else: 
            if type(self._nTupleFile) == list: self._initData = combineDataSetParts(self._nTupleFile, self._nTupleName)
            else:
                from ROOT import TFile
                self._initData = TFile.Open(self._nTupleFile).Get(self._nTupleName)
        

        #####################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################
        
        # helicity  amplitudes
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
        amplitudes = Amplitudes(    ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
                                  , ImApar  = sqrt(AparMag2Val  / A0Mag2Val) * sin(AparPhVal)
                                  , ReAperp = sqrt(AperpMag2Val / A0Mag2Val) * cos(AperpPhVal)
                                  , ImAperp = sqrt(AperpMag2Val / A0Mag2Val) * sin(AperpPhVal)
                                  , ReAS    = 0.
                                  , ImAS    = 0.
                                  )
        # B lifetime
        from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        lifetimeParams = LifetimeParams( Gamma = GammaVal, dGamma = dGammaVal, dM = dMVal )

        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
        tResArgs = { }
        tResArgs['time'] = trueTime
        timeResModel = TimeResolution( **tResArgs )

        # CP violation parameters
        from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
        lambdaCP    = CPParam( lambdaCP   = 1., phiCP = phiCPVal )
        lambdaCPVal = lambdaCP._lambdaCP.getVal()
      
        # tagging parameters
        from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
        taggingParams = TaggingParams()

        # coefficients for time functions
        from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        index = [ 'A0', 'Apar', 'Aperp', 'AS' ] if self._modelSwave else [ 'A0', 'Apar', 'Aperp' ]
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, amplitudes, lambdaCP, index )
        
        # build underlying physics PDF
        args = dict(    time            = trueTime
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
        self._pdf =  BTagDecay(self._pref + '_sig_t_angles_tagCat_iTag', **args )

        print self._pdf.ws()[self._pref + '_sig_t_angles_tagCat_iTag'].getVal() 
        
        if self._multiplyWithTimeAcc: 
            #from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder
            #Bs2Jpsiphi_PdfBuilder._multiplyByTimeAcceptance( self._pdf,  data  = self._initData, 
            #                                       timeEffHistFiles = self._timeEffFiles,
            #                                      timeEffType = 'paper2012'
            #                                             )
            obsDict = dict(helcosthetaK  = helcosthetaK, 
                           helcosthetaL  = helcosthetaL, 
                           helphi        = helphi, 
                           time          = trueTime, 
                           hlt1ExclB     = self._initData.get().find('hlt1_excl_biased_dec'),
                           hlt2B         = self._initData.get().find('hlt2_biased'),
                           hlt2UB        = self._initData.get().find('hlt2_unbiased')
                           )
        
            from P2VV.Parameterizations.FullPDFs import multiplyByTimeAcceptance
            timeAccArgs = dict( histFile           = self._timeEffFiles,
                                histUBName         = self._timeEffHistNames['UB'], histExclBName = self._timeEffHistNames['ExclB'],
                                timeEffParameters  = {},
                                timeEffType        = 'paper2012',
                                data               = self._initData, 
                                motherPdf          = self._pdf,
                                observables        = obsDict,
                                timeResModels      = timeResModel['model'],
                                timeResModelsOrig  = timeResModel['model'],
                                parNamePrefix      = self._pref
                                )
 
            multiplyByTimeAcceptance( self._pdf, self, **timeAccArgs )

        print self._pdf.ws()[self._pref + '_sig_t_angles_tagCat_iTag'].getVal() 
        
        self._AngAmpsParsVals = dict(A0Mag2     = A0Mag2Val, 
                                     AperpMag2  = AperpMag2Val,
                                     AparMag2   = AparMag2Val,
                                     A0Ph       = A0PhVal,
                                     AperpPh    = AperpPhVal,
                                     AparPh     = AparPhVal,
                                     phiCP      = phiCPVal,
                                     Gamma      = GammaVal,
                                     dGamma     = dGammaVal,
                                     dM         = dMVal,
                                     #tResSigma  = tResSigma,
                                     lambdaCP   = lambdaCPVal 
                                     )     
        
    def _createNtupleVars(self):
        # NOTE: All these objects have a reason to be out of P2VV, mainly to avoid conflicts with the K*_P
        # ntuple variables
        from ROOT import RooRealVar, RooNumber
        RooInf  = RooNumber.infinity()
        B_P        = RooRealVar( 'B_P',        'B_P',         0  , RooInf )
        B_PT       = RooRealVar( 'B_Pt',       'B_Pt',        0  , RooInf )   
        Kplus_P    = RooRealVar( 'Kplus_P',    'Kplus_P',     0  , RooInf )  
        Kplus_PX   = RooRealVar( 'Kplus_PX',   'Kplus_PX',   -RooInf, RooInf )   
        Kplus_PY   = RooRealVar( 'Kplus_PY',   'Kplus_PY',   -RooInf, RooInf )   
        Kplus_PZ   = RooRealVar( 'Kplus_PZ',   'Kplus_PZ',   -RooInf, RooInf )  
        Kminus_P   = RooRealVar( 'Kminus_P',   'Kminus_P',    0  , RooInf )  
        Kminus_PX  = RooRealVar( 'Kminus_PX',  'Kminus_PX',  -RooInf, RooInf )  
        Kminus_PY  = RooRealVar( 'Kminus_PY',  'Kminus_PY',  -RooInf, RooInf )  
        Kminus_PZ  = RooRealVar( 'Kminus_PZ',  'Kminus_PZ',  -RooInf, RooInf )  
        muplus_P   = RooRealVar( 'muplus_P',   'muplus_P',    0  , RooInf )  
        muplus_PX  = RooRealVar( 'muplus_PX',  'muplus_PX',  -RooInf, RooInf )  
        muplus_PY  = RooRealVar( 'muplus_PY',  'muplus_PY',  -RooInf, RooInf )   
        muplus_PZ  = RooRealVar( 'muplus_PZ',  'muplus_PZ',  -RooInf, RooInf )   
        muminus_P  = RooRealVar( 'muminus_P',  'muminus_P',   0  , RooInf )   
        muminus_PX = RooRealVar( 'muminus_PX', 'muminus_PX', -RooInf, RooInf )   
        muminus_PY = RooRealVar( 'muminus_PY', 'muminus_PY', -RooInf, RooInf )   
        muminus_PZ = RooRealVar( 'muminus_PZ', 'muminus_PZ', -RooInf, RooInf )   
        
        self._ntupleVars = [ Kplus_P, Kplus_PX, Kplus_PY, Kplus_PZ, Kminus_P, Kminus_PX, Kminus_PY, Kminus_PZ,\
                             muminus_P, muminus_PX, muminus_PY, muminus_PZ, muplus_P, muplus_PX, muplus_PY, muplus_PZ, B_P, B_PT ]

    
    def selectMonteCarloData(self):     
        bkgcatCut      = '(bkgcat == 0 || bkgcat == 50)'
        trackChiSqCuts = 'muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
        massCuts       = 'mass > 5200. && mass < 5550. && mdau1 > 3030. && mdau1 < 3150. && mdau2 > 990. && mdau2 < 1050.'
        timeCuts       = 'time > 0.3 && time < 14. && sigmat < 0.12'
        tagCuts        = '(tagdecision == 0 || tagdecision == -1 || tagdecision == +1)'
        cuts = bkgcatCut + ' && ' + trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
        cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
        
        from ROOT import RooFit, RooDataSet, RooArgSet, TFile       
        print 'P2VV - INFO: reading NTuple "%s" from file "%s"' % ( self._nTupleFile, self._nTupleName )
        mcT = TFile.Open(self._nTupleFile).Get(self._nTupleName)
        print 'P2VV - INFO: applying cuts: %s' % cuts
        junkFile = TFile.Open('/data/bfys/vsyropou/junk.root','recreate')
        junkFile.cd()
        mcT = mcT.CopyTree(cuts)
     
        self._initData = RooDataSet( self._nTupleName, self._nTupleName, RooArgSet(self._mcObsSet + self._ntupleVars), RooFit.Import(mcT) )


    def setMonteCarloParameters(self, pars=None):
        if not pars:
            mcPars = self._AngAmpsParsVals
            from math import sin, cos, sqrt
            pars = dict(  ReAperp    = sqrt( mcPars['AperpMag2'] /  mcPars['A0Mag2'] ) * cos(  mcPars['AperpPh'] ),
                          ImAperp    = sqrt( mcPars['AperpMag2'] /  mcPars['A0Mag2'] ) * sin(  mcPars['AperpPh'] ),
                          ReApar     = sqrt( mcPars['AparMag2']  /  mcPars['A0Mag2'] ) * cos(  mcPars['AparPh']  ),
                          ImApar     = sqrt( mcPars['AparMag2']  /  mcPars['A0Mag2'] ) * sin(  mcPars['AparPh']  ),
                          ReA0       = cos(  mcPars['A0Ph'] ),
                          ImA0       = sin(  mcPars['A0Ph'] ),
                          ReAS        = 0.,
                          ImAS        = 0.,
                          dM         = mcPars['dM'],
                          dGamma     = mcPars['dGamma'],
                          Gamma      = mcPars['Gamma'],
                          phiCP      = mcPars['phiCP'],
                          lambdaCP   = mcPars['lambdaCP'] 
                          )
        print 'P2VV - INFO: setMonteCarloParameters: Setting the following parameters as MC parameters.'
        displayOnlyDict = {}
        for var in self._pdf.Parameters():
            var.setVal( pars[ var.GetName().partition('_')[2] ] ) # remove prefix from the key. 
            displayOnlyDict[var.GetName()] = var.getVal()
        print displayOnlyDict
            
    def setDataFitParameters(self, dataPars, KKmassCat=None):
        print 'P2VV - INFO: setDataFitParameters: Setting the following parameters as data parameters.'
        displayOnlyDict = {}
        from math import sqrt,sin, cos
        if not KKmassCat:
            AparMag2 = 1. - dataPars['A0Mag2'] - dataPars['AperpMag2']
            ASMag2   = dataPars['f_S_bin0'] / (1 - dataPars['f_S_bin0'])
            for par in self._pdf.Parameters():
                name = par.GetName()
                if 'Re' in name :
                    if   'Aperp' in name: par.setVal(  sqrt(dataPars['AperpMag2']/dataPars['A0Mag2']) * cos(dataPars['AperpPhase'])  )
                    elif 'Apar'  in name:  par.setVal(  sqrt(     AparMag2 / dataPars['A0Mag2']      ) * cos(dataPars['AparPhase'] )  )
                    elif  'AS'   in name:    par.setVal(  sqrt(       ASMag2 / dataPars['A0Mag2']      ) * cos(dataPars['ASOddPhase_bin0'])  )
                    elif  'A0'   in name:    par.setVal(  cos(dataPars['A0Phase'])                                                      )                   
                elif 'Im' in name :
                    if  'Aperp' in name: par.setVal( sqrt(dataPars['AperpMag2']/dataPars['A0Mag2']) * sin(dataPars['AperpPhase'])  )
                    elif 'Apar' in name:  par.setVal( sqrt(     AparMag2/dataPars['A0Mag2']        ) * sin(dataPars['AparPhase'] )  )
                    elif 'AS'   in name:    par.setVal( sqrt(       ASMag2 / dataPars['A0Mag2']      ) * sin(dataPars['ASOddPhase_bin0'])  )
                    elif 'A0'   in name:    par.setVal( sin(dataPars['A0Phase'])                                                      )
                else: par.setVal( dataPars[ name.partition('_')[2] ] ) # remove prefix from the key.  
                displayOnlyDict[par.GetName()] = par.getVal()
            print displayOnlyDict
        else: print ' ' # TODO: LoopOver the KKbins and set f_S_i and ASOddPhase_i


    def calculateWeights(self, iterNumb, dataParameters):
        self._iterNumb = iterNumb

        from ROOT import RooArgSet
        normVars =  RooArgSet(obs._target_() for obs in self._mcObsSet)
        
        # Reweights MC verticaly to match the Physics of data.
        nominators, denominators, self._weights = [], [], []
            
        self._pdf.attachDataSet( self._initData ) # make the pdf directly dependant on data
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'    
        
        # if iterNumb==1: self.setMonteCarloParameters()
        #else: pass          
        #for var in self._pdf.Parameters(): var.Print() 

        self.setMonteCarloParameters()
        for nev in xrange(self._initData.numEntries()):
            self._initData.get(nev)
            denominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating nominators for phyisics matching weights'
        self.setDataFitParameters(dataParameters)
        #for var in self._pdf.Parameters(): var.Print()
        for nev in xrange(self._initData.numEntries()):
            self._initData.get(nev)
            nominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating phyisics matching weights'
        for n,d in zip(nominators,denominators): self._weights += [n/d]
        

    def writeWeights(self, weightsName='weightPhys'):
        from ROOT import RooArgSet,RooRealVar,RooDataSet,RooNumber
        from P2VV.RooFitWrappers import RealVar
        RooInf  = RooNumber.infinity()
        weightsVar     = RooRealVar( weightsName,  weightsName, -RooInf, RooInf   )
        weightsArgSet  = RooArgSet( weightsVar )
        weightsDataSet = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        
        for weight in self._weights:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        
        # remove previous weights column
        if self._initData.get().find(weightsName): 
            self._initData.get().remove( self._initData.get().find(weightsName) )
            self._initData = self._initData.reduce( self._initData.get() )

        self._initData.merge( weightsDataSet )
        self._initData.SetName('MC_physicsReweighted_%s_iteration'%self._iterNumb )
        print 'P2VV - INFO: Phyisics matching weights added to dataset: '+'MC_physicsReweighted_%s_iteration'%self._iterNumb

        self._weightsName = weightsName
        del self._weights
        del weightsDataSet
    

    def getPdf(self):               return self._pdf
    def getNtupleVars(self):        return self._ntupleVars
    def getAngleFunctions(self):    return self._angleFuncs
    def getDataset(self):           return self._initData
    def getAllWeights(self):        return self._allWeights
    def getMcObsSet(self)  :
        return  [o for o in self._mcObsSet if o.GetName().startswith('hel') or 'time' in o.GetName()]



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
            self._inTree  = self._inTree.buildTree()
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
            

    def mimicWeights(self):
        self._mimicedVars = dict( inDistr={} )
        for var in self._vars:
            print 'P2VV - INFO: Input Distribution: Mimicing weights of variable: ', var
            self._mimicedVars['inDistr'][var]  = self.MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins, 'in')

            if self._mimicedVars.has_key('outDistr'):
                print 'Entrer first if statement'
                print self._mimicedVars.kyes()
                if not self._mimicedVars['outDistr'].has_key(var):
                    print 'Enter second if statement'
                    print self._mimicedVars.keys()
                    print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
                    self._mimicedVars['outDistr'][var] = self.MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out') 

            else:
                print 'Enter else statement'
                print self._mimicedVars.keys()
                print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
                self._mimicedVars['outDistr'] = {}
                self._mimicedVars['outDistr'][var] = self.MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins, 'out')            


    def MimicWeightedDistribution( self, t, var, wPref, Nbins, varSpec ):
        #List of new mimiced distribution
        newDistribution = []

        #Create Binning
        varRange = [t.GetMinimum(var), t.GetMaximum(var)]
        binWidth = float( varRange[1]-varRange[0] )/ float(Nbins)

        # Calculate Bin boundaries
        bounds, binning = [],{}
        low_bin  = [varRange[0]]
        high_bin = [varRange[0] + binWidth]
        bounds   = [ [varRange[0], varRange[0]+binWidth] ]
        for b in xrange(Nbins-1):
            low_bin  = bounds[b][0] + binWidth
            high_bin = bounds[b][1] + binWidth
            bounds.append([low_bin,high_bin])

        # Create entry lists for the bins and prepare sumW dictionary
        from ROOT import gDirectory, TFile # I am not sure if a file is needed to create entrylists
        file = TFile.Open('junkDeleteThis.root','recreate')
        file.cd()
        entryLists, sumW = {},{}
        for b in xrange(Nbins):
            t.Draw('>>elist'+str(b), var+'>'+ str(bounds[b][0]) +'&&'+ var +'<='+ str(bounds[b][1])    )
            entryLists.update( {'bin'+str(b):gDirectory.Get('elist'+str(b))} )            
            sumW.update({ 'bin'+str(b) : {'bounds':bounds[b],'sumW':0}  })  

        # Loop over tree and sum the weights
        for idx in xrange(t.GetEntries()):
            t.GetEntry(idx)
            for k in entryLists.keys():
                if entryLists[k].Contains(idx):
                    sumW[k]['sumW']+=getattr(t,wPref)
                    break
        # Replace sWeighted distribution with an equivalent one
        # by binning and generating n = sumOfWeights random numbers per bin
        from ROOT import TRandom
        rdm = TRandom()
        for b in sumW.keys():
            for evnt in xrange( int(round(sumW[b]['sumW'])) ):
                number = rdm.Uniform( sumW[b]['bounds'][0], sumW[b]['bounds'][1] ) 
                newDistribution.append(number)
                
        # compare before and after mimic
        from ROOT import TH1F
        if 'in' in varSpec:
            hist = TH1F( varSpec + var + str(self._itNum), 
                         varSpec + var + str(self._itNum), 
                         100, self._inTree.GetMinimum(var), self._inTree.GetMaximum(var) )
            for i in newDistribution: hist.Fill(i)
            self._inTree.Draw(var, self._inWeightName)
            hist.Draw('same err')
            hist.Scale( self._inSumW /  hist.GetEntries() )
        if 'out' in varSpec:
            hist = TH1F( varSpec + var + str(self._itNum), 
                         varSpec + var + str(self._itNum), 
                         100, self._outTree.GetMinimum(var), self._outTree.GetMaximum(var) )
            for i in newDistribution: hist.Fill(i)
            self._outTree.Draw(var, self._outWeightName)
            hist.Draw('same err')
            hist.Scale( self._outSumW /  hist.GetEntries() )
        _P2VVPlotStash.append(hist)
               
        return newDistribution


    def reweight( self, itNum, data ):
        self._itNum = itNum
        self._physWeightsVar = data.get().find( self._inWeightName ) # get physics weights var   

        # get rid of the previous tree
        if self._inTree: del self._inTree
       
        # get dataset scales
        self._inSumW = data.sumEntries()
        self._dataSetsScale = self._outSumW / self._inSumW
        self._physRewScale  = self._inSumW / data.numEntries()
        assert self._physRewScale <1.01, 'P2VV - WARNING: Scale between original MC and physics reweighted MC is %s, need to remove the weights and then perform the  momentum reweighting.' %self._physRewScale 
        
        # convert roodataset to tree
        from ROOT import gROOT
        gROOT.cd('PyROOT:/')
        self._inTree = data.buildTree() # conver dataset to tree
        
        # mimic the weights 
        self.mimicWeights()

        # reweight and recalculate angles
        self.ReweightAndTransformAngles(self._inTree, self._mimicedVars['inDistr'][self._vars[0]],
                                                      self._mimicedVars['outDistr'][self._vars[0]],)


    def ReweightAndTransformAngles( self,t, pin, pout, Nbins= None ):
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
                                            RooFit.Import(recalculatedData), 
                                            RooFit.WeightVar(self._physWeightsVar) ) 

   
    def getDataSet(self, tree=False): 
        if tree:
            from ROOT import gROOT
            gROOT.cd('PyROOT:/')
            return self._recalculatedData.buildTree()
        else:    return self._recalculatedData



class BuildBs2JpsiKK2011sFit():
    def __init__(self,**kwargs):
        print 'P2VV - INFO: Initialised physics reweighting class: buildBs2JpsiKK2011sFit().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
        pdfConfig = PdfConfig()

        self._dataSetPath = kwargs.pop('dataSetPath', None)
        self._dataSetName = kwargs.pop('dataSetName', None)
        self._weightsName = kwargs.pop('weightsName', 'JpsiKK_sigSWeight')
    
        self._doUntaggedFit = kwargs.pop('doUntaggedFit', '')        

        pdfConfig['timeEffHistFiles'] = dict(  file      = kwargs.pop('timeEffHistFile')
                                             , hlt1UB    = kwargs.pop('timeEffHistUBName')
                                             , hlt1ExclB = kwargs.pop('timeEffHistExclBName')
                                            )
    
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
            dataSet = combineDataSetParts(self._dataSetPath, self._dataSetName, weightName=self._weightsName)
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
            self._pdfConfig['signalData'] = combineDataSetParts(self._dataSetPath, self._dataSetName, weightName=self._weightsName)
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

