###########################################################################################################################################
## Utilities.MCReweighting: P2VV utilities for reweighting Monte Carlo data to create desired distributions                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   VS, Vasilis Syropoulos, Nikhef, v.syropoulos@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################

def CompareWeightedDistributions(tree, sTree, var, **kwargs):
    sVar      =  kwargs.pop('sVar',      None        )
    cut       =  kwargs.pop('cut',       None        )
    sCut      =  kwargs.pop('sCut',      None        )
    weight    =  kwargs.pop('weight',    None        )
    sWeight   =  kwargs.pop('sWeight',   None        )
    rangeX    =  kwargs.pop('rangeX',    None        )
    bins      =  kwargs.pop('bins',      100         )
    assymPlot =  kwargs.pop('assymPlot', False       )
    save      =  kwargs.pop('Save',      [False,'_'] )

    from ROOT import RooDataSet
    if type(tree)  == RooDataSet: tree  = tree.buildTree()
    if type(sTree) == RooDataSet: sTree = sTree.buildTree()

    if rangeX:
        Xmin=str(rangeX[0])
        Xmax=str(rangeX[1])
    else:
        Xmin= str(min(tree.GetMinimum(var),sTree.GetMinimum(var)))
        Xmax= str(max(tree.GetMaximum(var),sTree.GetMaximum(var)))

    from ROOT import gPad
    if cut:
        if weight: tree.Draw(var + '>>hm('+str(bins)+','+Xmin+','+Xmax+')', weight +'*'+'('+cut+')' )
        else     : tree.Draw(var + '>>hm('+str(bins)+','+Xmin+','+Xmax+')',                 cut     )
    else         :
        print 'P2VV - INFO: Ploting first distribution (1st arguement) with weight: ' + weight + '.'
        if weight: tree.Draw(var + '>>hm('+str(bins)+','+Xmin+','+Xmax+')', weight)
        else:      tree.Draw(var + '>>hm('+str(bins)+','+Xmin+','+Xmax+')',    '' )
    hm = gPad.GetPrimitive('hm')

    if not sVar: sVar=var
    if sCut:
        if sWeight: sTree.Draw(sVar + '>>hs('+str(bins)+','+Xmin+','+Xmax+')', sWeight +'*'+'('+sCut+')', 'err' )
        else      : sTree.Draw(sVar + '>>hs('+str(bins)+','+Xmin+','+Xmax+')',                  sCut    , 'err' )
    else          :
        print 'P2VV - INFO: Ploting second distribution (2nd arguement) with weight: ' ,  sWeight
        if sWeight: sTree.Draw(sVar + '>>hs('+str(bins)+','+Xmin+','+Xmax+')', sWeight, 'err' )
        else:       sTree.Draw(sVar + '>>hs('+str(bins)+','+Xmin+','+Xmax+')',    ''  , 'err' )
    hs = gPad.GetPrimitive('hs')

    hm.SetFillColor(2)
    hm.SetStats(0)
    hm.SetTitle(var)

    hs.SetMarkerStyle(20)
    hs.SetMarkerSize(.5)
    hs.SetTitle(var)
    hs.SetStats()

    def getSumOfWeights(t,pref,cut):
        ## TODO: should return sumW of the selected events by the cut string.
        if cut: print 'WARNING: Returned number is sumW of the entire tree not the of the subseset selected by cut. '
        if pref=='': return t.GetEntries(cut)
        else:
            sumW=0
            for e in t:sumW+=getattr(e,pref)
            return sumW

    if cut==None:  cut=''
    if sCut==None: sCut=''
    if weight==None: weight=''
    if sWeight==None: sWeight=''
    if tree.GetEntries(cut)>sTree.GetEntries(sCut): hm.Scale(getSumOfWeights(sTree,sWeight,sCut) / getSumOfWeights(tree,weight,cut))
    else:                                           hs.Scale(getSumOfWeights(tree,weight,cut)    / getSumOfWeights(sTree,sWeight,sCut))

    if rangeX: hm.GetXaxis().SetRangeUser(rangeX[0], rangeX[1])
    if hm.GetMaximum() < hs.GetMaximum(): hm.GetYaxis().SetRangeUser(0, int( hs.GetMaximum() + .08* hs.GetMaximum() ))

    from ROOT import TCanvas
    c_distr = TCanvas(var,var)
    if assymPlot:
        from ROOT import TH1F, TMath
        c_asymt = TCanvas('asymetryPlot','asymetryPlot')
        asymPlot = TH1F('asymPlot','Assymetry Plot', bins, float(Xmin), float(Xmax))
        for b in xrange(1,hm.GetNbinsX()):
            try:asym=(hm.GetBinContent(b) - hs.GetBinContent(b)) / (hm.GetBinContent(b) + hs.GetBinContent(b))
            except ZeroDivisionError: asym=0
            ##TODO:: Impliment the errors
            #error = TMath.sqrt ( 1/hm.GetSumOfWeights() + 1/hs.GetSumOfWeights() )
            asymPlot.SetBinContent(b,asym)
            #asymPlot.SetBinError(b,error)
            c_asymt.cd()
            asymPlot.SetStats(0)
            asymPlot.Draw()

            c_distr.cd()
            hm.Draw()
            hs.Draw('same')
            if save[0]:
               c_distr.SaveAs('comp_' + save[1])
               asymPlot.SaveAs('assym_'  + save[1])
            return c_distr, asymPlot
    else:
        c_distr.cd()
        hm.Draw()
        hs.Draw('same')
        if save[0]: c_distr.SaveAs('comp_' + save[1])
        return c_distr

def HelicityAngles(**kwargs):
    # Calculation based on the ANA-2012-067-v3
    k1_P = kwargs.pop('Kpl_P')
    k2_P = kwargs.pop('Kmi_P')
    m1_P = kwargs.pop('mPl_P')
    m2_P = kwargs.pop('mMi_P')

    # Bs, KK, mm momenta 4 vectors 
    KK_P   = k1_P + k2_P 
    mm_P   = m1_P + m2_P
    KKmm_P = KK_P + mm_P

    # Unit vector along mumu direction in the KK mass r.f.
    m1_P.Boost( - KK_P.BoostVector() )
    m2_P.Boost( - KK_P.BoostVector() )
    e_KK = - (m1_P + m2_P).Vect().Unit()
    # Boost the muons back to lab frame
    m1_P.Boost( KK_P.BoostVector() )
    m2_P.Boost( KK_P.BoostVector() )

    # Unit vector along KK direction in the mm mass r.f.
    k1_P.Boost( - mm_P.BoostVector() )
    k2_P.Boost( - mm_P.BoostVector() )
    e_mm = - (k1_P+k2_P).Vect().Unit()
    # Boost the Kaons back to lab frame
    k1_P.Boost( mm_P.BoostVector() )
    k2_P.Boost( mm_P.BoostVector() )

    # Unit vector along KK direction in the mm mass r.f.
    k1_P.Boost( - KKmm_P.BoostVector() )
    k2_P.Boost( - KKmm_P.BoostVector() )
    m1_P.Boost( - KKmm_P.BoostVector() )
    m2_P.Boost( - KKmm_P.BoostVector() )
    e_KKmm = (m1_P + m2_P).Vect().Unit()

    # Perpenticular vectors to KK and mm planes in the KKmmm r.f.
    eta_KK = ( k1_P.Vect().Cross( k2_P.Vect()) ).Unit()
    eta_mm = ( m1_P.Vect().Cross( m2_P.Vect()) ).Unit()

    k1_P.Boost( KKmm_P.BoostVector() )
    k2_P.Boost( KKmm_P.BoostVector() )
    m1_P.Boost( KKmm_P.BoostVector() )
    m2_P.Boost( KKmm_P.BoostVector() )

    # Helicity angles. 
    from math import asin, acos, pi    
    k1_P.Boost( - KK_P.BoostVector() )
    m1_P.Boost( - mm_P.BoostVector() )

    costhetaK = ( k1_P.Vect().Unit() ).Dot(e_KK)
    costhetaL = ( m1_P.Vect().Unit() ).Dot(e_mm)
    
    cosphi = eta_KK.Dot(eta_mm)
    sinphi = eta_KK.Cross(eta_mm).Dot(e_KKmm)
    if sinphi>0: Phi = + acos( eta_KK.Dot(eta_mm) )
    else       : Phi = - acos( eta_KK.Dot(eta_mm) )

    return costhetaK, costhetaL, Phi


class UniFunc:
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


# Vertical reweighting class of MC to match the physics of sWeighted data.
class MatchMCphysics2Data():
    def __init__( self,nTupleFile, **kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'
        self._pref       = kwargs.pop('ParNamePrefix','')
        self._nTupleFile = nTupleFile
        self._nTupleName = kwargs.pop('nTupleName', 'DecayTree')
        self._allWeights = {}
        
        self._trueTime   = kwargs.pop('trueTime',  True) 
        self._modelSwave = kwargs.pop('modelSwave', True) 

    def buildMonteCarloPdf( self, dataPdfBuilder=None ):
        if not dataPdfBuilder: assert False, 'P2VV - ERORR: Build data pdf first and provide the pdf builder object.'
        from math import pi, sin, cos, sqrt

        # job parameters
        tResModel   = ''
        trigger     = ''
        timeInt     = False

        # transversity amplitudes
        A0Mag2Val    = 0.60
        AperpMag2Val = 0.16
        AparMag2Val  = 1. - A0Mag2Val - AperpMag2Val
        A0PhVal    =  0.
        AperpPhVal = -0.17
        AparPhVal  =  2.50
    
        # CP violation parameters
        phiCPVal      = -0.04

        # B lifetime parameters
        GammaVal  = 0.679
        dGammaVal = 0.060
        dMVal     = 17.8
        tResSigma = 0.045

        angleNames = ( 'cos#kern[0.1]{#theta_{K}}', 'cos#kern[0.1]{#theta_{l}}', '#varphi [rad]' )
        effLabels  = (  '#int d_{}cos#theta_{#mu} d#varphi #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
                      , '#int d_{}cos#theta_{K} d#varphi #varepsilon_{#Omega}(#Omega) / (4#pi #LT#varepsilon_{#Omega}#GT)'
                      , '#int d_{}cos#theta_{K} dcos#theta_{#mu} #varepsilon_{#Omega}(#Omega) / (4 #LT#varepsilon_{#Omega}#GT)'
                     )

        ####################################
        ## create variables and read data ##
        ####################################

        # import RooFit wrappers
        from P2VV.Load import RooFitOutput

        # Set global object prefix
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        setParNamePrefix( 'mc' )

        # angular functions
        #from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        #self._angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )
        self._angleFuncs = dataPdfBuilder['angleFuncs']

        # get observables.
        self._obsSet = dataPdfBuilder['obsSetP2VV']        
        for obs in self._obsSet: 
                if obs.GetName() == 'time': 
                    time = obs
                    break
        if self._trueTime: 
            from P2VV.RooFitWrappers import RealVar    
            trueTime = RealVar( 'truetime', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ))
            self._obsSet.remove( time )
            self._obsSet.append( trueTime )
     
        from  P2VV.RooFitWrappers import Category
        iTag = Category('iTag',     Title = 'Initial state flavour tag',    Observable = True, States = { 'Untagged' : 0 } )
        self._obsSet.append( iTag )

        # Read data
        self._readMonteCarloData()
        
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

        tResArgs = { }
        if tResModel == 'Gauss' :
            from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
            tResArgs['time']         = self._obsSettime
            tResArgs['timeResSigma'] = tResSigma
        elif tResModel == '3Gauss' :
            from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            tResArgs['time'] = time
        else :
            from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
            tResArgs['time'] = trueTime
        timeResModel = TimeResolution( **tResArgs )

        # CP violation parameters
        from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
        lambdaCP = CPParam( lambdaCPSq = 1., phiCP = phiCPVal )
   
        # tagging parameters
        from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
        taggingParams = TaggingParams()

        # coefficients for time functions
        from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        index = [ 'A0', 'Apar', 'Aperp', 'AS' ] if self._modelSwave else [ 'A0', 'Apar', 'Aperp' ]
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, amplitudes, lambdaCP, index )
        
        # build underlying physics PDF
        args = dict(    time            = time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime
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
        self._pdf = pdf = BTagDecay(self._pref + '_sig_t_angles_tagCat_iTag', **args )
        
        self._AngAmpsParsVals = dict(A0Mag2Val=A0Mag2Val, 
                                     AperpMag2Val=AperpMag2Val,
                                     AparMag2Val=AparMag2Val,
                                     A0PhVal=A0PhVal,
                                     AperpPhVal=AperpPhVal,
                                     AparPhVal=AparPhVal,
                                     phiCPVal=phiCPVal,
                                     GammaVal=GammaVal,
                                     dGammaVal=dGammaVal,
                                     dMVal=dMVal,
                                     tResSigma=tResSigma,
                                     lambdaCPSq=1)     
        

    def _readMonteCarloData(self):
        # ntuple variables
        from P2VV.RooFitWrappers import RealVar
        B_P        = RealVar( 'B_P',        Title='B_P',          Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        B_PT       = RealVar( 'B_Pt',       Title='B_Pt',         Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        Kplus_P    = RealVar( 'Kplus_P',    Title = 'Kplus_P',    Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        Kplus_PX   = RealVar( 'Kplus_PX',   Title = 'Kplus_PX',   Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        Kplus_PY   = RealVar( 'Kplus_PY',   Title = 'Kplus_PY',   Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        Kplus_PZ   = RealVar( 'Kplus_PZ',   Title = 'Kplus_PZ',   Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        Kminus_P   = RealVar( 'Kminus_P',   Title = 'Kminus_P',   Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        Kminus_PX  = RealVar( 'Kminus_PX',  Title = 'Kminus_PX',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        Kminus_PY  = RealVar( 'Kminus_PY',  Title = 'Kminus_PY',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        Kminus_PZ  = RealVar( 'Kminus_PZ',  Title = 'Kminus_PZ',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muplus_P   = RealVar( 'muplus_P',   Title = 'muplus_P',   Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        muplus_PX  = RealVar( 'muplus_PX',  Title = 'muplus_PX',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muplus_PY  = RealVar( 'muplus_PY',  Title = 'muplus_PY',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muplus_PZ  = RealVar( 'muplus_PZ',  Title = 'muplus_PZ',  Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muminus_P  = RealVar( 'muminus_P',  Title = 'muminus_P',  Unit = 'MeV/c',    Observable = False,  MinMax = ( 0  , 1e7 )   )
        muminus_PX = RealVar( 'muminus_PX', Title = 'muminus_PX', Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muminus_PY = RealVar( 'muminus_PY', Title = 'muminus_PY', Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        muminus_PZ = RealVar( 'muminus_PZ', Title = 'muminus_PZ', Unit = 'MeV/c',    Observable = False,  MinMax = (-1e7, 1e7 )   )
        
        self._ntupleVars = [ Kplus_P, Kplus_PX, Kplus_PY, Kplus_PZ, Kminus_P, Kminus_PX, Kminus_PY, Kminus_PZ,\
                             muminus_P, muminus_PX, muminus_PY, muminus_PZ, muplus_P, muplus_PX, muplus_PY, muplus_PZ, B_P, B_PT ]
        
        bkgcatCut      = '(bkgcat == 0 || bkgcat == 50)'
        trackChiSqCuts = 'muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
        massCuts       = 'mass > 5200. && mass < 5550. && mdau1 > 3030. && mdau1 < 3150. && mdau2 > 990. && mdau2 < 1050.'
        timeCuts       = 'time > 0.3 && time < 14. && sigmat < 0.12'
        tagCuts        = '(tagdecision == 0 || tagdecision == -1 || tagdecision == +1)'

        from P2VV.Utilities.DataHandling import readData
        cuts = bkgcatCut + ' && ' + trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
        cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
        data = readData(  self._nTupleFile, dataSetName = self._nTupleName, NTuple = True, observables = self._obsSet+self._ntupleVars, ntupleCuts = cuts )

        self._initData = data


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
                          ReS        = 0.,
                          ImS        = 0.,
                          dM         = mcPars['dM'],
                          dGamma     = mcPars['dGamma'],
                          Gamma      = mcPars['Gamma'],
                          phiCP      = mcPars['phiCP'],
                          lambdaCPSq = mcPars['lambdaCPSq'] 
                          )

        from ROOT import RooArgSet
        pdfParSet = RooArgSet(p._target_() for p in self._pdf.Parameters())
        for k in self._pdf.Parameters(): pdfParSet.find( k.GetName() ).setVal( pars[k.GetName() ])


    def setDataFitParameters(self, dataPars, KKmassCat=None):                                                              
        #TODO: Accomodate KKmass binning.
        from math import sqrt,sin, cos
        if not KKmassCat:
            AparMag2 = 1. - dataPars['A0Mag2'] - dataPars['AperpMag2']
            ASMag2   = dataPars['f_S'] / (1 - dataPars['f_S'])
            for par in self._pdf.Parameters():
                name = par.GetName()
                if name.__contains__('Re'):
                    if   name.__contains__('Aperp'): par.setVal(  sqrt(dataPars['AperpMag2']/dataPars['A0Mag2']) * cos(dataPars['AperpPhase'])  )
                    elif name.__contains__('Apar'):  par.setVal(  sqrt(     AparMag2 / dataPars['A0Mag2']      ) * cos(dataPars['AparPhase'] )  )
                    elif name.__contains__('AS'):    par.setVal(  sqrt(       ASMag2 / dataPars['A0Mag2']      ) * cos(dataPars['ASOddPhase'])  )
                    elif name.__contains__('A0'):    par.setVal(  cos(dataPars['A0Phase'])                                                      )
                    
                    
                elif name.__contains__('Im'):
                    if   name.__contains__('Aperp'): par.setVal( sqrt(dataPars['AperpMag2']/dataPars['A0Mag2']) * sin(dataPars['AperpPhase'])  )
                    elif name.__contains__('Apar'):  par.setVal( sqrt(     AparMag2/dataPars['A0Mag2']        ) * sin(dataPars['AparPhase'] )  )
                    elif name.__contains__('AS'):    par.setVal( sqrt(       ASMag2 / dataPars['A0Mag2']      ) * sin(dataPars['ASOddPhase'])  )
                    elif name.__contains__('A0'):    par.setVal( sin(dataPars['A0Phase'])                                                      )
                else: par.setVal( dataPars[par] )
        else: print ' ' # 'Impliment KK Mass cat parametrization of the MC pdf.'


    def calculateWeights(self, iterNumb, dataParameters, data):
        self._iterNumb = iterNumb
        self._currentDataSet = data

        from ROOT import RooArgSet
        normVars =  RooArgSet(obs._target_() for obs in self._obsSet)
        
        # Reweights MC verticaly to match the Physics of data.
        nominators, denominators,weights = [], [], []
            
        self._pdf.attachDataSet( data ) # make the pdf dependant directly on data
        
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'    
        self.setMonteCarloParameters()     
        for nev in xrange(data.numEntries()):
            data.get(nev)
            denominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating nominators for phyisics matching weights'
        self.setDataFitParameters(dataParameters) 
        for nev in xrange(data.numEntries()):
            data.get(nev)
            nominators.append( self._pdf.getVal(normVars) )
        
        print 'P2VV - INFO: Calculating phyisics matching weights'
        for n,d in zip(nominators,denominators): weights += [n/d]
        self._allWeights['weightsSet_%s'%iterNumb] =  weights

  
    def combineWeights(self):  
        currentWeights = self._allWeights['weightsSet_%s'%(self._iterNumb)]
        if self._iterNumb==1:
            self._combWeights = self._allWeights['weightsSet_1']
        else:
            assert len(self._combWeights)==len(currentWeights)
            prod = []
            for w1, w2 in zip(self._combWeights,currentWeights): prod.append(w1*w2)
            self._combWeights = prod
        

    def writeWeights(self, weightsName='weightPhys'):
        from ROOT import RooArgSet,RooRealVar,RooDataSet
        from P2VV.RooFitWrappers import RealVar
        physWeightVar = RealVar(weightsName, Title=weightsName, Observable = False, MinMax = (-1e5, 1e7 ) )
        #physWeightVar = RooRealVar( weightsName, weightsName,        -1e5, 1e3           )
        weightsSet    = RooDataSet( 'weightsSet', 'weightsSet', RooArgSet(physWeightVar) )
        
        for weight in self._combWeights:
            physWeightVar.setVal( weight )
            weightsSet.add( RooArgSet(physWeightVar) )
        
        #if self._iterNumb!=1:self._currentDataSet.reduce( RooArgSet(self._obsSet + self._ntupleVars) )
        self._currentDataSet.merge( weightsSet )
        self._currentDataSet.SetName('MC_physicsReweighted_%s_iteration'%self._iterNumb )
        print 'P2VV - INFO: Phyisics matching weights added to dataset: '+'MC_physicsReweighted_%s_iteration'%self._iterNumb

        self._weightsName = weightsName


    def getPdf(self):               return self._pdf
    def getAngleFunctions(self):    return self._angleFuncs
    def getObservables(self):       return self._obsSet
    def getInitialMCafterSel(self): return self._initData
    def getNtupleVars(self):        return self._ntupleVars
    def getWeightName(self):        return self._weightsName
    def getAllWeights(self):        
        self._allWeights['combWeights'] = self._combWeights
        return self._allWeights


# Match MC to sWeighted data with horizontal reweighting of B_P and recalculate angles.
class MatchWeightedDistributions():
    def __init__( self,  **kwargs ):
        print 'P2VV - INFO: Initialised kinematic reweighting class: matchWeightedDistributions().'
        self._inTree      = kwargs.pop('inTree')
        self._outTree     = kwargs.pop('outTree')
        self._nBins = kwargs.pop('nBins', '1000')
        self._vars  = kwargs.pop('whichVars','Kminus_P')
        self._itNum = kwargs.pop('itNum')

        self._inWeightName   = kwargs['inWeightName']
        self._outWeightName  = kwargs['outWeightName']
        self._PhysWeightName = kwargs['PhysWeightName']

        from ROOT import RooDataSet
        if type(self._inTree)  ==RooDataSet:self._inTree  = self._inTree.buildTree()
        if type(self._outTree) ==RooDataSet:self._outTree = self._outTree.buildTree()


    def mimicWeights(self):
        print 'P2VV - INFO: Mimicing weights of variables: ', self._vars
        self._mimicedVars = dict( inDistr={}, outDistr={} )
        for var in self._vars:
            self._mimicedVars['inDistr'][var]  = self.MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins )
            self._mimicedVars['outDistr'][var] = self.MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins )


    def MimicWeightedDistribution( self, t, var, wPref, Nbins=1000 ):
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
                newDistribution.append( rdm.Uniform(sumW[b]['bounds'][0], sumW[b]['bounds'][1]) )

        return newDistribution


    def reweightMC( self, obsSet, copyVars = None):
        self._obsSet = obsSet
        self._copyVars = copyVars
        self._muonSet, self._KaonSet, self._BmomSet = [], [],[]
        for var in self._copyVars:
            if var.GetName().startswith('mu'):
                self._muonSet.append(var)
            elif var.GetName().startswith('K'):
                self._KaonSet.append(var)
            elif var.GetName().startswith('B'):
                self._BmomSet.append(var)
            else: continue
            
        inDistrList  = self._mimicedVars['inDistr'][self._vars[0]]
        outDistrList = self._mimicedVars['outDistr'][self._vars[0]]
        self.TransformAnglesWithMomentumWeight(self._inTree, inDistrList, outDistrList, Nbins=self._nBins)


    def TransformAnglesWithMomentumWeight(self,t, pin, pout, Nbins= None):
        """ t: TTree, pin: original momentum distribution (python list), pout : the momentum distribution you want (python list)
        Nbins controls the number of points for the transformation functions
        """
        if Nbins==None: Nbins=self._nBins

        # Transformation of input and output distributions to uniform.
        Udat = UniFunc(pout, nbinsmax = Nbins)
        Umc = UniFunc(pin, nbinsmax = Nbins)

        # Put the newly recalculated angles plus time and true time in a RooDataSet.
        from ROOT import RooDataSet, RooArgSet, gROOT
        helcosthetaK = self._obsSet[1]
        helcosthetaL = self._obsSet[2]
        helphi       = self._obsSet[3]
        time         = self._obsSet[0]
        recalculatedVars = RooArgSet( self._obsSet[1:4] + self._KaonSet + self._BmomSet ) 
        physWeightsVar   = time.ws().obj(self._PhysWeightName)
        
        gROOT.cd('PyROOT:/') # This is necessary to create later a Ttree object later
        recalculatedData = RooDataSet( 'MomRewMC_%s_Iter'%self._itNum, 'MomRewMC_%s_Iter'%self._itNum, recalculatedVars )
        copiedData       = RooDataSet( 'copiedData', 'copiedData', self._inTree, RooArgSet(self._muonSet + [time,physWeightsVar])      )

        from ROOT import TDatabasePDG
	MeV = 1000 # TDatabasePDG is in GeV, this is the factor needed to go to MeV
        PDG = TDatabasePDG()
        Mmu = PDG.GetParticle('mu-').Mass()*MeV
        Mk  = PDG.GetParticle('K-').Mass()*MeV
        
        from ROOT import TVector3, TLorentzVector
        from math import sqrt
        _VM2LV  = lambda v,m : TLorentzVector( v, sqrt( m*m + v.Mag2() ) ) # TVector3,mass to TLorentzVector 
        _E2V    = lambda entry, label : TVector3( getattr(entry,label+'_PX'),getattr(entry,label+'_PY'),getattr(entry,label+'_PZ')) # entry to TVenctor3
        _B3PMAG = lambda K1,K2,mu1,mu2: (K1+K2+mu1+mu2).Mag() 
        _BPT    = lambda K1,K2,mu1,mu2: sqrt( (K1+K2+mu1+mu2).x()**2 + (K1+K2+mu1+mu2).y()**2 )
        
        def setKaonMomVals( k1 ,k2 ):
            for var in self._KaonSet:
                if var.GetName().__contains__('plus'):
                    if   var.GetName().__contains__('X'):var.setVal(  k1.x()  )
                    elif var.GetName().__contains__('Y'):var.setVal(  k1.y()  )
                    elif var.GetName().__contains__('Z'):var.setVal(  k1.z()  )
                    else:                                var.setVal( k1.Mag() )
                else:
                    if   var.GetName().__contains__('X'):var.setVal(  k2.x()  )
                    elif var.GetName().__contains__('Y'):var.setVal(  k2.y()  )
                    elif var.GetName().__contains__('Z'):var.setVal(  k2.z()  )
                    else:                                var.setVal( k2.Mag() )
        assert(self._BmomSet[0].GetName()=='B_P' and self._BmomSet[1].GetName()=='B_Pt')# Important check
      
        print 'P2VV - INFO: Recalculating decay angles after kinematic distributions matching.'
        for entry in t:
            k1_3P = _E2V(entry,'Kplus')
            k2_3P = _E2V(entry,'Kminus')
            k1_3P.SetMag( Udat.inverse(Umc(k1_3P.Mag())) )
            k2_3P.SetMag( Udat.inverse(Umc(k2_3P.Mag())) )
            
            mu1_3P = _E2V( entry , 'muplus' )
            mu2_3P = _E2V( entry , 'muminus')
           
            cThK, cThL, phi = HelicityAngles( Kpl_P = _VM2LV( k1_3P,  Mk  ),
                                              Kmi_P = _VM2LV( k2_3P,  Mk  ),
                                              mPl_P = _VM2LV( mu1_3P, Mmu ),
                                              mMi_P = _VM2LV( mu2_3P, Mmu )
                                               )

            helcosthetaK.setVal( cThK )
            helcosthetaL.setVal( cThL )
            helphi.setVal( phi )
            self._BmomSet[0].setVal( (k1_3P + k2_3P + mu1_3P + mu2_3P).Mag() )
            self._BmomSet[1].setVal( (k1_3P + k2_3P + mu1_3P + mu2_3P).Pt()  )
            setKaonMomVals( k1_3P, k2_3P )

            recalculatedData.add( recalculatedVars )
     
        gROOT.cd('PyROOT:/')
        recalculatedData.merge(copiedData )
        self._recalculatedData = recalculatedData 
        

    def getDataSet(self, tree=False): 
        if tree: return self._recalculatedData.buildTree()
        else:    return self._recalculatedData



class BuildBs2JpsiKK2011sFit():
    def __init__(self,**kwargs):
        print 'P2VV - INFO: Initialised physics reweighting class: buildBs2JpsiKK2011sFit().'
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
        pdfConfig = PdfConfig()

        pdfConfig['timeEffHistFile']      = kwargs.pop('timeEffHistFile')
        pdfConfig['timeEffHistUBName']    = kwargs.pop('timeEffHistUBName')
        pdfConfig['timeEffHistExclBName'] = kwargs.pop('timeEffHistExclBName')
    
        self._dataSetPath =  kwargs.pop('dataSetPath', None)
        self._dataSetName =  kwargs.pop('dataSetName', None)
        self._dataType    =  kwargs.pop('dataType', 'TTree')
   
        parFileIn  = kwargs.pop( 'parFileIn',  '' )
        parFileOut = kwargs.pop( 'parFileOut', '' )
        
        # fit options
        pdfConfig['fitOptions'] = kwargs.pop( 'fitOpts', None)
        if not pdfConfig['fitOptions']:
            pdfConfig['fitOptions'] = dict(  NumCPU    = 2
                                             , Optimize  = 2
                                             , Minimizer = 'Minuit2'
                                             , Offset    = True
                                             #               , Hesse     = False
                                             , Timer     = True
                                             #               , Verbose   = True
                                             )
        self._FitResults = {} # Save all the fit results

        self._corrSFitErr     = '' #  'sumWeight' 
        randomParVals   = ( ) # ( 1., 12345 )
        self._MinosPars = [#  'AparPhase'
                           #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
                           #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
                          ]

        # PDF options 
        KKmassParam = kwargs.pop( 'KKmassBins', None )
        pdfConfig['timeEffType']          = 'paper2012'
        pdfConfig['anglesEffType']        = '' 
        pdfConfig['KKMassBinBounds']      = [ 990., 1050. ] if  not KKmassParam else [ 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12. ]
        pdfConfig['CSPValues']            = [ 0.498]        if not KKmassParam else [ 0.959, 0.770, 0.824, 0.968 ] 

        KKMassPars = pdfConfig['obsDict']['KKMass']
        pdfConfig['obsDict']['KKMass'] = ( KKMassPars[0], KKMassPars[1], KKMassPars[2]
                                         , 1020., pdfConfig['KKMassBinBounds'][0], pdfConfig['KKMassBinBounds'][-1] )
        
        pdfConfig['constrainTagging']   = 'constrain'
        pdfConfig['timeResType']           = 'eventNoMean'
        pdfConfig['numTimeResBins']        = 40
        pdfConfig['constrainDeltaM'] = 'constrain'
        pdfConfig['lambdaCPParam'] = 'lambPhi'
        
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
        from P2VV.Utilities.DataHandling import readData
        dataSet = readData( filePath = self._dataSetPath, dataSetName = self._dataSetName,  NTuple = False )
        pdfConfig['signalData'] = dataSet
        pdfConfig['readFromWS'] = True

        # build pdf.
        from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
        pdfConfig['parNamePrefix'] =  'data' # give the data parameters a prefix.
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
            
            # data set with weights corrected for background dilution: for phi_s fit only!
        if self._corrSFitErr == 'sumWeight'\
                or ( type(self._corrSFitErr) != str and hasattr( self._corrSFitErr, '__iter__' ) and hasattr( self._corrSFitErr, '__getitem__' ) ) :
            from P2VV.Utilities.DataHandling import correctSWeights
            self._fitData = correctSWeights( dataSet, 'N_cbkgMass_sw'
                                       , 'KKMassCat' if KKmassParam  else ''
                                       , CorrectionFactors = None if self._corrSFitErr == 'sumWeight' else self._corrSFitErr )

        else :
            self._fitData = dataSet

        self._pdfConfig = pdfConfig
        

    def doFit(self,iterNumb, randomParVals=None):

        # get observables and parameters in PDF
        pdfObs  = self._pdf.getObservables(self._fitData)
        pdfPars = self._pdf.getParameters(self._fitData)

        # float/fix values of some parameters
        for CEvenOdds in self._pdfBuild['taggingParams']['CEvenOdds'] :
            if not self._pdfConfig['sameSideTagging'] :
                CEvenOdds.setConstant('avgCEven.*')
                CEvenOdds.setConstant( 'avgCOdd.*', True )
        else :
            for CEvenOdd in CEvenOdds :
                CEvenOdd.setConstant('avgCEven.*')
                CEvenOdd.setConstant( 'avgCOdd.*', True )

        self._pdfBuild['tagCatsOS'].parameter('wTagDelP1OS').setVal(0.)
        self._pdfBuild['tagCatsSS'].parameter('wTagDelP1SS').setVal(0.)
        self._pdfBuild['tagCatsOS'].setConstant('wTagDelP1')
        self._pdfBuild['tagCatsSS'].setConstant('wTagDelP1')

        self._pdfBuild['amplitudes'].setConstant('C_SP')


        if randomParVals :
            # give parameters random offsets
            import random
            print 'Bs2JpsiKK2011Fit: give floating parameters random offsets (scale = %.2f sigma; seed = %s)'\
                % ( randomParVals[0], str(randomParVals[1]) if randomParVals[1] else 'system time' )
            random.seed( randomParVals[1] if randomParVals[1] else None )
            for par in pdfPars :
                if not par.isConstant() : par.setVal( par.getVal() + 2. * ( random.random() - 0.5 ) * randomParVals[0] * par.getError() ) 
   
        print 120 * '='
        print 'Bs2JpsiKK2011Fit: fitting %d events (%s)' % ( self._fitData.numEntries(), 'weighted' if self._fitData.isWeighted() else 'not weighted' )
        
        RooMinPars = [ ]
        if self._MinosPars :
            print 'Bs2JpsiKK2011Fit: running Minos for parameters',
            for parName in self._MinosPars :
                RooMinPars.append( pdfPars.find(parName) )
                print '"%s"' % RooMinPars[-1],
            print

        fitResult = self._pdf.fitTo( self._fitData, SumW2Error = True if self._corrSFitErr == 'matrix' else False
                                     , Minos = RooMinPars, Save = True,  **self._pdfConfig['fitOptions']
                                     )

        # print parameter values
        from P2VV.Imports import parNames, parValues2011 as parValues
        print 'Bs2JpsiKK2011Fit: parameters:'
        fitResult.SetName('sFit_%s_Iteration'%iterNumb)
        fitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
        fitResult.covarianceMatrix().Print()
        fitResult.correlationMatrix().Print()
        self._FitResults['iter_%s'%iterNumb] = fitResult 
    
        print 120 * '=' + '\n'

        if parFileOut :
            # write parameters to file
            self._pdfConfig.getParametersFromPdf( self._pdf, self._fitData )
            self._pdfConfig.writeParametersToFile( filePath = parFileOut )

    def getFitResult(self,number):return self._FitResults['iter_%s'%iterNumb]


    def getDataSet(self):   return self._fitData
    def getPdf(self):       return self._pdf
    def getAngleFuncs(self) : return 0
    def getPdfBuilderObject(self) : return self._pdfBuild
 
