###########################################################################################################################################
## Utilities.MCReweighting: P2VV utilities for reweighting Monte Carlo data to create desired distributions                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   VS, Vasilis Syropoulos, Nikhef, v.syropoulos@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################


global _TreeStash # keep reference to trees, so that they are not created multiple times.
_TreeStash = []

global _PlotStash # keep reference to canvases, so that they are not created multiple times.
_PlotStash = []

def CompareWeightedDistributions(tree, sTree, var, **kwargs):
    sVar      = kwargs.pop('sVar',      None      )
    cut       = kwargs.pop('cut',       None      )
    sCut      = kwargs.pop('sCut',      None      )
    weight    = kwargs.pop('weight',    None      )
    sWeight   = kwargs.pop('sWeight',   None      )
    rangeX    = kwargs.pop('rangeX',    None      )
    bins      = kwargs.pop('bins',       50       )
    assymPlot = kwargs.pop('assymPlot', False     )
    save      = kwargs.pop('save',        ''      )
    title     = kwargs.pop('title',       ''      )

    from ROOT import RooDataSet, gROOT
    if type(tree)  == RooDataSet and tree not in _TreeStash:
        gROOT.cd('PyROOT:/')
        tree = tree.buildTree()
        _TreeStash.append(tree)
    if type(sTree) == RooDataSet and sTree not in _TreeStash: 
        gROOT.cd('PyROOT:/')
        sTree = sTree.buildTree()
        _TreeStash.append(sTree)

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
    hm.SetTitle(title)

    hs.SetMarkerStyle(20)
    hs.SetMarkerSize(.5)
    hs.SetTitle(title)
    hs.SetStats()

    def getSumOfWeights(t,pref,cut):
        ## TODO: should return sumW of the selected events by the cut string.
        if cut: print 'WARNING: Returned number is sumW of the entire tree, not of the subseset selected by cut. '
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
            if save:
               c_distr.SaveAs( save )
               asymPlot.SaveAs('assym_' + save)
            return c_distr, asymPlot
    else:
        c_distr.cd()
        hm.Draw()
        hs.Draw('same')
        if save: c_distr.SaveAs( save )
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


# Vertical reweighting class of MC to match the physics of sWeighted data.
class MatchMCphysics2Data():
    def __init__( self,nTupleFile, **kwargs ):
        print 'P2VV - INFO: Initialised physics reweighting class: matchMCphysics2Data().'
        self._nTupleFile = nTupleFile
        self._nTupleName = kwargs.pop('nTupleName', 'DecayTree')
        
        self._trueTime   = kwargs.pop('trueTime',  True  ) 
        self._modelSwave = kwargs.pop('modelSwave', True ) 
        self._selectData       = kwargs.pop('selectData', False)

        self._allWeights = {} # collect all the weights throughout the loop
        self._pref       = 'mc' # prefix for objects created by this instance
        
        # set the  prefix
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

        ####################################
        ## create variables and read data ##
        ####################################

        # import RooFit wrappers
        from P2VV.Load import RooFitOutput

        # get observables.
        self._mcObsSet = [ obs for obs in dataPdfBuilder['obsSetP2VV'] ]        
        for obs in self._mcObsSet: 
                if obs.GetName()   == 'time':        time = obs
                elif obs.GetName() =='helcosthetaK': helcosthetaK = obs
                elif obs.GetName() =='helcosthetaL': helcosthetaL = obs
                elif obs.GetName() =='helphi':       helphi = obs
        if self._trueTime: 
            from P2VV.RooFitWrappers import RealVar    
            trueTime = RealVar( 'truetime', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0., 20. ))
            self._mcObsSet.remove( time )
            self._mcObsSet.append( trueTime )
     
        from  P2VV.RooFitWrappers import Category
        iTag = Category('iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
        self._mcObsSet.append( iTag )

        # angular functions
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        self._angleFuncs = AngleFuncs( cpsi = helcosthetaK, ctheta = helcosthetaL, phi = helphi )
       
        # Read data        
        self._createNtupleVars()
        if self._selectData: 
            self._selectMonteCarloData()
        else: 
            from ROOT import TFile
            self._initData = TFile.Open(self._nTupleFile).Get(self._nTupleName)
        from ROOT import RooDataSet
        assert type(self._initData)==RooDataSet, 'P2VV - ERROR: MC data must be of RooDataSet type.'
        
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
            tResArgs['time']         = time
            tResArgs['timeResSigma'] = tResSigma
        elif tResModel == '3Gauss' :
            from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            tResArgs['time'] = time
        else :
            from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution as TimeResolution
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
        
        self._AngAmpsParsVals = dict(A0Mag2=A0Mag2Val, 
                                     AperpMag2  = AperpMag2Val,
                                     AparMag2   = AparMag2Val,
                                     A0Ph       = A0PhVal,
                                     AperpPh    = AperpPhVal,
                                     AparPh     = AparPhVal,
                                     phiCP      = phiCPVal,
                                     Gamma      = GammaVal,
                                     dGamma     = dGammaVal,
                                     dM         = dMVal,
                                     tResSigma  = tResSigma,
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

        for var in self._pdf.Parameters(): var.setVal( pars[ var.GetName().partition('_')[2] ] ) # remove prefix from the key. 


    def setDataFitParameters(self, dataPars, KKmassCat=None):                                                              
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
                else: par.setVal( dataPars[ name.partition('_')[2] ] ) # remove prefix from the key.
        else: print ' ' # TODO: LoopOver the KKbins and set f_S_i and ASOddPhase_i


    def calculateWeights(self, iterNumb, dataParameters, data):
        self._iterNumb = iterNumb
        self._currentDataSet = data

        from ROOT import RooArgSet
        normVars =  RooArgSet(obs._target_() for obs in self._mcObsSet)
        
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
        self._allWeights['weightsSet_%s'%iterNumb] = weights

  
    def combineWeights(self):  
        currentWeights = self._allWeights['weightsSet_%s'%(self._iterNumb)]
        if self._iterNumb==1:
            self._combWeights = self._allWeights['weightsSet_1']
        else:
            print 'P2VV - INFO: Combining phyisics matching weights'
            assert len(self._combWeights)==len(currentWeights)
            prod = []
            for w1, w2 in zip(self._combWeights,currentWeights): prod.append(w1*w2)
            self._combWeights = prod
        self._allWeights['combWeights'] = self._combWeights


    def writeWeights(self, weightsName='weightPhys'):
        from ROOT import RooArgSet,RooRealVar,RooDataSet,RooNumber
        from P2VV.RooFitWrappers import RealVar
        RooInf  = RooNumber.infinity()
        weightsVar     = RooRealVar( weightsName,  weightsName, -RooInf, RooInf   )
        weightsArgSet  = RooArgSet( weightsVar )
        weightsDataSet = RooDataSet( 'weightsSet', 'weightsSet', weightsArgSet )
        
        for weight in self._combWeights:
            weightsVar.setVal( weight )
            weightsDataSet.add( weightsArgSet )
        
        self._currentDataSet.merge( weightsDataSet )
        self._currentDataSet.SetName('MC_physicsReweighted_%s_iteration'%self._iterNumb )
        print 'P2VV - INFO: Phyisics matching weights added to dataset: '+'MC_physicsReweighted_%s_iteration'%self._iterNumb

        self._weightsName = weightsName
    

    def compareAngles(self, save=True):
        # NOTE: It is not necessary to scale the distributions, scale factor is 1.006 .
        from ROOT import gROOT
        gROOT.cd('PyROOT:/')       
        data    = self._currentDataSet.buildTree()        
  
        time = 'truetime' if self._trueTime else 'time'
        observables = [ 'helcosthetaK', 'helcosthetaL', 'helphi', time]

        from ROOT import TCanvas
        canvs = TCanvas('anglesPhysReweight_%s'%self._iterNumb, 'anglesPhysReweight_%s'%self._iterNumb)
        canvs.Divide(2,2)
        for i, obs in zip(range(1,5), observables):
            canvs.cd(i)
            data.Draw( obs )
            data.Draw( obs , self._weightsName, 'same err' )
        _PlotStash.append(canvs)
        
        if save: canvs.Print('anglesPhysReweight_%s.pdf'%self._iterNumb)

        
    def compareTrackMomenta(self, save=True):
        # NOTE: It is not necessary to scale the distributions, scale factor is 1.006 .
        from ROOT import gROOT
        gROOT.cd('PyROOT:/')        
        data = self._currentDataSet.buildTree() 
        
        kaons = ['Kplus_PX',  'Kplus_PY',  'Kplus_PZ',  'Kminus_PX',  'Kminus_PY',  'Kminus_PZ',  'Kplus_P',  'Kminus_P' ]
        muons = ['muplus_PX', 'muplus_PY', 'muplus_PZ', 'muminus_PX', 'muminus_PY', 'muminus_PZ', 'muplus_P', 'muminus_P' ]

        from ROOT import TCanvas
        canvKaon = TCanvas('KaonMomentaPhysReweight_%s'%self._iterNumb, 'KaonMomentaPhysReweight_%s'%self._iterNumb)
        canvmuon = TCanvas('muonMomentaPhysReweight_%s'%self._iterNumb, 'muonMomentaPhysReweight_%s'%self._iterNumb)
        canvmuon.Divide(4,2)
        canvKaon.Divide(4,2)
        for i, K, mu in zip(range(1,9), kaons, muons):
            canvKaon.cd(i)
            data.Draw( K )
            data.Draw( K , self._weightsName, 'same err' )
        
            canvmuon.cd(i)
            data.Draw( mu )
            data.Draw( mu , self._weightsName, 'same err' )
        
        _PlotStash.append(canvKaon)
        _PlotStash.append(canvmuon)

        if save: 
            canvKaon.Print('KaonMomentaPhysReweight_%s.pdf'%self._iterNumb)
            canvmuon.Print('muonMomentaPhysReweight_%s.pdf'%self._iterNumb)


    def getPdf(self):               return self._pdf
    def getNtupleVars(self):        return self._ntupleVars
    def getAngleFunctions(self):    return self._angleFuncs
    def getInitialMCafterSel(self): return self._initData
    def getAllWeights(self):        return self._allWeights
    def getMcObsSet(self)  :
        return  [o for o in self._mcObsSet if o.GetName().startswith('hel') or o.GetName().__contains__('time')]



# Match MC to sWeighted data with horizontal reweighting of B_P and recalculate angles.
class MatchWeightedDistributions():
    def __init__( self,  **kwargs ):
        print 'P2VV - INFO: Initialised kinematic reweighting class: matchWeightedDistributions().'
        self._inTree          = kwargs.pop('inTree')
        self._outTree         = kwargs.pop('outTree')
        self._inWeightName    = kwargs.pop('inWeightName')
        self._outWeightName   = kwargs.pop('outWeightName')
        self._mcObsSet        = kwargs.pop('observables')
        self._copyVars        = kwargs.pop('spectatorVars', None)
        self._physWeightsVar  = self._inTree.get().find( self._inWeightName )
        
        self._nBins = kwargs.pop('nBins', '1000')
        self._vars  = kwargs.pop('whichVars','Kminus_P')
        self._itNum = kwargs.pop('itNum')
       
        from ROOT import RooDataSet, gROOT
        if type(self._inTree)  == RooDataSet:
            gROOT.cd('PyROOT:/')
            self._inTree  = self._inTree.buildTree()
        if type(self._outTree) == RooDataSet:
            gROOT.cd('PyROOT:/')
            self._outTree = self._outTree.buildTree()

        def getSumOfWeights(t,pref):
            sumW=0
            for e in t:sumW+=getattr(e,pref)
            return sumW

        self._dataSetsScale = getSumOfWeights(self._outTree, self._outWeightName) / getSumOfWeights(self._inTree,self._inWeightName)

        self._physRewScale = getSumOfWeights(self._inTree,self._inWeightName) / self._inTree.GetEntries()
        assert self._physRewScale<1.01, 'P2VV - WARNING: Scale between original MC and physics reweighted MC is %s, need to propagate the errors in the momentum reweighting.' %self._physRewScale 
        
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
            self._mimicedVars['inDistr'][var]  = self.MimicWeightedDistribution( self._inTree,  var, self._inWeightName , self._nBins )

            if self._mimicedVars.has_key('outDistr'):
                if not self._mimicedVars['outDistr'].has_key(var):
                    print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
                    self._mimicedVars['outDistr'][var] = self.MimicWeightedDistribution( self._outTree, var, self._outWeightName, self._nBins ) 

            else:
                print 'P2VV - INFO: Output Distribution: Mimicing weights of variable: ', var
                self._mimicedVars['outDistr'] = {}
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


    def reweight( self ):
        # mimic the weights 
        self.mimicWeights()

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

        # put the newly recalculated angles plus time in a RooDataSet.
        from ROOT import RooDataSet, RooArgSet
        for obs in self._mcObsSet:
            if   obs.GetName()=='helcosthetaK': helcosthetaK = obs
            elif obs.GetName()=='helcosthetaL': helcosthetaL = obs
            elif obs.GetName()=='helphi':          helphi    = obs
            elif obs.GetName()=='time':             time     = obs
            elif obs.GetName()=='truetime':         time     = obs

        recalculatedVars = RooArgSet( [helcosthetaK,helcosthetaL,helphi] + self._KaonSet + self._BmomSet ) 
        recalculatedData = RooDataSet( 'MomRewMC_%s_Iter'%self._itNum, 'MomRewMC_%s_Iter'%self._itNum, recalculatedVars )
        copiedData       = RooDataSet( 'copiedData', 'copiedData', self._inTree, RooArgSet(self._muonSet + [time, self._physWeightsVar]) )

        from ROOT import TDatabasePDG
	MeV = 1000 # TDatabasePDG is in GeV, this is the factor needed to go to MeV.
        PDG = TDatabasePDG()
        Mmu = PDG.GetParticle('mu-').Mass()*MeV
        Mk  = PDG.GetParticle('K-').Mass()*MeV
        
        from ROOT import TVector3, TLorentzVector
        from math import sqrt
        _VM2LV  = lambda v,m : TLorentzVector( v, sqrt( m*m + v.Mag2() ) ) # TVector3,mass to TLorentzVector 
        _E2V    = lambda entry, label : TVector3( getattr(entry,label+'_PX'),getattr(entry,label+'_PY'),getattr(entry,label+'_PZ')) # entry to TVenctor3
        _B3PMAG = lambda K1,K2,mu1,mu2: (K1+K2+mu1+mu2).Mag() # B momentum.
        _BPT    = lambda K1,K2,mu1,mu2: sqrt( (K1+K2+mu1+mu2).x()**2 + (K1+K2+mu1+mu2).y()**2 ) # B transverse momentum.
       
        # get Kaon variables references. 
        for var in self._KaonSet: 
            if var.GetName().__contains__('plus'):
                if   var.GetName().__contains__('X'): KPX = var
                elif var.GetName().__contains__('Y'): KPY = var
                elif var.GetName().__contains__('Z'): KPZ = var
                elif var.GetName().endswith('P'):     KP  = var
            elif var.GetName().__contains__('minus'):
                if   var.GetName().__contains__('X'): KmX = var
                elif var.GetName().__contains__('Y'): KmY = var
                elif var.GetName().__contains__('Z'): KmZ = var
                elif var.GetName().endswith('P'):     Km  = var

        assert(self._BmomSet[0].GetName()=='B_P' and self._BmomSet[1].GetName()=='B_Pt')# important check
        BP, BPt = self._BmomSet[0], self._BmomSet[1]

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

            # fill the new dataset.
            helcosthetaK.setVal( cThK )
            helcosthetaL.setVal( cThL )
            helphi.setVal( phi )
            
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
        self._recalculatedData = recalculatedData 


    def compareAngles(self, save=True):
        from ROOT import gROOT
        gROOT.cd('PyROOT:/') 

        physRewData = self._inTree
        momRewData  = self._recalculatedData.buildTree()
        sWeightData = self._outTree

        for o in self._mcObsSet:
            if   o.GetName()=='truetime': observables = [ 'helcosthetaK', 'helcosthetaL', 'helphi', 'truetime'] 
            elif o.GetName()=='time':     observables = [ 'helcosthetaK', 'helcosthetaL', 'helphi', 'time'    ] 
        
        from ROOT import TCanvas, gPad
        canvsMomRewData = TCanvas( 'anglesMomReweight_%s'%self._itNum, 'anglesMomReweight_%s'%self._itNum )
        canvsRewVSData  = TCanvas( 'anglesRewVSsData_%s'%self._itNum,  'anglesRewVSsData_%s'%self._itNum  )
        canvsMomRewData.Divide(2,2)
        canvsRewVSData.Divide(2,2)

        for i, obs in zip(range(1,5), observables):
            canvsMomRewData.cd(i)
            physRewData.Draw( obs, self._inWeightName  )
            momRewData.Draw( obs , self._inWeightName, 'same err' )
            
            canvsRewVSData.cd(i)
            momRewData.Draw(  obs + '>>hMomAngles()', self._inWeightName )
            if obs.endswith('time'): sWeightData.Draw( 'time', self._outWeightName, 'same err' )
            else:                    sWeightData.Draw( obs,    self._outWeightName, 'same err' )
            gPad.GetPrimitive('hMomAngles').Scale(self._dataSetsScale)
            gPad.GetPrimitive('hMomAngles').SetName( obs + 'MomRew_%s' %self._itNum)

        _PlotStash.append(canvsMomRewData)
        _PlotStash.append(canvsRewVSData)

        print self._dataSetsScale
        if save: 
            canvsMomRewData.Print('anglesMomReweight_%s.pdf'%self._itNum)
            canvsMomRewData.Print('anglesRewVSsData_%s.pdf'%self._itNum )

        
    def compareTrackMomenta(self, save=True):
        from ROOT import gROOT
        gROOT.cd('PyROOT:/') 

        physRewData = self._inTree
        momRewData  = self._recalculatedData.buildTree()
        sWeightData = self._outTree

        kaons = ['Kplus_PX',  'Kplus_PY',  'Kplus_PZ',  'Kminus_PX',  'Kminus_PY',  'Kminus_PZ',  'Kplus_P',  'Kminus_P' ]
        muons = ['muplus_PX', 'muplus_PY', 'muplus_PZ', 'muminus_PX', 'muminus_PY', 'muminus_PZ', 'muplus_P', 'muminus_P' ]

        from ROOT import TCanvas, gPad
        canvKaon          = TCanvas('KaonMomentaMomReweight_%s'%self._itNum, 'KaonMomentaMomReweight_%s'%self._itNum)
        canvmuon          = TCanvas('muonMomentaMomReweight_%s'%self._itNum, 'muonMomentaMomReweight_%s'%self._itNum)
        canvKaonRewVsData = TCanvas('KaonMomentaRewVSsData_%s'%self._itNum, 'KaonMomentaRewVSsData_%s'%self._itNum)
        canvmuonRewVsData = TCanvas('muonMomentaRewVSsData_%s'%self._itNum, 'muonMomentaRewVSsData_%s'%self._itNum)
        
        for can in [canvKaon, canvmuon, canvKaonRewVsData, canvmuonRewVsData]: can.Divide(4,2)
        for i, K, mu in zip(range(1,9), kaons, muons):
            canvKaon.cd(i)
            physRewData.Draw( K, self._inWeightName             )
            momRewData.Draw(  K, self._inWeightName, 'same err' )
            
            canvKaonRewVsData.cd(i)
            sWeightData.Draw( K + '>>hPhysKaon()', self._outWeightName            )
            momRewData.Draw(  K + '>>hMomKaon()' , self._inWeightName, 'same err' )
            gPad.GetPrimitive('hMomKaon').Scale(self._dataSetsScale)
            gPad.GetPrimitive('hMomKaon').SetName( K + 'MomRew_%s' %self._itNum)

            canvmuon.cd(i)
            physRewData.Draw( mu, self._inWeightName              )
            momRewData.Draw(  mu, self._inWeightName , 'same err' )

            canvmuonRewVsData.cd(i)
            sWeightData.Draw( mu + '>>hPhysMuon()', self._outWeightName            )
            momRewData.Draw(  mu + '>>hMomMuon()' , self._inWeightName, 'same err' )
            gPad.GetPrimitive('hMomMuon').Scale(self._dataSetsScale)
            gPad.GetPrimitive('hMomMuon').SetName( mu + 'MomRew_%s' %self._itNum)

        _PlotStash.append(canvKaonRewVsData)
        _PlotStash.append(canvKaon)
        _PlotStash.append(canvmuonRewVsData)
        _PlotStash.append(canvmuon)

        if save: 
            canvKaon.Print('KaonMomentaMomReweight_%s.pdf'%self._itNum)
            canvmuon.Print('muonMomentaMomReweight_%s.pdf'%self._itNum)
            canvKaonRewVsData.Print('KaonMomentaRewVSsData_%s.pdf'%self._itNum)
            canvmuonRewVsData.Print('muonMomentaRewVSsData_%s.pdf'%self._itNum) 
        

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

        pdfConfig['timeEffHistFiles'] = dict(  file   = kwargs.pop('timeEffHistFile')
                                             , hlt1UB = kwargs.pop('timeEffHistUBName')
                                             , hlt1ExclB = kwargs.pop('timeEffHistExclBName')
                                            )
    
        self._dataSetPath = kwargs.pop('dataSetPath', None)
        self._dataSetName = kwargs.pop('dataSetName', None)
        self._weightsName = kwargs.pop('weightsName', 'JpsiKK_sigSWeight')
           
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
        self._FitResults = {} # collect all the fit results
        self._Moments    = {} # collect all moments
        self._PDFS       = {} # collect allpdfs

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
        pdfConfig['KKMassBinBounds']      = [ 990., 1050. ] if not KKmassParam else [990., 1020. - 12., 1020., 1020. + 12., 1050.]
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
        

    def _setWeightsToDataset(self):
        # reconstruct the RooDataSet object to set the weights. 
        if not self._fitData.isWeighted():
            from ROOT import RooDataSet, RooFit
            self._fitData = RooDataSet( self._dataSetName, self._dataSetName, self._fitData.get(), \
                                            RooFit.Import(self._fitData),  RooFit.WeightVar(self._weightsName)
                                        )


    def doFit(self,iterNumb, pdf=None, randomParVals=None):
        pref = self._pdfConfig['parNamePrefix'] + '_'
        
        if not pdf: self._currentPdf = self._pdf
        else: self._currentPdf = pdf
        self._PDFS['%s_iteration'%iterNumb] = self._currentPdf
    
        # set weights to the RooDataSet object
        self._setWeightsToDataset()         

        # get observables and parameters in PDF
        pdfObs  = self._currentPdf.getObservables(self._fitData)
        pdfPars = self._currentPdf.getParameters(self._fitData)

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

        fitResult = self._currentPdf.fitTo( self._fitData, SumW2Error = True if self._corrSFitErr == 'matrix' else False
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

        
    def multiplyPdfWithAcc( self, effFile, iterNumb=None ):
        # read moments file and multiply pure pdf with angular acceptance
        from P2VV.Utilities.DataMoments import angularMomentIndices, RealMomentsBuilder
        moments = RealMomentsBuilder()
        moments.appendPYList( self._pdfBuild['angleFuncs'].angles, angularMomentIndices( 'weights', self._pdfBuild['angleFuncs'] ) )
        moments.read( effFile )    
        self._Moments['%s_iteration'%iterNumb] = moments # collect all moments  
        if iterNumb: return moments.multiplyPDFWithEff( self._pdf, CoefName = 'effC%d' % iterNumb )
        else:        return moments.multiplyPDFWithEff( self._pdf, CoefName = 'effnom'            )
        

    def plotPdfs( self, Save=True ):
        self._setWeightsToDataset() # in case it is not set to weighteed.         
        from ROOT import TCanvas
        from P2VV.Utilities.Plotting import plot

        # rename pdfs and give a line color to each one of them.
        addPdfsOpts = []
        for key, name, color in zip(  self._PDFS.keys(),
                                     ['addPDF%s'%num for num in range( len(self._PDFS.keys()) ) ],     
                                     [     col       for col in range( len(self._PDFS.keys()) ) ] 
                                     ):
            self._PDFS[key].SetName(name)
            addPdfsOpts.append( dict(LineColor=color +1) )

        # superimporse all the pdfs 
        ObjNam = lambda obj: obj.GetName() + '_AllPdfs'
        canvs = [ TCanvas(ObjNam(o),ObjNam(o)) for o in self.getPdfBuilderObject()['obsSetP2VV'] ]
        for can, obs in zip( canvs, self.getPdfBuilderObject()['obsSetP2VV'] ):
                      plot(can, obs, data        = self._fitData, 
                                     pdf         = self._pdf.getPdf('bin0'), 
                                     addPDFs     = [self._PDFS[key] for key in self._PDFS.keys()], 
                                     addPDFsOpts = addPdfsOpts
                           )
                      if Save==True: can.Print( ObjNam(obs) + '.pdf')

          
    def makeEfficiencyPlots(self): # INCOMPLETE - Does not work
        basisMoments = RealMomentsBuilder()
        basisMoments.appendPYList( angleFuncs.angles, indices, PDF = pdf, IntSet = [ ], NormSet = angles )
        basisMoments.compute(data)

    # add column with efficiency weights to data set

        data = 7
        from ROOT import RooArgSet, RooArgList, RooFormulaVar, RooConstVar
        rooIntSet  = RooArgSet()
        rooNormSet = RooArgSet( var._var for var in angles )
        pdfInt = pdf.createIntegral( rooIntSet, rooNormSet )
        nEvents = RooConstVar( 'nEvents', 'Number of events', data.sumEntries() )
        effWeightList = RooArgList( nEvents, pdfInt )
        effWeight = RooFormulaVar( 'effWeight', 'Efficiency weight', '1./@0/@1', effWeightList )
        data.addColumn(effWeight)
        
    # apply efficiency weights to events in data set
        from ROOT import RooDataSet
        weightedData = RooDataSet( 'effWeightData', 'effWeightData', data.get(), Import = data, WeightVar = ( 'effWeight', True ) )
        
        effTerms = basisMoments.buildPDFTerms()
        effFunc = effTerms.buildAddition('effFunc')
        
        from ROOT import RooArgSet
        ctkSet = RooArgSet( angles[1], angles[2] )
        ctlSet = RooArgSet( angles[0], angles[2] )
        phiSet = RooArgSet( angles[0], angles[1] )
        effInts = (  effFunc.createIntegral( ctkSet, RooArgSet() )
                     , effFunc.createIntegral( ctlSet, RooArgSet() )
                     , effFunc.createIntegral( phiSet, RooArgSet() )
                     )
      

    # plot (binned) efficiency function integrals
        canvs += [  TCanvas( 'cpsiDataIntCanv',   'Angular Efficiency' )
                    , TCanvas( 'cthetaDataIntCanv', 'Angular Efficiency' )
                    , TCanvas( 'phiDataIntCanv',    'Angular Efficiency' )
                    ]
        for ( pad, func, angle, xTitle, yTitle, yScale, norm, nBins )\
                in zip(  canvs[ : 3 ]
                         , effInts
                         , angles
                         , angleNames
                         , effLabels
                         , [ ( 0.88, 1.12 ), ( 0.9328, 1.1872 ), ( 0.88, 1.12 ) ]
                         , [ 1. / 4. / pi, 1. / 4. / pi, 1. / 4. ]
                         , numEffBins
                         ) :
                pad.SetLeftMargin(0.28)
                plot(  pad, angle, weightedData, func
                       , xTitle       = xTitle
                       , yTitle       = yTitle
                       , yScale       = yScale
                       , yTitleOffset = 1.0
                       , frameOpts    = dict( Title = angle.GetTitle(), Bins = nBins )
                       , dataOpts     = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.9, LineWidth = 3
                                              , Rescale = norm * float(nBins) / ( angle.getMax() - angle.getMin() ) )
                       , pdfOpts      = dict( LineColor = kBlue, LineWidth = 4, Normalization = norm )
                       )
                
    
    def printAllMoments( self, nominamMomFile=None ):
        if nominalMomFile:
            from P2VV.Utilities.DataMoments import angularMomentIndices, RealMomentsBuilder
            moments = RealMomentsBuilder()
            moments.appendPYList( self._pdfBuild['angleFuncs'].angles, angularMomentIndices( 'weights', self._pdfBuild['angleFuncs'] ) )
            moments.read( nominamMomFile )
            self._Moments['nominal'] = moments
            
        # ???? scaleFactor = 1 / 16. / sqrt(pi)
        for mom in self._Moments.keys(): self._Moments[mom].Print()

    def printAllFitResults(self):
        for result in self._FitResults.keys():
            self._FitResults[result].Print() # Try Jeroen's printSpecial.
          

    def setFitParameters(self, iterNumb, oldPars):
        fitResult = self._FitResults['iter_%s'%iterNumb]
        parkeys = oldPars.keys()
        parkeys.remove('A0Phase')
        for par in parkeys:
            oldPars[par] = fitResult.floatParsFinal().find( self._pdfConfig['parNamePrefix'] + '_' + par).getVal()

    def getDataSet(self):           return self._fitData
    def getPdf(self):               return self._pdf
    def getPdfBuilderObject(self) : return self._pdfBuild
    def getFitResults(self):        return self._FitResults
    def getMomentsDict(self):       return self._Moments['%s_iteration'%itNum]
    def getObservables(self,which=None):
        if   which=='angles': return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName().startswith('hel') ]
        elif which=='time':   return [o for o in self._pdfBuild['obsSetP2VV'] if o.GetName()==('time')         ]
        else:                 return             self._pdfBuild['obsSetP2VV']


 
