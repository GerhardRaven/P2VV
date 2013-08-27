###########################################################################################################################################
## Utilities.MCReweighting: P2VV utilities for reweighting Monte Carlo data to create desired distributions                              ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   VS, Vasilis Syropoulos, Nikhef, v.syropoulos@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################

def compareWeightedDistributions(tree, sTree, var, **kwargs):
    sVar      =  kwargs.pop('sVar',      None        )
    cut       =  kwargs.pop('cut',       None        )
    sCut      =  kwargs.pop('sCut',      None        )
    weight    =  kwargs.pop('weight',    None        )
    sWeight   =  kwargs.pop('sWeight',   None        )
    rangeX    =  kwargs.pop('rangeX',    None        )
    bins      =  kwargs.pop('bins',      100         )
    assymPlot =  kwargs.pop('assymPlot', False       )
    save      =  kwargs.pop('Save',      [False,'_'] )

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


# Vertical reweighting class of MC to match the physics of sWeighted data.
class matchMCphysics2Data():
    def __init__( self,nTupleFile, nTupleName = 'DecayTree' ):
        ## TODO::Add code that configures the MC_pdf building upon initilisation.
        # i.e mimic the pdfConfig and pdfBuild stracture.
        print 'P2VV - INFO: Initialised physics reweighting class GeneralUtilities.matchMCphysics2Data()'
        self._nTupleFile = nTupleFile
        self._nTupleName = nTupleName

    def buildMonteCarloPdf(self,TIME=True):
        # Build Mc pdf
        from math import pi, sin, cos, sqrt

        # job parameters
        #makePlots   = True
        physPdf     = True
        tResModel   = ''
        trigger     = ''
        timeInt     = False
        self._TIME = TIME

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

        # workspace
        from P2VV.RooFitWrappers import RooObject
        worksp = RooObject( workspace = 'angEff' ).ws()

        # angular functions
        from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
        angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

        # variables in PDF
        from P2VV.RooFitWrappers import RealVar, Category
        time     = RealVar(  'time',     Title = 'Decay time',      Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. ) )
        trueTime = RealVar(  'truetime', Title = 'True decay time', Unit = 'ps', Observable = True, Value = 0.,  MinMax = ( 0.,  20. ) )
        iTag     = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'Untagged' : 0 } )
        angles   = [ angleFuncs.angles['cpsi'], angleFuncs.angles['ctheta'], angleFuncs.angles['phi'] ]

        obsSet = [ trueTime if TIME else time ] + angles

        # read ntuple
        bkgcatCut      = '(bkgcat == 0 || bkgcat == 50)'
        trackChiSqCuts = 'muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
        massCuts       = 'mass > 5200. && mass < 5550. && mdau1 > 3030. && mdau1 < 3150. && mdau2 > 990. && mdau2 < 1050.'
        timeCuts       = 'time > 0.3 && time < 14. && sigmat < 0.12'
        tagCuts        = '(tagdecision == 0 || tagdecision == -1 || tagdecision == +1)'

        from P2VV.Utilities.DataHandling import readData
        cuts = bkgcatCut + ' && ' + trackChiSqCuts + ' && ' + massCuts + ' && ' + timeCuts + ' && ' + tagCuts
        if trigger == 'ExclBiased' :
            cuts  = 'sel == 1 && sel_cleantail==1 && hlt1_excl_biased_dec == 1 && hlt2_biased == 1 && ' + cuts
            data = readData( self._nTupleFile, dataSetName = self._nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )

        elif trigger == 'Unbiased' :
            cuts = 'sel == 1 && sel_cleantail==1 && hlt1_unbiased_dec == 1 && hlt2_biased == 1 && ' + cuts
            data = readData(  self._nTupleFile, dataSetName = self._nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )

        else :
            cuts = 'sel == 1 && sel_cleantail==1 && (hlt1_unbiased_dec == 1 || hlt1_biased == 1) && hlt2_biased == 1 && ' + cuts
            data = readData(  self._nTupleFile, dataSetName = self._nTupleName, NTuple = True, observables = obsSet, ntupleCuts = cuts )

        #####################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################

        if physPdf :
            # transversity amplitudes
            from P2VV.Parameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(  ReApar  = sqrt(AparMag2Val  / A0Mag2Val) * cos(AparPhVal)
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
            from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
            lambdaCP = CPParam( lambdaCPSq = 1., phiCP = phiCPVal )

            # tagging parameters
            from P2VV.Parameterizations.FlavourTagging import Trivial_TaggingParams as TaggingParams
            taggingParams = TaggingParams()

            # coefficients for time functions
            from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
            timeBasisCoefs = TimeBasisCoefs( angleFuncs.functions, amplitudes, lambdaCP, [ 'A0', 'Apar', 'Aperp' ] )

            # build underlying physics PDF
            args = dict(  time            = time if tResModel in [ 'Gauss', '3Gauss' ] else trueTime
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
            self._pdf = pdf  = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )
            self._angleFuncs = angleFuncs
            self._obsSet = obsSet
            self._data = data

            self._helcosthetaK = angles[0]
            self._helcosthetaL = angles[1]
            self._helphi = angles[2]
            self._trueTime = trueTime
            self._time = time

            self._amplitudes = amplitudes
            self._dMVal = dMVal
            self._dGammaVal = dGammaVal
            self._GammaVal = GammaVal
            self._phiCPVal = phiCPVal
            self._lambdaCP = lambdaCP
            self._data = data
            self._cuts = cuts

    def getPdf(self):            return self._pdf
    def getAngleFunctions(self): return self._angleFuncs
    def getObservables(self):    return self._obsSet
    def getInitialMCafterSel(self):  return self._data

    def setMonteCarloParameters(self, pars=None):
        if not pars:
            pars = dict(  ReAperp     = self._amplitudes['Aperp'].Re.getVal()
                          ,ImAperp    = self._amplitudes['Aperp'].Im.getVal()
                          ,ReApar     = self._amplitudes['Apar'].Re.getVal()
                          ,ImApar     = self._amplitudes['Apar'].Im.getVal()
                          ,ReA0       = self._amplitudes['A0'].Re.getVal()
                          ,ImA0       = self._amplitudes['A0'].Im.getVal()
                          ,dM         = self._dMVal
                          ,dGamma     = self._dGammaVal
                          ,Gamma      = self._GammaVal
                          ,phiCP      = self._phiCPVal
                          ,lambdaCPSq = self._lambdaCP._lambdaCPSq.getVal()
                          )
            from ROOT import RooArgSet
            pdfParSet = RooArgSet(p._target_() for p in self._pdf.Parameters())
            for k in self._pdf.Parameters(): pdfParSet.find( k.GetName() ).setVal( pars[k.GetName() ])

    def setDataFitParameters(self, pars):
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
        amps = Amplitudes(AmbiguityParameters=False, ASParameterization='deltaPerp', AparParameterization='phase'
                          ,prefix     = pars['prefix']
                          ,A0Mag2     = pars['A0Mag2']
                          ,A0Phase    = pars['A0Phase']
                          ,AperpMag2  = pars['AperpMag2']
                          ,AperpPhase = pars['AperpPhase']
                          ,AparPhase  = pars['AparPhase']
                          ,ASOddPhase = pars['ASOddPhase']
                          ,C_SP       = pars['C_SP']
                          ,f_S        = pars['f_S']
                          )

        for p in self._pdf.Parameters():
            key = p.GetName()
            if   key.startswith('Re'):p.setVal( amps[ pars['prefix']+key[2:] ].Re.getVal() )
            elif key.startswith('Im'):p.setVal( amps[ pars['prefix']+key[2:] ].Im.getVal() )
            else:                     p.setVal(       pars[key]                            )

    def calculateWeights(self,dataParameters):
        from ROOT import RooArgSet
        normVars =  RooArgSet(obs._target_() for obs in self._obsSet)
        # Reweights MC according to match the Physics of the sFit to data
        nominators, denominators,weights = [], [], []
        print 'P2VV - INFO: Calculating denominators for phyisics matching weights'
        self.setMonteCarloParameters()
        for event in self._data:
            if self._TIME :self._trueTime.setVal    ( event.find('truetime').getVal()     )
            else:          self._time.setVal        ( event.find('time').getVal()         )
            self._helcosthetaK.setVal( event.find('helcosthetaK').getVal() )
            self._helcosthetaL.setVal( event.find('helcosthetaL').getVal() )
            self._helphi.setVal      ( event.find('helphi').getVal()       )
            denominators.append( self._pdf.getVal(normVars) )

            # Set Monte carlo parameters to pdf and catch the pdf value for each event
        print 'P2VV - INFO: Calculating nominators for phyisics matching weight'
        self.setDataFitParameters(dataParameters) # dataParameters dict defined constructMCpdf
        for event in self._data:
            if self._TIME :self._trueTime.setVal    ( event.find('truetime').getVal()     )
            else:          self._time.setVal        ( event.find('time').getVal()         )
            self._helcosthetaK.setVal( event.find('helcosthetaK').getVal() )
            self._helcosthetaL.setVal( event.find('helcosthetaL').getVal() )
            self._helphi.setVal      ( event.find('helphi').getVal()       )
            nominators.append( self._pdf.getVal(normVars) )
        print 'P2VV - INFO: Calculating phyisics matching weights'
        for n,d in zip(nominators,denominators): weights += [n/d]
        self._weights = weights

    def writeWeightsToFile(self,path, weightsName='weightPhys'):
        # TODO:: Use Roels ROOT function that writes weights to file, your way is not optimal
        
        # Fill the MC tree with the new weights column.
        from ROOT import TFile
        initFile = TFile.Open(self._nTupleFile,'READ')
        initTree = initFile.Get(self._nTupleName) 

        # TODO:This is not the correct way to copy big trees fix this
        outFile = TFile.Open(path,'RECREATE')
        outTree = initTree.CopyTree(self._cuts)

        # Check tuple alignment
        try: assert outTree.GetEntries()==self._data.numEntries()== len(self._weights)
        except AssertionError: print 'P2VV - ERROR: Source and target ntuple files are not aligned'

        #Create new branch fro the weights
        from array import array
        address = array('f',[0])
        branch = outTree.Branch( weightsName, address, weightsName + '/F' )

        for w in self._weights:
            address[0] = w
            branch.Fill()

        self._weightedNtuplePath = path
        self._weightedNtupleName = outTree.GetName()
        self._physWeightsName = weightsName

        outFile.cd()
        outTree.Write()
        outFile.Close()
        initFile.Close()
        del outFile
        del initFile
        print 'P2VV - INFO: Phyisics matching weights written to file: ' + path


# Match MC to sWeighted data with horizontal reweighting of B_P and recalculate angles.
class matchWeightedDistributions():
    def __init__( self, outputname,  **kwargs ):
        print 'P2VV - INFO: Initialised kinematic reweighting class GeneralUtilities.matchWeightedDistributions()'
        mcInfo      = kwargs.pop('mcInfo')
        sDInfo      = kwargs.pop('sDInfo')
        self._nBins = kwargs.pop('nBins', '1000')
        self._vars  = kwargs.pop('whichVars')
        self._itNum = kwargs.pop('itNum')

        self._tmc       = mcInfo['path']
        self._tsD       = sDInfo['path']
        self._tmcName   = mcInfo['name']
        self._tsDName   = sDInfo['name']
        self._tmcWeight = mcInfo['weight']
        self._tsDWeight = sDInfo['weight']
        self._outFile   = outputname

        if self._vars=='KaonMomenta':
            self._vars = dict( mc    = { 'vars':['Kminus_P'], 'where':mcInfo }
                              ,sData = { 'vars':['Kminus_P'], 'where':sDInfo}
                               )
        elif self._vars:
            self._vars['mc'].update( {'where':mcInfo} )
            self._vars['sData'].update( {'where':sDInfo} )
        else: print 'P2VV - ERROR: Do not know where to get input variables for calss matchWeightedDistributions'

    def mimicWeights(self):
        # Warning: Mimicinc might increase the stat error on the acceptance determination. Validate this
        from ROOT import TFile
        mimicedVars = dict( mc={}, sData={} )
        print 'P2VV - INFO: Mimicing weights of variables: ', self._vars['mc']['vars']
        for varList in self._vars.keys():
            t = TFile.Open( self._vars[varList]['where']['path'] ).Get( self._vars[varList]['where']['name'] )
            w = self._vars[varList]['where']['weight']
            b = self._nBins
            mimicedVars[varList].update(dict( 
                    (  v , self.MimicWeightedDistribution(t,v,w,b)  )for v in self._vars[varList]['vars']
                ))
        self._mimicedVars = mimicedVars

    def MimicWeightedDistribution(self,t,var,wPref,Nbins=1000):
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

    def TransformAnglesWithMomentumWeight(self,t, pin,pout, outputname = "/tmp/test", Nbins= None):
        """ t: TTree, pin: original momentum distribution (python list), pout : the momentum distribution you want (python list)
        Nbins controls the number of points for the transformation functions
        """

        from ROOT import TDatabasePDG
	MeV = 1000 # TDatabasePDG is in GeV, this is the factor needed to go to MeV
        PDG = TDatabasePDG()
        Mmu = PDG.GetParticle('mu-').Mass()*MeV
        Mk  = PDG.GetParticle('K-').Mass()*MeV

        if Nbins==None: Nbins=self._nBins
        names, labels = [], []
        a = t.GetListOfBranches()

        from SomeUtils.GLBasic import UniFunc
        #/cvmfs/lhcb.cern.ch/lib/lhcb/URANIA/URANIA_v1r1/InstallArea/x86_64-slc6-gcc46-opt/python/SomeUtils
        print 'P2VV - INFO: Matching kinematic distributions.'
        Udat = UniFunc(pout, nbinsmax = Nbins)
        Umc = UniFunc(pin, nbinsmax = Nbins)

        # for branch in a:
        #         name = branch.GetName()
        #         names.append(branch.GetName())
        #         labels.append(branch.GetName() + "/F")
        # num = self._itNum
        # labels += ["Kplus_P_mod%s/F" %num, "Kminus_P_mod%s/F"%num, \
        #            "helcosthetaK_mod%s/F"%num, "helcosthetaL_mod%s/F"%num, "helphi_mod%s/F"%num]
        #from RTuple import RTuple
        #tup = RTuple(outputname, labels)

        from math import sqrt, cos
        from SomeUtils.alyabar import vunit, vector, P_VV_angles, vmod

        # from ROOT import TVector3, TLorentzVector
        # _LV2L  = lambda lv : [ lv[3], [ lv[0],lv[1],lv[2] ]  ]
        # _VM2LV = lambda v,m : TLorentzVector( v, sqrt( m*m + v.Mag2() ) )
        # _E2V   = lambda entry, label : TVector3( getattr(entry,label+'_PX'),getattr(entry,label+'_PY'),getattr(entry,label+'_PZ'))
        # _VM2L  = lambda v,m : _LV2L( _VM2LV(v,m) )

        # Put the newly recalculated angles plus time and true time in a RooDataSet.
        from ROOT import RooDataSet, RooArgSet
        helcosthetaK = self._obsSet[1]
        helcosthetaL = self._obsSet[2]
        helphi = self._obsSet[3]
        time = self._obsSet[0]
        obsSet = RooArgSet(helcosthetaK,helcosthetaL,helphi,time)
        RewData = RooDataSet('MomRewMC_%s_Iter'%self._itNum, 'MomRewMC_%s_Iter'%self._itNum, obsSet)

        print 'P2VV - INFO: Recalculating decay angles after kinematic distributions matching.'

        for entry in t:
            p01 = sqrt(entry.Kplus_PX**2 + entry.Kplus_PY**2 + entry.Kplus_PZ**2)
            p02 = sqrt(entry.Kminus_PX**2 + entry.Kminus_PY**2 + entry.Kminus_PZ**2)

            pmod1 = Udat.inverse(Umc(p01))
            pmod2 = Udat.inverse(Umc(p02))

            #### Modify the momentum scale, not the slop
            p1 = vunit( vector ( entry.Kplus_PX, entry.Kplus_PY, entry.Kplus_PZ))
            p2 = vunit( vector ( entry.Kminus_PX, entry.Kminus_PY, entry.Kminus_PZ))

            p1 = pmod1*p1
            p2 = pmod2*p2

            pmu1 = vector(entry.muplus_PX, entry.muplus_PY, entry.muplus_PZ)
            pmu2 = vector(entry.muminus_PX, entry.muminus_PY, entry.muminus_PZ)

            Ek1= sqrt( Mk**2 + pmod1**2)
            Ek2= sqrt( Mk**2 + pmod2**2)

            Emu1 = sqrt(Mmu**2 + (vmod(pmu1))**2)
            Emu2 = sqrt(Mmu**2 + (vmod(pmu2))**2)

            l0 = [ Ek1, p1]
            l1 = [Ek2, p2]

            l2 = [Emu1, pmu1]
            l3 = [Emu2, pmu2]
            # print l0
            # print l1
            # print l2
            # print l3
            # print  P_VV_angles(l0,l1,l2,l3)

            ### for candidate in t:
            # pmu1 = _E2V( entry , 'muplus' )
            # pmu2 = _E2V( entry , 'muminus')
            # pK1  = _E2V( entry , 'Kplus' )
            # pK2  = _E2V( entry , 'Kminus' )
            # #### Modify the momentum scale, keep the direction
            # pK1.SetMag( Udat.inverse(Umc(pK1.Mag())) )
            # pK2.SetMag( Udat.inverse(Umc(pK2.Mag())) )
            # l0 = _VM2L( pK1, Mk )
            # l1 = _VM2L( pK2, Mk )
            # l2 = _VM2L( pmu1, Mmu )
            # l3 = _VM2L( pmu2, Mmu )
            # print l0
            # print l1
            # print l2
            # print l3
            # print  P_VV_angles(l0,l1,l2,l3)

            Th1,Th2,Phi = P_VV_angles(l0,l1,l2,l3)

            helcosthetaK.setVal(cos(Th1))
            helcosthetaL.setVal(cos(Th2))
            helphi.setVal(Phi)
            if self._obsSet[0].GetName().startswith('true'): time.setVal(entry.truetime )
            else: time.setVal(entry.time )
            RewData.add(obsSet)

            # for name in names: tup.fillItem(name,float(getattr(t,name)))
            # tup.fillItem("Kplus_P_mod%s" %num,pmod1)
            # tup.fillItem("Kminus_P_mod%s" %num,pmod2)
            # tup.fillItem("helcosthetaK_mod%s"%num,cos(Th1))
            # tup.fillItem("helcosthetaL_mod%s"%num,cos(Th2))
            # tup.fillItem("helphi_mod%s"%num,Phi)
            #tup.fill()
        #tup.close()

        from ROOT import TFile
        f = TFile.Open(outputname + '_RDS.root','RECREATE')
        f.cd()
        RewData.Write()
        f.Close()

    def reweightMC(self, outPath=None, var=None, obsSet=None):
        self._obsSet = obsSet
        if not outPath:outPath= self._outFile
        if not var:var='Kminus_P'

        from ROOT import TFile
        tmc = TFile.Open(self._tmc).Get(self._tmcName)
        mcList = self._mimicedVars['mc'][var]
        sDList = self._mimicedVars['sData'][var]

        self.TransformAnglesWithMomentumWeight(tmc,mcList,sDList,outputname=outPath,Nbins=self._nBins)
