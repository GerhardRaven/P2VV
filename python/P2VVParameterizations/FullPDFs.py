###########################################################################################################################################
## P2VVParameterizations.FullPDFs: Parameterizations of complete PDFs that are used in an analysis                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

class PdfConfiguration( dict ) :
    def __init__( self, parameters = None ) :
        self._parameters = { }
        if parameters != None : self.addParameters(parameters)
        self['parameters'] = self._parameters

    def __getitem__( self, key ) :
        if key not in self and key in self._parameters : return self._parameters[key]
        return dict.__getitem__( self, key )

    def parameters( self ) : return self._parameters

    def addParameter( self, name, values ) :
        if hasattr( values, '__iter__' ) and len(values) == 5 :
            self._parameters[name] = tuple( [ val for val in values ] )
        else :
            raise KeyError('PdfConfiguration.addParameter(): format of "values" argument should be ( value, error, min, max, floating? )')

    def addParameters( self, parameters ) :
        if type(parameters) == dict :
            for name, vals in parameters.iteritems() : self.addParameter( name, vals )
        else :
            raise KeyError('PdfConfiguration.addParameters(): argument "parameters" should be a dictionary')


class Bs2Jpsiphi_Winter2012( PdfConfiguration ) :
    def __init__( self ) :
        # job parameters
        self['makePlots']  = True
        self['SFit']       = False
        self['blind']      = False
        self['nominalPdf'] = True

        self['numEvents'] = 30000

        self['nTupleName'] = 'DecayTree'
        self['nTupleFile'] = 'Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'

        self['timeEffHistFile'] = 'BuBdBdJPsiKsBsLambdab0Hlt2DiMuonDetachedJPsiAcceptance_sPlot_20110120.root'
        self['timeEffHistName'] = 'BsHlt2DiMuonDetachedJPsiAcceptance_Data_Reweighted_sPlot_40bins'

        self['angEffMomentsFile'] = 'effMoments'

        # fit options
        self['fitOptions'] = dict( NumCPU = 1, Timer = 1, Save = True )

        # PDF parameters
        self['tagCats'] = [ ]

        # PDF options
        self['components']         = 'all'              # 'all' / 'signal' / 'background'
        self['transversityAngles'] = False
        self['bkgAnglePdf']        = 'histPdf'
        self['bkgTaggingPdf']      = 'histPdf'
        self['multiplyByTimeEff']  = ''                 # 'all' / 'signal'

        self['taggingConditionals'] = ''                # 'all' / 'tagCat' / 'estWTag' / 'iTag'
        self['numEstWTagBins']      = 100

        self['iTagZeroTrick'] = False
        self['iTagStates'] = { }                        # { } / { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }

        self['eventTimeResolution'] = True
        self['numTimeResBins']      = 100

        self['signalFraction'] = 0.67
        self['massRangeBackground'] = False

        self['amplitudeParam'] = 'phasesSWaveFrac'
        self['polarSWave']     = False

        self['carthLambdaCP'] = False

        self['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                              , ( 'trcostheta', 'cos(#theta_{tr})' )
                              , ( 'trphi',      '#phi_{tr}'        )
                             )

        self['numBMassBins'] = [ 50, 10, 10 ]
        self['numTimeBins']  = 30
        self['numAngleBins'] = ( 5, 7, 9 )

        # initialize PdfConfiguration object
        PdfConfiguration.__init__( self )


class PdfBuilder ( object ) :
    def __init__( self, pdf, observables, parameters ) :
        self._pdf         = pdf
        self._observables = observables
        self._parameters  = parameters

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf(self)         : return self._pdf
    def observables(self) : return self._observables
    def parameters(self)  : return self._parameters


class Bs2Jpsiphi_PdfBuilder ( PdfBuilder ) :
    """builds the PDF for the measurement of phi_s in B_s -> J/psi(->mu^+ mu^-) phi(->K^+ K^-)
    """

    def __init__( self, **kwargs ) :
        ###################################################################################################################################
        ## get/set parameters ##
        ########################

        # copy configuration arguments
        pdfConfig = kwargs.copy()

        # job parameters
        SFit       = pdfConfig.pop('SFit')
        blind      = pdfConfig.pop('blind')
        nominalPdf = pdfConfig.pop('nominalPdf')
        makePlots  = pdfConfig.pop('makePlots')

        numEvents = pdfConfig.pop('numEvents')

        angEffMomentsFile = pdfConfig.pop('angEffMomentsFile')

        nTupleName = pdfConfig.pop('nTupleName')
        nTupleFile = pdfConfig.pop('nTupleFile')

        timeEffHistFile = pdfConfig.pop('timeEffHistFile')
        timeEffHistName = pdfConfig.pop('timeEffHistName')

        fitOpts = pdfConfig.pop('fitOptions')

        angleNames   = pdfConfig.pop('angleNames')
        numTimeBins  = pdfConfig.pop('numTimeBins')
        numAngleBins = pdfConfig.pop('numAngleBins')

        # PDF parameters
        parameters = pdfConfig.pop('parameters')
        tagCats    = pdfConfig.pop('tagCats')

        # PDF options
        components        = pdfConfig.pop('components')
        transAngles       = pdfConfig.pop('transversityAngles')
        bkgAnglePdf       = pdfConfig.pop('bkgAnglePdf')
        bkgTaggingPdf     = pdfConfig.pop('bkgTaggingPdf')
        multiplyByTimeEff = pdfConfig.pop('multiplyByTimeEff')
        numBMassBins      = pdfConfig.pop('numBMassBins')

        self._iTagZeroTrick = pdfConfig.pop('iTagZeroTrick')
        iTagStates = pdfConfig.pop('iTagStates')
        if not iTagStates : iTagStates = { 'B' : +1, 'Bbar' : -1 }
        if +1 not in iTagStates.values() or -1 not in iTagStates.values() or 0 in iTagStates.values() : self._iTagZeroTrick = True
        if self._iTagZeroTrick : iTagStatesDecision = iTagStates
        else                   : iTagStatesDecision = { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }

        tagConds = pdfConfig.pop('taggingConditionals')
        numEstWTagBins = pdfConfig.pop('numEstWTagBins')

        eventTimeRes = pdfConfig.pop('eventTimeResolution')
        numTimeResBins = pdfConfig.pop('numTimeResBins')

        sigFrac = pdfConfig.pop('signalFraction')
        massRangeBackground = pdfConfig.pop('massRangeBackground')

        amplitudeParam = pdfConfig.pop('amplitudeParam')
        polarSWave     = pdfConfig.pop('polarSWave')

        carthLambdaCP = pdfConfig.pop('carthLambdaCP')

        if makePlots :
            # import plotting tools
            from P2VVLoad import ROOTStyle
            from P2VVGeneralUtils import plot
            from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed


        ###################################################################################################################################
        ## create variables (except for tagging category) ##
        ####################################################

        # RooObject wrappers
        from RooFitWrappers import RooObject, RealVar, Category
        ws = RooObject().ws()

        # angular functions
        if nominalPdf or transAngles :
            from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
            self._angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
        else :
            from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
            self._angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

        # variables in PDF (except for tagging category)
        time = RealVar( 'time', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. )
                       , Ranges = dict( Bulk = ( None, 5. ) ), nBins = numTimeBins )
        timeRes = RealVar(  'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, Value = 0.007, MinMax = (0.007, 0.12)
                          , nBins = numTimeResBins )
        timeRes.setBins( numTimeResBins, 'cache' )

        cpsi   = self._angleFuncs.angles['cpsi']
        ctheta = self._angleFuncs.angles['ctheta']
        phi    = self._angleFuncs.angles['phi']

        if nominalPdf or not self._iTagZeroTrick :
            iTag = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = iTagStates )
        estWTag = RealVar( 'tagomega_os', Title = 'Estimated wrong tag probability', Observable = True
                          , Value = 0.25, MinMax = ( 0., 0.50001 ), nBins = numEstWTagBins )
        estWTag.setBins( numEstWTagBins, 'cache' )

        BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True
                        , Value = 5368., MinMax = ( 5200., 5550. ), nBins = numBMassBins[0] + numBMassBins[1] + numBMassBins[2]
                        ,  Ranges = dict(  LeftSideBand  = ( None,  5330. )
                                         , Signal        = ( 5330., 5410. )
                                         , RightSideBand = ( 5410., None  )
                                        )
                       )

        angles = [ cpsi, ctheta, phi ]
        obsSetP2VV = [ time ] + angles
        if not nominalPdf and not ( tagConds in [ 'all', 'iTag' ] or self._iTagZeroTrick ) : obsSetP2VV.append(iTag)
        if not SFit : obsSetP2VV.append(BMass)

        # ntuple variables
        mpsi = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. ), nBins =  32 )
        mphi = RealVar( 'mdau2', Title = 'M(KK)',     Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. ), nBins =  16 )

        tagDecision = Category( 'tagdecision_os', Title = 'Tag decision', Observable = True, States = iTagStatesDecision )
        tagCat = Category( 'tagcat_os',   Title = 'Tagging Category', Observable = True
                          , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                         )
        if not nominalPdf and not tagConds in [ 'all', 'iTag' ] and self._iTagZeroTrick : obsSetP2VV.append(tagDecision)

        sel   = Category( 'sel',             Title = 'Selection',        Observable = True, States = { 'selected' : +1 } )
        trig  = Category( 'triggerDecision', Title = 'Trigger Decision', Observable = True, States = { 'selected' : +1 } )

        observables = dict(  time        = time
                           , cpsi        = angles[0]
                           , ctheta      = angles[1]
                           , phi         = angles[2]
                           , iTag        = iTag if nominalPdf or not self._iTagZeroTrick else tagDecision
                           , tagDecision = tagDecision
                           , estWTag     = estWTag
                           , tagCat      = tagCat
                           , BMass       = BMass
                           , mpsi        = mpsi
                           , mphi        = mphi
                           , timeRes     = timeRes
                           , sel         = sel
                           , trig        = trig
                          )

        obsSetNTuple = [ time ] + angles +  [ BMass, mpsi, mphi, timeRes ] + [ tagDecision, estWTag, tagCat ] + [ sel, trig ]


        ###################################################################################################################################
        ## read data ##
        ###############

        if nTupleFile :
            from P2VVGeneralUtils import readData
            self._data = readData(  filePath = nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSetNTuple
                                  , Rename = 'JpsiphiData' )

        else :
            self._data = None


        ###################################################################################################################################
        ## initialize PDF component objects ##
        ######################################

        nEvents     = self._data.sumEntries() if self._data else numEvents
        nSignal     = nEvents * sigFrac
        nBackground = nEvents * ( 1. - sigFrac )

        from RooFitWrappers import Component
        if SFit :
            signalComps  = Component( 'signal',  [ ]                                       )
            sigMassComps = Component( 'sigMass', [ ], Yield = ( nSignal,     0., nEvents ) )
            bkgMassComps = Component( 'bkgMass', [ ], Yield = ( nBackground, 0., nEvents ) )
        else :
            signalComps     = Component( 'signal', [ ], Yield = ( nSignal,     0., nEvents ) )
            backgroundComps = Component( 'bkg'   , [ ], Yield = ( nBackground, 0., nEvents ) )


        ###################################################################################################################################
        ## build mass PDFs ##
        #####################

        # build the signal and background mass PDFs
        from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass, LP2011_Background_Mass as BackgroundBMass
        signalBMass     = SignalBMass(     Name = 'sig_m', mass = BMass )
        backgroundBMass = BackgroundBMass( Name = 'bkg_m', mass = BMass )

        from RooFitWrappers import buildPdf
        if SFit :
            sigMassComps += signalBMass.pdf()
            bkgMassComps += backgroundBMass.pdf()
            self._massPdf = buildPdf( [ sigMassComps, bkgMassComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )

        else :
            signalComps     += signalBMass.pdf()
            backgroundComps += backgroundBMass.pdf()
            self._massPdf = buildPdf( [ signalComps, backgroundComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )


        ###################################################################################################################################
        ## compute S-weights and create signal and background data sets ##
        ##################################################################

        if self._data :
            print 120 * '='
            print 'Bs2Jpsiphi_PdfBuilder: computing S-weights'

            # compute S-weights
            from P2VVGeneralUtils import SData, splot
            self._massPdf.fitTo( self._data, **fitOpts )
            #for par in self._massPdf.Parameters() : par.setConstant( not par.getAttribute('Yield') )
            self._SWeightData = SData( Pdf = self._massPdf, Data = self._data, Name = 'massSData' )

            # create signal and background data sets with S-weights
            self._sigSWeightData = self._SWeightData.data( 'sigMass' if SFit else 'signal' )
            self._bkgSWeightData = self._SWeightData.data( 'bkgMass' if SFit else 'bkg'    )

            self._sigSWeightData.Print()
            self._bkgSWeightData.Print()
            print 120 * '=' + '\n'

            # create signal and background data sets with side band ranges
            self._sigRangeData = self._data.reduce( CutRange = 'Signal'       )
            self._bkgRangeData = self._data.reduce( CutRange = 'LeftSideBand' )
            self._bkgRangeData.append( self._data.reduce( CutRange = 'RightSideBand' ) )

            if makePlots :
                # plot mass distributions
                self._massCanv = TCanvas( 'massCanv', 'B mass' )
                for ( pad, frameRange, nBins, plotTitle )\
                      in zip(  self._massCanv.pads( 2, 2, lambda pad : pad != 2 )
                             , [ 'Signal', 'LeftSideBand', 'RightSideBand' ]
                             , numBMassBins
                             , [  BMass.GetTitle() + ' mass fit - signal'
                                , BMass.GetTitle() + ' mass fit - left side band'
                                , BMass.GetTitle() + ' mass fit - right side band'
                               ]
                            ) :
                    plot(  pad, BMass, self._data, self._massPdf
                         , frameOpts  = dict( Range = frameRange, Bins = nBins, Title = plotTitle )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                   )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2                    )
                         , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                                         , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                                        }
                        )

        else :
            self._sigSWeightData = None
            self._bkgSWeightData = None
            self._sigRangeData   = None
            self._bkgRangeData   = None


        ###################################################################################################################################
        ## build tagging categories ##
        ##############################

        if nominalPdf or not self._iTagZeroTrick :
            # build tagging categories
            from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
            if nominalPdf or tagConds in [ 'estWTag', 'all' ] :
                self._tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = self._sigSWeightData, estWTag = estWTag
                                                  , wTagP0Constraint = True, wTagP1Constraint = True )
            else :
                self._tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = self._sigSWeightData, estWTagName = estWTag.GetName()
                                                  , TagCats = tagCats, NumSigmaTagBins = 1.
                                                  , wTagP0Constraint = True, wTagP1Constraint = True )
                if tagConds == 'tagCat' : self._tagCats.addConditional( self._tagCats['tagCat'] )

            if nominalPdf or tagConds in [ 'all', 'iTag' ] : self._tagCats.addConditional( iTag )

            tagCatP2VV = self._tagCats['tagCat']
            tagCatP2VV.setIndex(1)
            observables[tagCatP2VV.GetName()] = tagCatP2VV
            if not nominalPdf and not tagConds in [ 'all', 'tagCat', 'estWTag' ] : obsSetP2VV.append(tagCatP2VV)

            # add tagging category to data sets
            from P2VVGeneralUtils import addTaggingObservables
            for data in [ self._data, self._sigSWeightData, self._bkgSWeightData, self._sigRangeData, self._bkgRangeData ] :
                if data:
                    addTaggingObservables( data, iTag.GetName(), tagCatP2VV.GetName(), tagDecision.GetName(), estWTag.GetName()
                                          , self._tagCats['tagCats'] )


        ###################################################################################################################################
        ## time acceptance function ##
        ##############################

        if multiplyByTimeEff in [ 'all', 'signal', 'background' ] :
            from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
            timeAcceptance = TimeAcceptance( time = time, Input = timeEffHistFile, Histogram = timeEffHistName )


        ###################################################################################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################

        # transversity amplitudes
        if nominalPdf or amplitudeParam == 'phasesSWaveFrac' :
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(polarSWave = True if nominalPdf else polarSWave)

        else :
            from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes()

        # B lifetime
        from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        dGammaVar = dict( Name = 'dGamma' )
        if blind : dGammaVar['Blind'] = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
        lifetimeParams = LifetimeParams( dGamma = dGammaVar, dMConstraint = True )

        if nominalPdf or eventTimeRes :
            from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as TimeResolution
            timeResModel = TimeResolution( time = time, sigmat = timeRes, timeResSFConstraint = True )
        else :
            from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            timeResModel = TimeResolution( time = time, timeResSFConstraint = True )

        # CP violation parameters
        if carthLambdaCP : 
            from P2VVParameterizations.CPVParams import LambdaCarth_CPParam as CPParam
            ReLambdaCPVar = dict( Name = 'ReLambdaCP' )
            ImLambdaCPVar = dict( Name = 'ImLambdaCP' )
            if blind: ReLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsGoofyMoriond2012', 0.3 )
            if blind: ImLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsPlutoMoriond2012', 0.3 )
            self._lambdaCP = CPParam( ReLambdaCP = ReLambdaCPVar, ImLambdaCP = ImLambdaCPVar )

        else :
            from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
            phiCPVar = dict( Name = 'phiCP' )
            if blind: phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
            self._lambdaCP = CPParam( phiCP = phiCPVar )

        # coefficients for time functions
        from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

        # tagging parameters
        if not nominalPdf and self._iTagZeroTrick :
            from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams(  estWTag = estWTag, p0 = dict( Name = 'wTagP0' ), p1 = dict( Name = 'wTagP1' )
                                                , p0Constraint = True, p1Constraint = True )

            args = dict(  iTag     = tagDecision
                        , dilution = self._taggingParams['dilution']
                        , ADilWTag = self._taggingParams['ADilWTag']
                        , avgCEven = self._taggingParams['avgCEven']
                        , avgCOdd  = self._taggingParams['avgCOdd']
                       )

        else :
            from P2VVParameterizations.FlavourTagging import CatDilutionsCoefAsyms_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams( AProd = 0., ANorm = -self._lambdaCP['C'].getVal(), **self._tagCats.tagCatsDict() )

            args = dict(  tagCat      = tagCatP2VV
                        , iTag        = iTag
                        , dilutions   = self._taggingParams['dilutions']
                        , ADilWTags   = self._taggingParams['ADilWTags']
                        , avgCEvens   = self._taggingParams['avgCEvens']
                        , avgCOdds    = self._taggingParams['avgCOdds']
                        , tagCatCoefs = self._taggingParams['tagCatCoefs']
                       )

        args = dict(  time                   = time
                    , tau                    = lifetimeParams['MeanLifetime']
                    , dGamma                 = lifetimeParams['dGamma']
                    , dm                     = lifetimeParams['dM']
                    , coshCoef               = timeBasisCoefs['cosh']
                    , sinhCoef               = timeBasisCoefs['sinh']
                    , cosCoef                = timeBasisCoefs['cos']
                    , sinCoef                = timeBasisCoefs['sin']
                    , resolutionModel        = timeResModel['model']
                    , ConditionalObservables = timeResModel.conditionalObservables() + self._taggingParams.conditionalObservables()
                    , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                               + timeResModel.externalConstraints()\
                                               + self._taggingParams.externalConstraints()
                    , **args
                   )

        # build signal PDF
        from RooFitWrappers import BTagDecay
        sig_t_angles_tagCat_iTag = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

        if angEffMomentsFile :
            # multiply signal PDF with angular efficiency
            print 'Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with angular efficiency moments from file "%s"' % angEffMomentsFile

            from P2VVGeneralUtils import RealMomentsBuilder
            moments = RealMomentsBuilder()
            moments.appendPYList( self._angleFuncs.angles, [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not nominalPdf \
                                                else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
                                )
            moments.read(angEffMomentsFile)
            moments.Print()
            sig_t_angles_tagCat_iTag = moments * sig_t_angles_tagCat_iTag


        if multiplyByTimeEff in [ 'all', 'signal' ] :
            # multiply signal PDF with time acceptance
            sig_t_angles_tagCat_iTag = timeAcceptance * sig_t_angles_tagCat_iTag

        signalComps += sig_t_angles_tagCat_iTag

        # print tagging category distribution for signal and background
        if self._sigSWeightData and ( nominalPdf or not self._iTagZeroTrick ) :
            print 'Bs2Jpsiphi_PdfBuilder: distribution in tagging category for signal:'
            self._sigSWeightData.table(tagCatP2VV).Print('v')
            print 'Bs2Jpsiphi_PdfBuilder: distribution in tagging category for background:'
            self._bkgSWeightData.table(tagCatP2VV).Print('v')


        ###################################################################################################################################
        ## build PDF for estimated wrong-tag probability ##
        ###################################################

        if self._data and self._sigSWeightData and self._bkgSWeightData and makePlots :
            tempTagCat = tagCat if not nominalPdf and self._iTagZeroTrick else tagCatP2VV

            # build PDF for estimated wrong-tag probability
            print 'Bs2Jpsiphi_PdfBuilder: building PDF for estimated wrong-tag probability'
            from RooFitWrappers import HistPdf
            self._estWTagData = self._sigSWeightData.reduce( '%s > 0' % tempTagCat.GetName() )
            self._sig_bkg_estWTag = HistPdf(  Name = 'sig_bkg_estWTag'
                                            , Observables = [ estWTag ]
                                            , Binning = { estWTag : numEstWTagBins }
                                            , Data = self._estWTagData
                                           )

            # get normalization correction for tagged events
            untagFrac    = self._data.table(tempTagCat).getFrac('Untagged')
            untagFracSig = self._sigSWeightData.table(tempTagCat).getFrac('Untagged')
            untagFracBkg = self._bkgSWeightData.table(tempTagCat).getFrac('Untagged')

            # plot estimated wrong-tag probability for signal and for background
            self._estWTagCanv = TCanvas( 'estWTagCanv', 'Estimated wrong-tag probability' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanv.pads( 1, 1 ) if SFit else self._estWTagCanv.pads( 2, 2, lambda pad : pad != 2 )
                         , [ self._data, self._sigSWeightData, self._bkgSWeightData ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - untagFrac, 1. - untagFracSig, 1. - untagFracBkg ]
                        ) :
                plot(  pad, estWTag, data, self._sig_bkg_estWTag
                     , frameOpts  = dict( Bins = nBins, Title = estWTag.GetTitle() + plotTitle, Range = ( 0., 0.499999 ) )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                              )
                     , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2, Normalization = norm                         )
                    )


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        if not SFit or makePlots :
            from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
            backgroundTime = BackgroundTime( Name = 'bkg_t', time = time, resolutionModel = timeResModel['model'] )
            self._bkg_t = backgroundTime.pdf()

            if multiplyByTimeEff in [ 'all', 'background' ] :
                # multiply background time PDF with time acceptance
                self._bkg_t = timeAcceptance * self._bkg_t

            backgroundComps += self._bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit or makePlots :
            print 'Bs2Jpsiphi_PdfBuilder: determining angular shape of background from %s data'\
                  % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )

            bkgAngleData = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
            if bkgAngleData and bkgAnglePdf == 'histPdf' :
                # create a HistPdf from background data
                from RooFitWrappers import HistPdf
                nBins = [ 5, 7, 9 ] if nominalPdf else numAngleBins
                self._bkg_angles = HistPdf(  Name = 'bkg_angles'
                                           , Observables = angles
                                           , Binning = { cpsi : nBins[0], ctheta : nBins[1], phi : nBins[2] }
                                           , Data = bkgAngleData
                                          )

            else :
                # define angular bins
                from math import pi
                from array import array
                if nominalPdf or transAngles :
                    cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
                    cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
                    phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )
                else :
                    cpsBinBounds = array( 'd', [ -1.,                      -0.6,      -0.2,      0.2,      0.6,                   1. ] )
                    cthBinBounds = array( 'd', [ -1., -0.95, -0.90, -0.85, -0.6,      -0.2,      0.2,      0.6, 0.85, 0.90, 0.95, 1. ] )
                    phiBinBounds = array( 'd', [ -pi,                      -0.6 * pi, -0.2 * pi, 0.2 * pi, 0.6 * pi,              pi ] )

                cpsNumBins = len(cpsBinBounds) - 1
                cthNumBins = len(cthBinBounds) - 1
                phiNumBins = len(phiBinBounds) - 1

                from ROOT import RooBinning
                cpsBins = RooBinning( cpsNumBins, cpsBinBounds, 'bkg_cpsBins' )
                cthBins = RooBinning( cthNumBins, cthBinBounds, 'bkg_cthBins' )
                phiBins = RooBinning( phiNumBins, phiBinBounds, 'bkg_phiBins' )
                cpsi.setBinning(   cpsBins, 'bkg_cpsBins' )
                ctheta.setBinning( cthBins, 'bkg_cthBins' )
                phi.setBinning(    phiBins, 'bkg_phiBins' )

                # create bin coefficients
                self._bkgAngCoefs = [ RealVar(  'bkg_angBin_%d_%d_%d' % ( bin0, bin1, bin2 )
                                              , Title    = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                                              , Value    = 1. / cpsNumBins / cthNumBins / phiNumBins
                                              , MinMax   = ( 0., 1. )
                                              , Constant = True
                                             )\
                                      if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                                      for bin2 in range( phiNumBins ) for bin1 in range( cthNumBins ) for bin0 in range( cpsNumBins )
                                    ]
                del self._bkgAngCoefs[0]

                # create a BinnedPdf
                from RooFitWrappers import BinnedPdf
                self._bkg_angles = BinnedPdf(  Name = 'bkg_angles'
                                             , Observables = angles
                                             , Binnings = [ cpsBins, cthBins, phiBins ]
                                             , Coefficients = self._bkgAngCoefs
                                             , BinIntegralCoefs = True
                                            )

                sumWeights = 0.
                sumBinWeights = ( cpsNumBins * cthNumBins * phiNumBins - 1 ) * [ 0. ]
                if bkgAngleData :
                    # determine bin coefficient values
                    angleInitVals = [ angle.getVal() for angle in angles ]
                    for obsSet in bkgAngleData :
                        for angle in angles : angle.setVal( obsSet.getRealValue( angle.GetName() ) )
                        sumWeights += bkgAngleData.weight()
                        cpsBin = cpsi.getBin('bkg_cpsBins')
                        cthBin = ctheta.getBin('bkg_cthBins')
                        phiBin = phi.getBin('bkg_phiBins')
                        bin = cpsBin + cthBin * cpsNumBins + phiBin * cpsNumBins * cthNumBins - 1
                        if bin >= 0 : sumBinWeights[ bin ] += bkgAngleData.weight()
                    for angle, val in zip( angles, angleInitVals ) : angle.setVal(val)

                # set bin coefficient values
                for coef, weight in zip( self._bkgAngCoefs, sumBinWeights ) :
                    if bkgAngleData :
                        value = weight / sumWeights
                        assert value >= 0.,\
                            'Bs2Jpsiphi_PdfBuilder: background angular PDF coefficient \"%s\" has negative value: %f'\
                            % (coef.GetName(), value)
                        coef.setVal(value)

            if not SFit : backgroundComps += self._bkg_angles

            if self._bkgRangeData and self._bkgSWeightData and makePlots :
                # plot background angles
                self._bkgAnglesCanv = TCanvas( 'bkgAnglesCanv', 'Background Decay Angles' )
                for ( pad, obs, data, nBins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesCanv.pads( 3, 2 )
                             , 2 * angles
                             , 3 * ( self._bkgRangeData, ) + 3 * ( self._bkgSWeightData, )
                             , 2 * numAngleBins
                             , [ angle.GetTitle() + ' - mass side bands' for angle in angles ]
                             + [ angle.GetTitle() + ' - mass S-weights'  for angle in angles ]
                             , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle
                         , frameOpts  = dict( Bins = nBins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
                        )


        ###################################################################################################################################
        ## build background tagging PDF ##
        ##################################

        if not SFit :
            if not nominalPdf and self._iTagZeroTrick :
                if tagConds not in [ 'iTag', 'all' ] :
                    self._bkg_tagCat_iTag = None
                    backgroundComps[tagDecision] = self._bkg_tagCat_iTag

            else :
                bkgTaggingData = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
                if not nominalPdf and bkgTaggingData and bkgTaggingPdf == 'histPdf' :
                    print 'Bs2Jpsiphi_PdfBuilder: creating background tagging PDFs from %s data'\
                          % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )
                    from RooFitWrappers import HistPdf
                    self._bkg_tagCat_iTag = HistPdf( Name = 'bkg_tagCat_iTag', Observables = [ tagCatP2VV, iTag ], Data = bkgTaggingData )
                    backgroundComps += self._bkg_tagCat_iTag

                else :
                    print 'Bs2Jpsiphi_PdfBuilder: determining tagged category coefficients from %s'\
                           % ( 'signal PDF' if bkgTaggingPdf == 'signal'\
                               else ( 'B mass side band data' if massRangeBackground else 'B mass S-weight data' ) )

                    # get number of events in each category
                    from RooFitWrappers import ArgSet, FormulaVar
                    if bkgTaggingData :
                        catTagTable = bkgTaggingData.table( ArgSet( 'catTagSet', [ tagCatP2VV, iTag ] ) )
                        nUntB    = catTagTable.get( '{Untagged;B}'    )
                        nUntBbar = catTagTable.get( '{Untagged;Bbar}' )

                        nTagB    = [ 0. for cat in range( tagCatP2VV.numTypes() ) ]
                        nTagBbar = [ 0. for cat in range( tagCatP2VV.numTypes() ) ]
                        for catIter in range( 1, tagCatP2VV.numTypes() ) :
                            taggedCatName = ( 'TagCat%d' % catIter ) if tagCatP2VV.numTypes() > 2 else 'Tagged'
                            nTagB[catIter]    += catTagTable.get( '{%s;B}'    % taggedCatName )
                            nTagBbar[catIter] += catTagTable.get( '{%s;Bbar}' % taggedCatName )
                            nTagB[0]    += nTagB[catIter]
                            nTagBbar[0] += nTagBbar[catIter]

                    # check order of types in tag
                    tagRevOrder = False
                    for iter, type in enumerate(iTag) :
                        if iter == 0 and type.getVal() < 0 : tagRevOrder = True

                    if nominalPdf or bkgTaggingPdf == 'signal' :
                        tagCatCoefNames = [ coef.GetName() for coef in sig_t_angles_tagCat_iTag.tagCatCoefs() ]
                        AUntagged = RealVar(  'bkgAUnt'
                                            , Title  = 'Untagged-tagged asymmetry in background tagging category coefficients'
                                            , Value  = bkgTaggingData.table(tagCatP2VV).getFrac('Untagged')\
                                                       / ws[ tagCatCoefNames[0] ].getVal() - 1. if bkgTaggingData\
                                                       else 0.
                                            , Error  = 0.1
                                            , MinMax = ( -1., 1. )
                                           )

                        if bkgTaggingData :
                            AUnt = float( nUntB - nUntBbar ) / float( nUntB + nUntBbar )
                            ATag = float( nTagB[0] - nTagBbar[0]) / float( nTagB[0] + nTagBbar[0] )

                        else :
                            AUnt = 0.
                            ATag = 0.

                        AUntaggedBBbar = RealVar(  'bkgAUntBBbar'
                                                 , Title  = 'B-Bbar asymmetry in untagged background'
                                                 , Value  = -AUnt if tagRevOrder else AUnt
                                                 , Error  = 0.1
                                                 , MinMax = ( -1., 1. )
                                                )

                        ATaggedBBbar = RealVar(  'bkgATagBBbar'
                                               , Title  = 'B-Bbar asymmetry in tagged background'
                                               , Value  = -ATag if tagRevOrder else ATag
                                               , Error  = 0.1
                                               , MinMax = ( -1., 1. )
                                              )

                        untaggedBCoef  = FormulaVar( 'bkgUntaggedBCoef',  '0.5*(1.+@0)*(1.-@1)'
                                                    , [ AUntagged, AUntaggedBBbar                         ] )
                        taggedBCoef    = FormulaVar( 'bkgTaggedBCoef',    '0.5*(1.-@2/(1.-@2)*@0)*(1.+@1)'
                                                    , [ AUntagged, ATaggedBBbar, ws[ tagCatCoefNames[0] ] ] )
                        taggedBbarCoef = FormulaVar( 'bkgTaggedBbarCoef', '0.5*(1.-@2/(1.-@2)*@0)*(1.-@1)'
                                                    , [ AUntagged, ATaggedBBbar, ws[ tagCatCoefNames[0] ] ] )

                        tagBinCoefs = [ ]
                        from RooFitWrappers import Product
                        for binIter in range( 1, tagCatP2VV.numTypes() * 2 ) :
                            cat = binIter % tagCatP2VV.numTypes()
                            tagBinCoefs.append( Product(  'bkgTagBinCoef%d' % binIter
                                                        , [ ws[ tagCatCoefNames[cat] ], untaggedBCoef if cat == 0\
                                                            else ( taggedBbarCoef if binIter < tagCatP2VV.numTypes() else taggedBCoef ) ]
                                                       )
                                              )

                    else :
                        ABBbars = [ ]
                        for catIter in range( tagCatP2VV.numTypes() ) :
                            if bkgTaggingData and catIter == 0 :
                                ABBbar = float( nUntB - nUntBbar ) / float( nUntB + nUntBbar )
                            elif bkgTaggingData :
                                ABBbar = float( nTagB[catIter] - nTagBbar[catIter] ) / float( nTagB[catIter] + nTagBbar[catIter] )
                            else :
                                ABBbar = 0.

                            ABBbars.append( RealVar(  'bkgABBbar%d' % catIter
                                                    , Title    = 'B-Bbar asymmetry in background %d' % catIter
                                                    , Value    = -ABBbar if tagRevOrder else ABBbar
                                                    , MinMax   = ( -1., 1. )
                                                    , Constant = True
                                                   )
                                          )

                        tagCatCoefs = [ ]
                        for catIter in range( tagCatP2VV.numTypes() ) :
                            if bkgTaggingData and catIter == 0. :
                                tagCatCoef = bkgTaggingData.table(tagCatP2VV).getFrac('Untagged')
                            elif bkgTaggingData :
                                taggedCatName = ( 'TagCat%d' % catIter ) if tagCatP2VV.numTypes() > 2 else 'Tagged'
                                tagCatCoef = bkgTaggingData.table(tagCatP2VV).getFrac(taggedCatName)
                            else :
                                tagCatCoef = 0.

                            tagCatCoefs.append( RealVar(  'bkgTagCatCoef%d' % catIter
                                                        , Title    = 'Background tagging category coefficient %d' % catIter
                                                        , Value    = tagCatCoef
                                                        , MinMax   = ( -1., 1. )
                                                        , Constant = True
                                                       )
                                              )

                        tagBinCoefs = [ ]
                        for binIter in range( 1, tagCatP2VV.numTypes() * 2 ) :
                            cat = binIter % tagCatP2VV.numTypes()
                            tagBinCoefs.append( FormulaVar(  'bkgTagBinCoef%d' % binIter
                                                           , '0.5*@0*(1%s@1)' % ( '+' if binIter < tagCatP2VV.numTypes() else '-' )
                                                           , [ tagCatCoefs[cat], ABBbars[cat] ]
                                                          )
                                              )

                    # build binned PDF
                    from RooFitWrappers import BinnedPdf
                    self._bkg_tagCat_iTag = BinnedPdf( 'bkg_tagCat_iTag', Categories = [ tagCatP2VV, iTag ], Coefficients = tagBinCoefs )
                    backgroundComps += self._bkg_tagCat_iTag


        ###################################################################################################################################
        ## build full PDF ##
        ####################

        from RooFitWrappers import buildPdf
        if SFit :
            pdf = buildPdf( [ signalComps ], Observables = obsSetP2VV, Name = 'Jpsiphi' )

        else :
            if components == 'signal' :
                pdf = buildPdf( [ signalComps                  ], Observables = obsSetP2VV, Name = 'JpsiphiSig' )
            elif components == 'background' :
                pdf = buildPdf( [ backgroundComps              ], Observables = obsSetP2VV, Name = 'JpsiphiBkg' )
            else :
                pdf = buildPdf( [ signalComps, backgroundComps ], Observables = obsSetP2VV, Name = 'Jpsiphi'    )

        assert not pdfConfig, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: superfluous arguments found: %s' % pdfConfig
        PdfBuilder.__init__( self, pdf, observables, { } )
