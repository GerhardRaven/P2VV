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

    def getParametersFromPdf( self, pdf, data ) :
        for par in pdf.getParameters(data) :
            self.addParameter( par.GetName(), (  par.getVal()
                                               , par.getError()
                                               , par.getMin()
                                               , par.getMax()
                                               , False if par.isConstant() else True
                                              )
                             )

    def setParametersInPdf( self, pdf ) :
        for par in pdf.getVariables() :
            if par.GetName() in self._parameters.keys() :
                par.setVal(      self._parameters[ par.GetName() ][0]                    )
                par.setError(    self._parameters[ par.GetName() ][1]                    )
                par.setMin(      self._parameters[ par.GetName() ][2]                    )
                par.setMax(      self._parameters[ par.GetName() ][3]                    )
                par.setConstant( False if self._parameters[ par.GetName() ][4] else True )

    def readParametersFromFile( self, filePath = 'parameters', **kwargs ) :
        # get file path
        filePath = filePath.strip()

        # open file
        try :
          parFile = open( filePath, 'r' )
        except :
          raise RuntimeError( 'P2VV - ERROR: PdfConfiguration.readParametersFromFile: unable to open file \"%s\"' % filePath )

        # get name requirements
        import re
        nameExpr = re.compile( kwargs.pop('Names') ) if 'Names' in kwargs else None

        # loop over lines and read parameters
        numPars = 0
        while True :
            # read next line
            line = parFile.readline()
            if not line : break

            # check for empty or comment lines
            line = line.strip()
            if not line or line[0] == '#' : continue

            # check moment format
            line = line.split()
            if len(line) != 6 : continue

            # check name
            if nameExpr and not nameExpr.match(line[0]) : continue

            try :
              parVal   = float(line[1])
              parErr   = float(line[2])
              parMin   = float(line[3])
              parMax   = float(line[4])
              parFloat = bool( 1 if line[5] == 'True' else 0 )
            except :
              continue

            # set parameter values
            self.addParameter( line[0], ( parVal, parErr, parMin, parMax, parFloat ) )
            numPars += 1

        parFile.close()

        print 'P2VV - INFO: PdfConfiguration.readParametersFromFile: %d parameter%s read from file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )

    def writeParametersToFile( self, filePath = 'parameters', **kwargs ) :
        # get file path and name
        filePath = filePath.strip()
        fileName = filePath.split('/')[-1]

        # open file
        try :
            parFile = open( filePath, 'w' )
        except :
            raise RuntimeError( 'P2VV - ERROR: PdfConfiguration.writeParametersToFile: unable to open file \"%s\"' % filePath )

        # get maximum length of parameter name
        maxLenName = 13
        for parName in self._parameters.keys() : maxLenName = max( len(parName), maxLenName )

        # get name requirements
        import re
        names = kwargs.pop( 'Names', None )
        nameExpr = re.compile(names) if names else None

        # get floating/fixed
        floating = kwargs.pop( 'Floating', None )
        if floating not in [ True, False ] : floating = None

        # write parameters to content string
        cont = '# %s: parameters\n' % fileName\
             + '# name requirement: \'{0}\'\n'.format( names if names else '' )\
             + '# floating:         \'{0}\'\n'.format( 'True' if floating == True else ( 'False' if floating == False else '' ) )\
             + '#\n'\
             + '# ' + '-' * (79 + maxLenName) + '\n'\
             + ( '# {0:<%s}   {1:<14}   {2:<13}   {3:<14}   {4:<14}   {5:<}\n' % maxLenName )\
                 .format( 'parameter', 'value', 'error', 'min', 'max', 'floating?' )\
             + '# ' + '-' * (79 + maxLenName) + '\n'

        numPars = 0
        for parName in sorted( self._parameters.keys() ) :
            if nameExpr and not nameExpr.match(parName) : continue

            parVals = self._parameters[parName]
            if ( floating == True and not parVals[4] ) or ( floating == False and parVals[4] ) : continue

            cont += ( '  {0:<%s}   {1:<+14.8g}   {2:<13.8g}   {3:<+14.8g}   {4:<+14.8g}   {5:<}\n' % maxLenName )\
                      .format( parName, parVals[0], parVals[1], parVals[2], parVals[3], 'True' if parVals[4] else 'False' )
            numPars += 1

        cont += '# ' + '-' * (79 + maxLenName) + '\n'

        # write content to file
        parFile.write(cont)
        parFile.close()

        print 'P2VV - INFO: PdfConfiguration.writeParametersToFile: %d parameter%s written to file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )


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
        self['transversityAngles']   = False
        self['bkgAnglePdf']          = 'histPdf'
        self['sigTaggingPdf']        = 'tagUntag'          # 'histPdf' / 'tagUntag' / 'tagCats'
        self['bkgTaggingPdf']        = 'tagUntagRelative'  # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        self['multiplyByTimeEff']    = ''                  # 'all' / 'signal'
        self['parameterizeKKMass']   = ''  # '' / 'functions' / 'simultaneous'
        self['ambiguityParameters']  = False
        self['KKMassBinBounds']      = [ 1020. - 12., 1020. + 12. ]
        self['SWaveAmplitudeValues'] = (  [ 0.026 ], [ 0. ] )
        self['CSPValues']            = [ 0.4976 ]

        self['sameSideTagging']    = True
        self['conditionalTagging'] = False
        self['continuousEstWTag']  = False
        self['numEstWTagBins']     = 100
        self['constrainTagging']   = True

        self['iTagZeroTrick'] = False
        self['iTagStates'] = { }                         # { } / { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }

        self['eventTimeResolution'] = True
        self['numTimeResBins']      = 100

        self['signalFraction'] = 0.67
        self['massRangeBackground'] = False

        self['amplitudeParam'] = 'phasesSWaveFrac'       # 'phases' / 'phasesSWaveFrac' / 'ReIm' / 'bank'
        self['ASParam']        = 'deltaPerp'             # 'delta0' / 'deltaPerp' / 'ReIm' / 'Mag2ReIm' / 'Mag2ReImPerp'
        self['AparParam']      = 'cos'                   # 'phase' / 'ReIm' / 'Mag2ReIm' / 'cos' / 'real'

        self['constrainDeltaM'] = True

        self['lambdaCPParam'] = 'lambSqPhi'

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
        transAngles       = pdfConfig.pop('transversityAngles')
        bkgAnglePdf       = pdfConfig.pop('bkgAnglePdf')
        sigTaggingPdf     = pdfConfig.pop('sigTaggingPdf')
        bkgTaggingPdf     = pdfConfig.pop('bkgTaggingPdf')
        multiplyByTimeEff = pdfConfig.pop('multiplyByTimeEff')
        paramKKMass       = pdfConfig.pop('parameterizeKKMass')
        numBMassBins      = pdfConfig.pop('numBMassBins')
        ambiguityPars     = pdfConfig.pop('ambiguityParameters')
        KKMassBinBounds   = pdfConfig.pop('KKMassBinBounds')
        SWaveAmpVals      = pdfConfig.pop('SWaveAmplitudeValues')
        CSPValues         = pdfConfig.pop('CSPValues')

        self._iTagZeroTrick = pdfConfig.pop('iTagZeroTrick')
        iTagStates = pdfConfig.pop('iTagStates')
        if not iTagStates : iTagStates = { 'B' : +1, 'Bbar' : -1 }
        if +1 not in iTagStates.values() or -1 not in iTagStates.values() or 0 in iTagStates.values() : self._iTagZeroTrick = True
        if self._iTagZeroTrick : iTagStatesDecision = iTagStates
        else                   : iTagStatesDecision = { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }

        SSTagging        = pdfConfig.pop('sameSideTagging')
        condTagging      = pdfConfig.pop('conditionalTagging')
        contEstWTag      = pdfConfig.pop('continuousEstWTag')
        numEstWTagBins   = pdfConfig.pop('numEstWTagBins')
        constrainTagging = pdfConfig.pop('constrainTagging')
        condTagging = True if contEstWTag else condTagging

        eventTimeRes = pdfConfig.pop('eventTimeResolution')
        numTimeResBins = pdfConfig.pop('numTimeResBins')

        sigFrac = pdfConfig.pop('signalFraction')
        massRangeBackground = pdfConfig.pop('massRangeBackground')

        amplitudeParam = pdfConfig.pop('amplitudeParam')
        ASParam        = pdfConfig.pop('ASParam')
        AparParam      = pdfConfig.pop('AparParam')

        if not paramKKMass :
            if not KKMassBinBounds : KKMassBinBounds = [ 1020. - 12., 1020., 1020. + 12. ]
        elif ambiguityPars :
            if amplitudeParam == 'bank' and ASParam != 'ReIm' :
                from math import pi
                for phaseIter, phase in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][phaseIter] = pi - phase
            elif amplitudeParam == 'phases' and ASParam in [ 'Mag2ReIm', 'ReIm' ] :
                for ImIter, Im in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][ImIter] = -Im

        constrainDeltaM = pdfConfig.pop('constrainDeltaM')

        lambdaCPParam = pdfConfig.pop('lambdaCPParam')

        if makePlots :
            # import plotting tools
            from P2VVLoad import ROOTStyle
            from P2VVGeneralUtils import plot
            from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed


        ###################################################################################################################################
        ## create variables (except for tagging category) ##
        ####################################################

        # RooObject wrappers
        from RooFitWrappers import RooObject, ConstVar, RealVar, Category
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
        timeRes = RealVar(  'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, Value = 0.10, MinMax = (0.0001, 0.12)
                          , nBins = numTimeResBins )
        timeRes.setBins( numTimeResBins, 'cache' )

        cpsi   = self._angleFuncs.angles['cpsi']
        ctheta = self._angleFuncs.angles['ctheta']
        phi    = self._angleFuncs.angles['phi']

        if nominalPdf or not self._iTagZeroTrick :
            iTagOS = Category( 'iTagOS', Title = 'Initial state flavour tag opposite side', Observable = True, States = iTagStates )
            iTagSS = Category( 'iTagSS', Title = 'Initial state flavour tag same side',     Observable = True, States = iTagStates )
        estWTagComb = RealVar( 'tagomega',    Title = 'Estimated wrong tag probability OS/SS combination', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ), nBins = numEstWTagBins )
        estWTagOS   = RealVar( 'tagomega_os', Title = 'Estimated wrong tag probability opposite side', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ), nBins = numEstWTagBins )
        estWTagSS   = RealVar( 'tagomega_ss', Title = 'Estimated wrong tag probability same side', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ), nBins = numEstWTagBins )
        estWTagOS.setBins( numEstWTagBins, 'cache' )
        estWTagSS.setBins( numEstWTagBins, 'cache' )

        BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True
                        , Value = 5368., MinMax = ( 5200., 5550. ), nBins = numBMassBins[0] + numBMassBins[1] + numBMassBins[2]
                        ,  Ranges = dict(  LeftSideBand  = ( None,  5330. )
                                         , Signal        = ( 5330., 5410. )
                                         , RightSideBand = ( 5410., None  )
                                        )
                       )

        angles = [ cpsi, ctheta, phi ]
        obsSetP2VV = [ time ] + angles
        if not self._iTagZeroTrick :
            obsSetP2VV.append(iTagOS)
            #obsSetP2VV.append(iTagSS)
        if not SFit :
            obsSetP2VV.append(BMass)

        # ntuple variables
        mumuMass = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. )
                           , nBins =  32 )
        KKMass = RealVar( 'mdau2', Title = 'M(KK)', Unit = 'MeV', Observable = True, MinMax = ( KKMassBinBounds[0], KKMassBinBounds[-1] )
                         , nBins =  32 )
        if paramKKMass == 'functions' : obsSetP2VV.append(KKMass)

        tagDecisionComb = Category( 'tagdecision',    Title = 'Tag decision OS/SS combination', Observable = True, States = iTagStatesDecision )
        tagDecisionOS   = Category( 'tagdecision_os', Title = 'Tag decision opposite side',     Observable = True, States = iTagStatesDecision )
        tagDecisionSS   = Category( 'tagdecision_ss', Title = 'Tag decision same side',         Observable = True, States = iTagStatesDecision )
        tagCatOS = Category( 'tagcat_os',   Title = 'Tagging category opposite side', Observable = True
                            , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                           )

        sel   = Category( 'sel',                     Title = 'Selection',                 Observable = True, States = { 'selected' : +1 } )
        trig  = Category( 'triggerDecisionUnbiased', Title = 'Trigger Decision Unbiased', Observable = True, States = { 'selected' : +1 } )

        muPlusTrackChi2 = RealVar( 'muplus_track_chi2ndof',  Title = 'mu+ track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
        muMinTrackChi2  = RealVar( 'muminus_track_chi2ndof', Title = 'mu- track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
        KPlusTrackChi2  = RealVar( 'Kplus_track_chi2ndof',   Title = 'K+ track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )
        KMinTrackChi2   = RealVar( 'Kminus_track_chi2ndof',  Title = 'K- track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )

        observables = dict(  time            = time
                           , cpsi            = angles[0]
                           , ctheta          = angles[1]
                           , phi             = angles[2]
                           , iTagOS          = iTagOS if nominalPdf or not self._iTagZeroTrick else tagDecisionOS
                           , iTagSS          = iTagSS if nominalPdf or not self._iTagZeroTrick else tagDecisionSS
                           , tagDecisionComb = tagDecisionComb
                           , tagDecisionOS   = tagDecisionOS
                           , estWTagComb     = estWTagComb
                           , estWTagOS       = estWTagOS
                           , estWTagSS       = estWTagSS
                           , tagCatOS        = tagCatOS
                           , BMass           = BMass
                           , mumuMass        = mumuMass
                           , KKMass          = KKMass
                           , timeRes         = timeRes
                           , sel             = sel
                           , trig            = trig
                           , muPlusTrackChi2 = muPlusTrackChi2
                           , muMinTrackChi2  = muMinTrackChi2
                           , KPlusTrackChi2  = KPlusTrackChi2
                           , KMinTrackChi2   = KMinTrackChi2
                          )

        obsSetNTuple = [ time ] + angles +  [ BMass, mumuMass, KKMass, timeRes ] + [ tagDecisionComb, estWTagComb ]\
                       + [ tagDecisionOS, estWTagOS, tagCatOS ] + [ tagDecisionSS, estWTagSS ]\
                       + [ sel, trig, muPlusTrackChi2, muMinTrackChi2, KPlusTrackChi2, KMinTrackChi2 ]


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
            self._signalComps  = Component( 'signal',  [ ]                                       )
            self._sigMassComps = Component( 'sigMass', [ ], Yield = ( nSignal,     0., nEvents ) )
            self._bkgMassComps = Component( 'bkgMass', [ ], Yield = ( nBackground, 0., nEvents ) )
        else :
            self._signalComps     = Component( 'signal', [ ], Yield = ( nSignal,     0., nEvents ) )
            self._backgroundComps = Component( 'bkg'   , [ ], Yield = ( nBackground, 0., nEvents ) )


        ###################################################################################################################################
        ## build B mass PDFs ##
        #######################

        # build the signal and background mass PDFs
        from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass, LP2011_Background_Mass as BackgroundBMass
        self._signalBMass     = SignalBMass(     Name = 'sig_m', mass = BMass )
        self._backgroundBMass = BackgroundBMass( Name = 'bkg_m', mass = BMass )

        from RooFitWrappers import buildPdf
        if SFit :
            self._sigMassComps += self._signalBMass.pdf()
            self._bkgMassComps += self._backgroundBMass.pdf()
            self._massPdf = buildPdf( [ self._sigMassComps, self._bkgMassComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )

        else :
            self._signalComps     += self._signalBMass.pdf()
            self._backgroundComps += self._backgroundBMass.pdf()
            self._massPdf = buildPdf( [ self._signalComps, self._backgroundComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )


        ###################################################################################################################################
        ## compute S-weights and create signal and background data sets ##
        ##################################################################

        if self._data :
            print 120 * '='
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: computing S-weights'

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
        ## build KK mass PDFs ##
        ########################

        if paramKKMass or makePlots :
            # create binning
            from array import array
            KKMassBinsArray = array( 'd', KKMassBinBounds )
            KKMassNBins = len(KKMassBinBounds) - 1

            from ROOT import RooBinning
            self._KKMassBinning = RooBinning( KKMassNBins, KKMassBinsArray, 'KKMassBinning' )
            KKMass.setBinning( self._KKMassBinning, 'KKMassBinning' )

            # build the signal and background KK mass PDFs
            from P2VVParameterizations.MassPDFs import Binned_MassPdf
            self._signalKKMass =     Binned_MassPdf( 'sig_mKK', KKMass, Binning = self._KKMassBinning, Data = self._sigSWeightData )
            self._backgroundKKMass = Binned_MassPdf( 'bkg_mKK', KKMass, Binning = self._KKMassBinning, Data = self._bkgSWeightData )
            if paramKKMass == 'functions' :
                self._signalComps += self._signalKKMass.pdf()
                if not SFit: self._backgroundComps += self._backgroundKKMass.pdf()

            if makePlots :
                self._KKMassCanv = TCanvas( 'KKMassCanv', 'KK Mass' )
                for ( pad, data, pdf, plotTitle )\
                      in zip(  self._KKMassCanv.pads( 2, 2 )
                             , [ self._sigSWeightData, self._bkgSWeightData ]
                             , [ self._signalKKMass.pdf(), self._backgroundKKMass.pdf() ]
                             , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                            ) :
                    plot(  pad, KKMass, data, None #pdf, logy = True
                         , frameOpts  = dict( Title = KKMass.GetTitle() + plotTitle )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )#, MarkerColor = kBlue, LineColor = kBlue   )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2    )
                        )

        ###################################################################################################################################
        ## build tagging categories ##
        ##############################

        if nominalPdf or not self._iTagZeroTrick :
            # build tagging categories opposite side
            from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
            if nominalPdf or contEstWTag :
                self._tagCatsOS = TaggingCategories(  tagCat = 'tagCatP2VVOS', DataSet = self._sigSWeightData, estWTag = estWTagOS
                                                    , wTagP0Constraint = True if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = True if nominalPdf else constrainTagging )
                self._tagCatsSS = TaggingCategories(  tagCat = 'tagCatP2VVSS', DataSet = self._sigSWeightData, estWTag = estWTagSS
                                                    , SameSide = True
                                                    , wTagP0Constraint = True if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = True if nominalPdf else constrainTagging )
            else :
                self._tagCatsOS = TaggingCategories(  tagCat = 'tagCatP2VVOS', DataSet = self._sigSWeightData
                                                    , estWTagName = estWTagOS.GetName()
                                                    , TagCats = tagCats, NumSigmaTagBins = 1.
                                                    , wTagP0Constraint = True if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = True if nominalPdf else constrainTagging )
                self._tagCatsSS = TaggingCategories(  tagCat = 'tagCatP2VVSS', DataSet = self._sigSWeightData
                                                    , SameSide = True
                                                    , estWTagName = estWTagSS.GetName()
                                                    , TagCats = tagCats, NumSigmaTagBins = 1.
                                                    , wTagP0Constraint = True if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = True if nominalPdf else constrainTagging )

            tagCatP2VVOS = self._tagCatsOS['tagCat']
            tagCatP2VVOS.setIndex(1)
            observables[tagCatP2VVOS.GetName()] = tagCatP2VVOS
            obsSetP2VV.append(tagCatP2VVOS)

            tagCatP2VVSS = self._tagCatsSS['tagCat']
            tagCatP2VVSS.setIndex(1)
            observables[tagCatP2VVSS.GetName()] = tagCatP2VVSS
            #obsSetP2VV.append(tagCatP2VVSS)

            if nominalPdf or condTagging :
                self._tagCatsOS.addConditional(tagCatP2VVOS)
                self._tagCatsOS.addConditional(iTagOS)
                self._tagCatsSS.addConditional(tagCatP2VVSS)
                self._tagCatsSS.addConditional(iTagSS)

            # add tagging categories to data sets
            from P2VVGeneralUtils import addTaggingObservables
            for data in [ self._data, self._sigSWeightData, self._bkgSWeightData, self._sigRangeData, self._bkgRangeData ] :
                if data:
                    addTaggingObservables( data, iTagOS.GetName(), tagCatP2VVOS.GetName(), tagDecisionOS.GetName(), estWTagOS.GetName()
                                          , self._tagCatsOS['tagCats'] )
                    addTaggingObservables( data, iTagSS.GetName(), tagCatP2VVSS.GetName(), tagDecisionSS.GetName(), estWTagSS.GetName()
                                          , self._tagCatsSS['tagCats'] )

            from P2VVParameterizations.FlavourTagging import Combined_TaggingCategories as CombTaggingCategories
            self._tagCatsComb = CombTaggingCategories( self._tagCatsOS, self._tagCatsSS )

            if self._sigSWeightData :
                # print tagging categories distribution for signal and background
                from RooFitWrappers import ArgSet
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in opposite side tagging category for signal:'
                self._sigSWeightData.table( ArgSet( 'sigOSTagSet', [ tagCatP2VVOS, iTagOS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in opposite side tagging category for background:'
                self._bkgSWeightData.table( ArgSet( 'bkgOSTagSet', [ tagCatP2VVOS, iTagOS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in same side tagging category for signal:'
                self._sigSWeightData.table( ArgSet( 'sigSSTagSet', [ tagCatP2VVSS, iTagSS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in same side tagging category for background:'
                self._bkgSWeightData.table( ArgSet( 'bkgSSTagSet', [ tagCatP2VVSS, iTagSS ] ) ).Print('v')


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
        commonArgs = dict( AmbiguityParameters = ambiguityPars )
        if paramKKMass == 'functions' :
            commonArgs[ 'KKMass' ]        = KKMass
            commonArgs[ 'KKMassBinning' ] = self._KKMassBinning
        if not nominalPdf and not ASParam.startswith('Mag2ReIm') :
            commonArgs[ 'C_SP' ] = CSPValues[0]

        if nominalPdf or amplitudeParam == 'phasesSWaveFrac' :
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = 'deltaPerp' if nominalPdf else ASParam
                                          , AparParameterization = 'phase' if nominalPdf else AparParam
                                          , **commonArgs )

        elif amplitudeParam == 'phases' :
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = 'deltaPerp' if nominalPdf else ASParam, **commonArgs )

        elif amplitudeParam == 'bank' :
            from P2VVParameterizations.DecayAmplitudes import JpsiVBank_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = ASParam, AparParameterization = AparParam
                                          , **commonArgs )

        else :
            raise RuntimeError('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no valid amplitude parameterization specified')

        # B lifetime
        from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        dGammaVar = dict( Name = 'dGamma' )
        if blind : dGammaVar['Blind'] = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
        self._lifetimeParams = LifetimeParams( dGamma = dGammaVar, dMConstraint = True if nominalPdf else constrainDeltaM )
        if ambiguityPars : self._lifetimeParams['dGamma'].setVal( -self._lifetimeParams['dGamma'].getVal() )

        if nominalPdf or eventTimeRes :
            from P2VVParameterizations.TimeResolution import Moriond2012_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution( time = time, sigmat = timeRes, timeResSFConstraint = True )
        else :
            from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution( time = time, timeResSFConstraint = True )

        # CP violation parameters
        if lambdaCPParam == 'ReIm' : 
            from P2VVParameterizations.CPVParams import LambdaCarth_CPParam as CPParam
            ReLambdaCPVar = dict( Name = 'ReLambdaCP' )
            ImLambdaCPVar = dict( Name = 'ImLambdaCP' )
            if blind: ReLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsGoofyMoriond2012', 0.1 )
            if blind: ImLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsPlutoMoriond2012', 0.1 )
            self._lambdaCP = CPParam( ReLambdaCP = ReLambdaCPVar, ImLambdaCP = ImLambdaCPVar )

        else :
            if lambdaCPParam == 'lambPhi' :
                from P2VVParameterizations.CPVParams import LambdaArg_CPParam as CPParam
            else :
                from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam

            phiCPVar = dict( Name = 'phiCP' )
            if blind: phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
            self._lambdaCP = CPParam( phiCP = phiCPVar )
            if ambiguityPars :
                from math import pi
                self._lambdaCP['phiCP'].setVal( pi - self._lambdaCP['phiCP'].getVal() )

        # coefficients for time functions
        from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, self._amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

        # tagging parameters
        if not nominalPdf and self._iTagZeroTrick :
            from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams(  estWTag = estWTagOS, p0 = dict( Name = 'wTagP0' ), p1 = dict( Name = 'wTagP1' )
                                                , p0Constraint = True if nominalPdf else constrainTagging
                                                , p1Constraint = True if nominalPdf else constrainTagging )

            args = dict(  iTag     = tagDecisionOS
                        , dilution = self._taggingParams['dilution']
                        , ADilWTag = self._taggingParams['ADilWTag']
                        , avgCEven = self._taggingParams['avgCEven']
                        , avgCOdd  = self._taggingParams['avgCOdd']
                       )

        elif not nominalPdf and sigTaggingPdf == 'histPdf' :
            raise RuntimeError('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: A histogram tagging PDF can only be used when tagging observables are conditional')

        else :
            # get tagging category parameters dictionary/dictionaries
            tagCatsDictOS = self._tagCatsOS.tagCatsDict()
            if SSTagging: tagCatsDictSS = self._tagCatsSS.tagCatsDict()

            if not nominalPdf and sigTaggingPdf.startswith('tagUntag') :
                # assume products of asymmetries are small and B-Bbar asymmetries are equal for all tagged categories
                if tagCatP2VVOS.numTypes() > 2 or ( SSTagging and tagCatP2VVSS.numTypes() > 2 ) :
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: tagging in signal PDF:\n'\
                        + '    * assuming B-Bbar asymmetries are equal for all tagged categories'

                # provide the same asymmetry for all tagged categories
                from math import sqrt
                asymVal = -self._lambdaCP['C'].getVal()
                asymErr = ( 10. / sqrt( self._sigSWeightData.sumEntries() ) ) if self._sigSWeightData else 0.1
                avgCEvenSum = RealVar( 'avgCEvenSum'    , Title = 'Sum of CP average even coefficients'
                                                        , Value = 1., MinMax = (  0., 2. ) )
                avgCOddSum  = RealVar( 'avgCOddSum'     , Title = 'Sum of CP average odd coefficients'
                                                        , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )
                avgCEvenOS = RealVar( 'avgCEvenOSTagged', Title = 'CP average even coefficients OS tagged categories'
                                                        , Value = 1., MinMax = (  0., 2. ) )
                avgCOddOS  = RealVar( 'avgCOddOSTagged' , Title = 'CP average odd coefficients OS tagged categories'
                                                        , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )
                if not SSTagging :
                    # only opposite side tagging: replace tagging efficiency asymmetry parameters by asymmetry coefficients
                    tagCatsDictOS[ 'AvgCEvenSum' ] = avgCEvenSum
                    tagCatsDictOS[ 'AvgCOddSum' ]  = avgCOddSum
                    for catOS in range( 1, tagCatsDictOS['NumTagCats'] ) :
                        tagCatsDictOS.pop( 'ATagEff%d' % catOS )
                        tagCatsDictOS[ 'AvgCEven%d' % catOS ] = avgCEvenOS
                        tagCatsDictOS[ 'AvgCOdd%d'  % catOS ] = avgCOddOS

                    tagCatsDict = tagCatsDictOS

                else :
                    # both same side tagging and opposite side tagging: build coefficients for SS tagged categories
                    avgCEvenSS = RealVar( 'avgCEvenSSTagged', Title = 'CP average even coefficients SS tagged categories'
                                                            , Value = 1., MinMax = (  0., 2. ) )
                    avgCOddSS  = RealVar( 'avgCOddSSTagged' , Title = 'CP average odd coefficients SS tagged categories'
                                                            , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )
                    avgCEven   = RealVar( 'avgCEvenTagged'  , Title = 'CP average even coefficients OS+SS tagged categories'
                                                            , Value = 1., MinMax = (  0., 2. ) )
                    avgCOdd    = RealVar( 'avgCOddTagged'   , Title = 'CP average odd coefficients OS+SS tagged categories'
                                                            , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )

                    # build dictionary with both opposite side and same side tagging categories parameters
                    tagCatsDict = dict(  NumTagCats0  = tagCatsDictOS['NumTagCats']
                                       , NumTagCats1  = tagCatsDictSS['NumTagCats']
                                       , AvgCEvenSum  = avgCEvenSum
                                       , AvgCOddSum   = avgCOddSum
                                       , Conditionals = tagCatsDictOS['Conditionals'] + tagCatsDictSS['Conditionals']
                                       , Constraints  = tagCatsDictOS['Constraints']  + tagCatsDictSS['Constraints']
                                      )

                    for catOS in range( tagCatsDictOS['NumTagCats'] ) :
                        if catOS > 0 :
                            tagCatsDict[ 'tagCatCoef0_%d'  % catOS ] = tagCatsDictOS[ 'tagCatCoef%d' % catOS ]
                            tagCatsDict[ 'AvgCEven%d-0'    % catOS ] = avgCEvenOS
                            tagCatsDict[ 'AvgCOdd%d-0'     % catOS ] = avgCOddOS
                            tagCatsDict[ 'tagDilution0_%s' % catOS ] = tagCatsDictOS[ 'tagDilution%d' % catOS ]
                            tagCatsDict[ 'ADilWTag0_%s'    % catOS ] = tagCatsDictOS[ 'ADilWTag%d'    % catOS ]

                        for catSS in range( 1, tagCatsDictSS['NumTagCats'] ) :
                            if catOS == 0 :
                                tagCatsDict[ 'tagCatCoef1_%d'  % catSS ] = tagCatsDictSS[ 'tagCatCoef%d' % catSS ]
                                tagCatsDict[ 'AvgCEven0-%s'    % catSS ] = avgCEvenSS
                                tagCatsDict[ 'AvgCOdd0-%s'     % catSS ] = avgCOddSS
                                tagCatsDict[ 'tagDilution1_%s' % catSS ] = tagCatsDictSS[ 'tagDilution%d' % catSS ]
                                tagCatsDict[ 'ADilWTag1_%s'    % catSS ] = tagCatsDictSS[ 'ADilWTag%d'    % catSS ]
                            else :
                                tagCatsDict[ 'AvgCEven%s-%s' % ( catOS, catSS ) ] = avgCEven
                                tagCatsDict[ 'AvgCOdd%s-%s'  % ( catOS, catSS ) ] = avgCOdd

            else :
                # use independent asymmetry for each category
                if not SSTagging :
                    # only opposite side tagging
                    tagCatsDict = tagCatsDictOS

                else :
                    # both opposite side and same side tagging
                    tagCatsDict = dict(  NumTagCats0  = tagCatsDictOS['NumTagCats']
                                       , NumTagCats1  = tagCatsDictSS['NumTagCats']
                                       , Conditionals = tagCatsDictOS['Conditionals'] + tagCatsDictSS['Conditionals']
                                       , Constraints  = tagCatsDictOS['Constraints']  + tagCatsDictSS['Constraints']
                                      )

                    for catOS in range( 1, tagCatsDictOS['NumTagCats'] ) :
                        tagCatsDict[ 'tagCatCoef0_%d'  % catOS ] = tagCatsDictOS[ 'tagcatCoef%d'  % catOS ]
                        tagCatsDict[ 'ATagEff0_%d'     % catOS ] = tagCatsDictOS[ 'ATagEff%d'     % catOS ]
                        tagCatsDict[ 'tagDilution0_%s' % catOS ] = tagCatsDictOS[ 'tagDilution%d' % catOS ]
                        tagCatsDict[ 'ADilWTag0_%s'    % catOS ] = tagCatsDictOS[ 'ADilWTag%d'    % catOS ]

                    for catSS in range( 1, tagCatsDictSS['NumTagCats'] ) :
                        tagCatsDict[ 'tagCatCoef1_%d'  % catSS ] = tagCatsDictSS[ 'tagcatCoef%d'  % catSS ]
                        tagCatsDict[ 'ATagEff1_%d'     % catSS ] = tagCatsDictSS[ 'ATagEff%d'     % catSS ]
                        tagCatsDict[ 'tagDilution1_%s' % catSS ] = tagCatsDictSS[ 'tagDilution%d' % catSS ]
                        tagCatsDict[ 'ADilWTag1_%s'    % catSS ] = tagCatsDictSS[ 'ADilWTag%d'    % catSS ]

                # add production asymmetry and normalization asymmetry to tagging categories dictionary
                tagCatsDict['AProd'] = 0.
                tagCatsDict['ANorm'] = -self._lambdaCP['C'].getVal()

            from P2VVParameterizations.FlavourTagging import CatDilutionsCoefAsyms_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams( **tagCatsDict )

            if not SSTagging :
                args = dict(  tagCat      = tagCatP2VVOS
                            , iTag        = iTagOS
                            , dilutions   = self._taggingParams['dilutions']
                            , ADilWTags   = self._taggingParams['ADilWTags']
                            , avgCEvens   = self._taggingParams['avgCEvens']
                            , avgCOdds    = self._taggingParams['avgCOdds']
                            , tagCatCoefs = self._taggingParams['tagCatCoefs']
                           )
            else :
                args = dict(  tagCat0     = tagCatP2VVOS
                            , tagCat1     = tagCatP2VVSS
                            , iTag0       = iTagOS
                            , iTag1       = iTagSS
                            , dilutions0  = self._taggingParams['dilutions'][0]
                            , dilutions1  = self._taggingParams['dilutions'][1]
                            , ADilWTags0  = self._taggingParams['ADilWTags'][0]
                            , ADilWTags1  = self._taggingParams['ADilWTags'][1]
                            , avgCEvens   = self._taggingParams['avgCEvens']
                            , avgCOdds    = self._taggingParams['avgCOdds']
                            , tagCatCoefs = self._taggingParams['tagCatCoefs']
                           )

        args = dict(  time                   = time
                    , tau                    = self._lifetimeParams['MeanLifetime']
                    , dGamma                 = self._lifetimeParams['dGamma']
                    , dm                     = self._lifetimeParams['dM']
                    , coshCoef               = timeBasisCoefs['cosh']
                    , sinhCoef               = timeBasisCoefs['sinh']
                    , cosCoef                = timeBasisCoefs['cos']
                    , sinCoef                = timeBasisCoefs['sin']
                    , resolutionModel        = self._timeResModel['model']
                    , ConditionalObservables = self._amplitudes.conditionalObservables()\
                                               + self._timeResModel.conditionalObservables()\
                                               + self._taggingParams.conditionalObservables()
                    , ExternalConstraints    = self._lifetimeParams.externalConstraints()\
                                               + self._timeResModel.externalConstraints()\
                                               + self._taggingParams.externalConstraints()
                    , **args
                   )

        # build signal PDF
        from RooFitWrappers import BTagDecay
        if nominalPdf or condTagging : sigPdf = BTagDecay( 'sig_t_angles',             **args )
        else :                         sigPdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

        if angEffMomentsFile :
            # multiply signal PDF with angular efficiency
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with angular efficiency moments from file "%s"'\
                  % angEffMomentsFile

            from P2VVGeneralUtils import RealMomentsBuilder
            moments = RealMomentsBuilder()
#            angMomInds = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
#                          for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            angMomInds = [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not ( nominalPdf or transAngles ) \
                          else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
            moments.appendPYList( self._angleFuncs.angles, angMomInds )
            moments.read(angEffMomentsFile)
            moments.Print()
            sigPdf = moments * sigPdf


        if multiplyByTimeEff in [ 'all', 'signal' ] :
            # multiply signal PDF with time acceptance
            sigPdfTimeAcc = timeAcceptance * sigPdf
        else :
            sigPdfTimeAcc = sigPdf

        self._signalComps += sigPdfTimeAcc
        if nominalPdf or condTagging : self._sig_t_angles = sigPdfTimeAcc
        else :                         self._sig_t_angles_tagCat_iTag = sigPdfTimeAcc


        ###################################################################################################################################
        ## build signal tagging PDF ##
        ##############################

        if nominalPdf or condTagging :
            if not nominalPdf and self._iTagZeroTrick :
                # no implementation for signal tagging PDF with tag = { B, Bbar, Untagged } (yet)
                pass

            else :
                # tagCat = { Untagged, TagCat1, TagCat2, ... }, tag = { B, Bbar }
                sigTaggingData = self._sigRangeData if massRangeBackground else self._sigSWeightData
                if not nominalPdf and sigTaggingData and sigTaggingPdf == 'histPdf' :
                    # create histogram from signal data and use the (fixed) bin coefficients for the PDF
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: creating signal tagging PDFs from %s data'\
                          % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )
                    from RooFitWrappers import HistPdf
                    if not SSTagging :
                        self._sig_tagCat_iTag = HistPdf(  Name = 'sig_tagCat_iTag'
                                                        , Observables = [ tagCatP2VVOS, iTagOS ]
                                                        , Data = sigTaggingData
                                                       )
                    else :
                        self._sig_tagCat_iTag = HistPdf(  Name = 'sig_tagCat_iTag'
                                                        , Observables = [ tagCatP2VVOS, tagCatP2VVSS, iTagOS, iTagSS ]
                                                        , Data = sigTaggingData
                                                       )
                    self._signalComps += self._sig_tagCat_iTag

                else :
                    # use a PDF with variable bin coefficients
                    if nominalPdf or sigTaggingPdf.startswith('tagUntag')\
                            or ( tagCatP2VVOS.numTypes() == 2 and ( not SSTagging or tagCatP2VVSS.numTypes() == 2 ) ) :
                        # assume B-Bbar asymmetry is equal for all tagged categories
                        if tagCatP2VVOS.numTypes() > 2 or (SSTagging and tagCatP2VVSS.numTypes() > 2 ) :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: signal tagging PDF:\n'\
                                + '    * assuming B-Bbar asymmetries are equal for all tagged categories'
                        else :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a binned signal tagging PDF'

                        from P2VVParameterizations.FlavourTagging import TagUntag_BinnedTaggingPdf as TaggingPdf

                    else :
                        # create independent tagging bin coefficients
                        if sigTaggingData :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining signal tagging coefficients from %s'\
                                  % ( 'B mass signal data' if massRangeBackground else 'B mass S-weight data' )
                        else :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: WARNING:'\
                                + ' no signal data available to determine signal tagging coefficients:\n'\
                                + '    * assuming absence of B-Bbar asymmetries'

                        from P2VVParameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

                    # build PDF
                    self._sigTaggingPdf = TaggingPdf(  'tagCat_iTag', tagCatP2VVOS, iTagOS
                                                  , NamePrefix    = 'sig'
                                                  , TagCatCoefs   = sigPdf.tagCatCoefs()
                                                  , TaggedCatName = 'TagCat' if tagCatP2VVOS.numTypes() > 2 else 'Tagged'
                                                  , Data          = sigTaggingData
                                                  , RelativeCoefs = False
                                                 )

                    self._sig_tagCat_iTag = self._sigTaggingPdf.pdf()
                    self._signalComps += self._sig_tagCat_iTag


        ###################################################################################################################################
        ## build PDFs for estimated wrong-tag probabilities ##
        ######################################################

        if self._data and self._sigSWeightData and self._bkgSWeightData and makePlots :
            tempTagCatOS = tagCatOS if not nominalPdf and self._iTagZeroTrick else tagCatP2VVOS
            tempTagCatSS = None     if not nominalPdf and self._iTagZeroTrick else tagCatP2VVSS

            # build PDFs for estimated wrong-tag probabilities
            from RooFitWrappers import HistPdf
            self._estWTagDataOS = self._sigSWeightData.reduce( '%s > 0' % tempTagCatOS.GetName() )
            if tempTagCatSS :
                self._estWTagDataSS = self._sigSWeightData.reduce( '%s > 0' % tempTagCatSS.GetName() )
            else :
                self._estWTagDataSS = self._sigSWeightData

            self._sig_bkg_estWTagOS = HistPdf(  Name = 'sig_bkg_estWTagOS'
                                              , Observables = [ estWTagOS ]
                                              , Binning = { estWTagOS : numEstWTagBins }
                                              , Data = self._estWTagDataOS
                                             )

            self._sig_bkg_estWTagSS = HistPdf(  Name = 'sig_bkg_estWTagSS'
                                              , Observables = [ estWTagSS ]
                                              , Binning = { estWTagSS : numEstWTagBins }
                                              , Data = self._estWTagDataSS
                                             )

            # get normalization correction for tagged events
            untagFracOS    = self._data.table(tempTagCatOS).getFrac('Untagged')
            untagFracSigOS = self._sigSWeightData.table(tempTagCatOS).getFrac('Untagged')
            untagFracBkgOS = self._bkgSWeightData.table(tempTagCatOS).getFrac('Untagged')
            if tempTagCatSS :
                untagFracSS    = self._data.table(tempTagCatSS).getFrac('Untagged')
                untagFracSigSS = self._sigSWeightData.table(tempTagCatSS).getFrac('Untagged')
                untagFracBkgSS = self._bkgSWeightData.table(tempTagCatSS).getFrac('Untagged')
            else :
                untagFracSS    = 1.
                untagFracSigSS = 1.
                untagFracBkgSS = 1.

            # plot estimated wrong-tag probabilities for signal and for background
            self._estWTagCanvOS = TCanvas( 'estWTagCanvOS', 'Estimated wrong-tag probability OS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvOS.pads( 1, 1 ) if SFit else self._estWTagCanvOS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ self._sigSWeightData if SFit else self._data, self._sigSWeightData, self._bkgSWeightData ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigOS if SFit else untagFracOS ), 1. - untagFracSigOS, 1. - untagFracBkgOS ]
                        ) :
                plot(  pad, estWTagOS, data, self._sig_bkg_estWTagOS
                     , frameOpts  = dict( Bins = nBins, Title = estWTagOS.GetTitle() + plotTitle, Range = ( 0., 0.499999 ) )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                )
                     , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2, Normalization = norm                           )
                    )

            self._estWTagCanvSS = TCanvas( 'estWTagCanvSS', 'Estimated wrong-tag probability SS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvSS.pads( 1, 1 ) if SFit else self._estWTagCanvSS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ self._sigSWeightData if SFit else self._data, self._sigSWeightData, self._bkgSWeightData ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigSS if SFit else untagFracSS ), 1. - untagFracSigSS, 1. - untagFracBkgSS ]
                        ) :
                plot(  pad, estWTagSS, data, self._sig_bkg_estWTagSS
                     , frameOpts  = dict( Bins = nBins, Title = estWTagSS.GetTitle() + plotTitle, Range = ( 0., 0.499999 ) )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                )
                     , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2, Normalization = norm                           )
                    )


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        if not SFit or makePlots :
            from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
            self._backgroundTime = BackgroundTime( Name = 'bkg_t', time = time, resolutionModel = self._timeResModel['model'] )
            self._bkg_t = self._backgroundTime.pdf()

            if multiplyByTimeEff in [ 'all', 'background' ] :
                # multiply background time PDF with time acceptance
                self._bkg_t = timeAcceptance * self._bkg_t

            if not SFit : self._backgroundComps += self._bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit or makePlots :
            bkgAngleData = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
            if bkgAngleData and bkgAnglePdf == 'histPdf' :
                # create a histogram from background data and use it for background angular PDF
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining angular shape of background from %s data'\
                      % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )

                from RooFitWrappers import HistPdf
                nBins = [ 5, 7, 9 ] if nominalPdf else numAngleBins
                self._bkg_angles = HistPdf(  Name = 'bkg_angles'
                                           , Observables = angles
                                           , Binning = { cpsi : nBins[0], ctheta : nBins[1], phi : nBins[2] }
                                           , Data = bkgAngleData
                                          )

            else :
                # create a binned PDF for background angular shape
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a binned angular PDF for background'

                # define angular bins
                from math import pi
                from array import array
                if nominalPdf or transAngles :
                    nBins = [ 5, 7, 9 ]
                    cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
                    cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
                    phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )
                else :
                    nBins = [ 5, 40, 5 ]
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
                            'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: background angular PDF coefficient \"%s\" has negative value: %f'\
                            % (coef.GetName(), value)
                        coef.setVal(value)

            if not SFit : self._backgroundComps += self._bkg_angles

            if self._bkgRangeData and self._bkgSWeightData and makePlots :
                # plot background angles with S-weights
                self._bkgAnglesSWeightCanv = TCanvas( 'bkgAnglesSWeightCanv', 'Background Decay Angles with S-Weights' )
                for ( pad, obs, data, bins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesSWeightCanv.pads( 3, 2 )
                             , angles
                             , 3 * ( self._bkgSWeightData, )
                             , nBins
                             , [ angle.GetTitle() for angle in angles ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
                        )

                # plot background angles with side bands
                self._bkgAnglesSideBandCanv = TCanvas( 'bkgAnglesSideBandCanv', 'Background Decay Angles with Side Bands' )
                for ( pad, obs, data, bins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesSideBandCanv.pads( 3, 2 )
                             , angles
                             , 3 * ( self._bkgRangeData, )
                             , nBins
                             , [ angle.GetTitle() for angle in angles ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
                        )

                # plot 2-D angular distributions
                from P2VVGeneralUtils import _P2VVPlotStash
                for angle0, angle1, data, canv, padNr in [  ( 1, 0, self._bkgSWeightData, self._bkgAnglesSWeightCanv,  4 )
                                                          , ( 2, 0, self._bkgSWeightData, self._bkgAnglesSWeightCanv,  5 )
                                                          , ( 1, 2, self._bkgSWeightData, self._bkgAnglesSWeightCanv,  6 )
                                                          , ( 1, 0, self._bkgRangeData,   self._bkgAnglesSideBandCanv, 4 )
                                                          , ( 2, 0, self._bkgRangeData,   self._bkgAnglesSideBandCanv, 5 )
                                                          , ( 1, 2, self._bkgRangeData,   self._bkgAnglesSideBandCanv, 6 )
                                                         ] :
                    bkgAngHist = data.createHistogram( angles[angle0]._var, angles[angle1]._var, nBins[angle0], nBins[angle1] )
                    _P2VVPlotStash.append(bkgAngHist)
                    bkgAngHist.SetStats(False)
                    bkgAngHist.SetTitle( '%s vs. %s' % ( angleNames[angle0][1], angleNames[angle1][1] ) )
                    bkgAngHist.SetMinimum(0.)
                    bkgAngHist.GetXaxis().SetTitle( angleNames[angle0][1] )
                    bkgAngHist.GetYaxis().SetTitle( angleNames[angle1][1] )
                    bkgAngHist.GetXaxis().SetLabelOffset(0.01)
                    bkgAngHist.GetYaxis().SetLabelOffset(0.008)
                    bkgAngHist.GetXaxis().SetTitleOffset(1.8)
                    bkgAngHist.GetYaxis().SetTitleOffset(1.8)
                    bkgAngPad = canv.cd(padNr)
                    bkgAngPad.SetLeftMargin(0.08)
                    bkgAngPad.SetRightMargin(0.05)
                    bkgAngPad.SetBottomMargin(0.05)
                    bkgAngPad.SetTopMargin(0.05)
                    bkgAngHist.Draw('lego2')


        ###################################################################################################################################
        ## build background tagging PDF ##
        ##################################

        if not SFit :
            if not nominalPdf and self._iTagZeroTrick :
                # no implementation for background tagging PDF with tag = { B, Bbar, Untagged } (yet)
                pass

            else :
                # tagCat = { Untagged, TagCat1, TagCat2, ... }, tag = { B, Bbar }
                bkgTaggingData = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
                if not nominalPdf and bkgTaggingData and bkgTaggingPdf == 'histPdf' :
                    # create histogram from background data and use the (fixed) bin coefficients for the PDF
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: creating background tagging PDFs from %s data'\
                          % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )
                    from RooFitWrappers import HistPdf
                    self._bkg_tagCat_iTag = HistPdf( Name = 'bkg_tagCat_iTag', Observables = [tagCatP2VVOS, iTagOS], Data = bkgTaggingData )
                    self._backgroundComps += self._bkg_tagCat_iTag

                else :
                    # use a PDF with variable bin coefficients
                    if nominalPdf or bkgTaggingPdf.startswith('tagUntag') or tagCatP2VVOS.numTypes() == 2 :
                        # couple background tagging category coefficients to signal tagging category coefficients
                        # and assume B-Bbar asymmetry is equal for all tagged categories
                        if tagCatP2VVOS.numTypes() > 2 :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: background tagging PDF:\n'\
                                + '    * assuming signal tagging category distribution for tagged categories\n'\
                                + '    * assuming B-Bbar asymmetries are equal for all tagged categories'
                        else :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a binned background tagging PDF'

                        from P2VVParameterizations.FlavourTagging import TagUntag_BinnedTaggingPdf as TaggingPdf

                    else :
                        # create independent tagging bin coefficients
                        if bkgTaggingData :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining background tagging coefficients from %s'\
                                  % ( 'B mass side band data' if massRangeBackground else 'B mass S-weight data' )
                        else :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: WARNING:'\
                                + ' no background data available to determine background tagging coefficients:\n'\
                                + ( '    * setting tagging category coefficients to signal values\n'\
                                  if bkgTaggingPdf.endswith('Relative') else '' )\
                                + '    * assuming absence of B-Bbar asymmetries'

                        from P2VVParameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

                    # build PDF
                    self._bkgTaggingPdf = TaggingPdf(  'tagCat_iTag', tagCatP2VVOS, iTagOS
                                                     , NamePrefix    = 'bkg'
                                                     , TagCatCoefs   = sigPdf.tagCatCoefs()\
                                                                       if bkgTaggingPdf.endswith('Relative') else None
                                                     , TaggedCatName = 'TagCat' if tagCatP2VVOS.numTypes() > 2 else 'Tagged'
                                                     , Data          = bkgTaggingData
                                                    )

                    self._bkg_tagCat_iTag = self._bkgTaggingPdf.pdf()
                    self._backgroundComps += self._bkg_tagCat_iTag


        ###################################################################################################################################
        ## build full PDF ##
        ####################

        from RooFitWrappers import buildPdf
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: requesting %s PDF for observables [%s]'\
              % ( 'signal' if SFit else 'signal + background', ', '.join( str(obs) for obs in obsSetP2VV ) )
        if SFit :
            self._fullPdf = buildPdf( [ self._signalComps ], Observables = obsSetP2VV, Name = 'Jpsiphi' )
        else :
            self._fullPdf = buildPdf( [ self._signalComps, self._backgroundComps ], Observables = obsSetP2VV, Name = 'Jpsiphi' )


        ###################################################################################################################################
        ## split PDF for different KK mass bins ##
        ##########################################

        if paramKKMass == 'simultaneous' :
            # create split category
            from RooFitWrappers import BinningCategory
            self._KKMassCat = BinningCategory( 'KKMassCat'
                                              , Observable = KKMass
                                              , Binning = self._KKMassBinning
                                              , Fundamental = True
                                              , Data = [  self._data
                                                        , self._sigSWeightData
                                                        , self._bkgSWeightData
                                                        , self._sigRangeData
                                                        , self._bkgRangeData
                                                       ]
                                              , CatTypeName = 'bin'
                                             )

            # specify parameters that are different in simultaneous categories
            splitParams = [ ]
            for amp in self._amplitudes.parameters() :
                if any( name in amp.GetName() for name in [ 'AS', 'A_S', 'fS', 'f_S' ] ) : splitParams.append(amp)
            if nominalPdf or condTagging :
                for par in self._sigTaggingPdf.parameters() :
                    if not par.isConstant() : splitParams.append(par)
            if not SFit :
                splitParams.append( self._signalComps.getYield() )
                splitParams.append( self._backgroundComps.getYield() )
                for par in self._bkgTaggingPdf.parameters() :
                    if not par.isConstant() : splitParams.append(par)

            # build simultaneous PDF
            from RooFitWrappers import SimultaneousPdf
            self._simulPdf = SimultaneousPdf(  self._fullPdf.GetName() + '_KKMassBins'
                                             , MasterPdf       = self._fullPdf
                                             , SplitCategory   = self._KKMassCat
                                             , SplitParameters = splitParams
                                            )

            # set values for splitted parameters
            splitCatIter = self._simulPdf.indexCat().typeIterator()
            splitCatState = splitCatIter.Next()
            while splitCatState :
                splitCatPars = self._simulPdf.getPdf( splitCatState.GetName() ).getVariables()

                if amplitudeParam == 'bank' and ASParam != 'ReIm' :
                    ASOddMag2  = splitCatPars.find( 'ASOddMag2_'  + splitCatState.GetName() )
                    ASOddPhase = splitCatPars.find( 'ASOddPhase_' + splitCatState.GetName() )

                    ASOddMag2.setVal( SWaveAmpVals[0][ splitCatState.getVal() ] )
                    ASOddMag2.setMax(5.)
                    ASOddPhase.setVal( SWaveAmpVals[1][ splitCatState.getVal() ] )

                elif amplitudeParam == 'phases' and ASParam in [ 'ReIm', 'Mag2ReIm' ] :
                    if ASParam == 'Mag2ReIm' :
                        ASMag2 = splitCatPars.find( 'ASMag2_' + splitCatState.GetName() )
                    ReAS = splitCatPars.find( 'ReAS_' + splitCatState.GetName() )
                    ImAS = splitCatPars.find( 'ImAS_' + splitCatState.GetName() )

                    if ASParam == 'Mag2ReIm' :
                        ASMag2.setVal( SWaveAmpVals[0][ splitCatState.getVal() ]**2 + SWaveAmpVals[1][ splitCatState.getVal() ]**2 )
                    ReAS.setVal( SWaveAmpVals[0][ splitCatState.getVal() ] )
                    ImAS.setVal( SWaveAmpVals[1][ splitCatState.getVal() ] )

                if not SFit :
                    sigYield = splitCatPars.find( 'N_signal_' + splitCatState.GetName() )
                    bkgYield = splitCatPars.find( 'N_bkg_'    + splitCatState.GetName() )

                    from math import sqrt
                    nSigEvBin = self._signalKKMass['numEventsBins'][ splitCatState.getVal() ]
                    nBkgEvBin = self._backgroundKKMass['numEventsBins'][ splitCatState.getVal() ]
                    sigYield.setVal( nSigEvBin )
                    sigYield.setError( sqrt(nSigEvBin) )
                    sigYield.setMin(0.)
                    sigYield.setMax( nSigEvBin + nBkgEvBin )
                    bkgYield.setVal( nBkgEvBin )
                    bkgYield.setError( sqrt(nBkgEvBin) )
                    bkgYield.setMin(0.)
                    bkgYield.setMax( nSigEvBin + nBkgEvBin )

                splitCatState = splitCatIter.Next()


        assert not pdfConfig, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: superfluous arguments found: %s' % pdfConfig
        PdfBuilder.__init__( self, self._simulPdf if paramKKMass == 'simultaneous' else self._fullPdf, observables, { } )
