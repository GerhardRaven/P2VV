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
        from math import pi

        # job parameters
        self['dataSample'] = ''              # '' / 'Summer2011'
        self['selection']  = 'HLT1Unbiased'  # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'timeEffFit'
        self['makePlots']  = True
        self['SFit']       = False
        self['blind']      = False
        self['nominalPdf'] = True    # nominal PDF option does not work at the moment

        self['numEvents'] = 30000

        self['nTupleName']     = 'DecayTree'
        self['nTupleFile']     = 'Bs2JpsiPhi_ntupleB_for_fitting_20120203.root'
        self['nominalDataSet'] = False

        self['timeEffHistFile']      = 'BuBdBdJPsiKsBsLambdab0_HltPropertimeAcceptance_20120504.root'
        self['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1UB_40bins'
        self['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_Data_Hlt2BHlt1ExclB_40bins'
        self['angEffMomentsFile']    = 'effMoments'

        # fit options
        self['fitOptions'] = dict( NumCPU = 1, Optimize = 1, Timer = True, Save = True )

        # PDF parameters
        self['tagCatsOS'] = [ ]
        self['tagCatsSS'] = [ ]

        # PDF options
        self['transversityAngles']   = False
        self['angularRanges']        = dict( cpsi = [ ( -1., +1. ) ], ctheta = [ ( -1., +1. ) ], phi = [ ( -pi, +pi ) ] )
        self['sigMassModel']         = 'doubleGauss'       # '' / 'doubleGauss' / 'box'
        self['bkgMassModel']         = 'exponential'       # '' / 'exponential' / 'linear'
        self['bkgAnglePdf']          = 'binned'            # '' / 'histPdf' / 'binned' / 'basis' / 'hybrid'
        self['sigTaggingPdf']        = 'tagUntag'          # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        self['bkgTaggingPdf']        = 'tagUntagRelative'  # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        self['multiplyByTagPdf']     = False
        self['multiplyByTimeEff']    = ''                  # '' / 'all' / 'signal'
        self['timeEffType']          = 'HLT1Unbiased'      # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'fit'
        self['multiplyByAngEff']     = ''                  # '' / 'basis012' / 'basis0123' / 'basisSig3' / 'basisSig6' / 'weights'
        self['parameterizeKKMass']   = ''                  # '' / 'functions' / 'simultaneous'
        self['ambiguityParameters']  = False
        self['SWeightsType']         = ''                  # '' / 'simultaneous' / 'simultaneousFixed' / 'simultaneousFreeBkg'
        self['KKMassBinBounds']      = [ 1020. - 12., 1020. + 12. ]
        self['SWaveAmplitudeValues'] = (  [ ], [ ] )
        self['CSPValues']            = [ 0.4976 ]
        self['lifetimeRange']        = ( 0.3, 14. )

        self['sameSideTagging']    = True
        self['conditionalTagging'] = False
        self['continuousEstWTag']  = False
        self['numEstWTagBins']     = 100
        self['constrainTagging']   = 'constrain'  # '' / 'constrain' / 'fixed'

        self['iTagZeroTrick'] = False
        self['iTagStates'] = { }        # { } / { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }

        self['timeResType']           = 'event'      # '' / 'event' / 'eventNoMean' / 'eventConstMean' / '3Gauss'
        self['numTimeResBins']        = 20
        self['constrainTimeResScale'] = 'constrain'  # '' / 'constrain' / 'fixed'

        self['signalFraction'] = 0.67
        self['massRangeBackground'] = False

        self['amplitudeParam'] = 'phasesSWaveFrac'  # 'phases' / 'phasesSWaveFrac' / 'ReIm' / 'bank'
        self['ASParam']        = 'deltaPerp'        # 'delta0' / 'deltaPerp' / 'ReIm' / 'Mag2ReIm' / 'Mag2ReImPerp'
        self['AparParam']      = 'cos'              # 'phase' / 'ReIm' / 'Mag2ReIm' / 'cos' / 'real'

        self['constrainDeltaM'] = 'constrain'  # '' / 'constrain' / 'fixed'

        self['lambdaCPParam'] = 'lambSqPhi'         # 'ReIm' / 'lambSqPhi' / 'lambPhi' / 'lambPhi_CPVDecay' / 'lambPhiRel_CPVDecay'

        self['angleNames'] = (  ( 'trcospsi',   'cos(#psi_{tr})'   )
                              , ( 'trcostheta', 'cos(#theta_{tr})' )
                              , ( 'trphi',      '#phi_{tr}'        )
                             )

        self['numBMassBins'] = [ 40, 20, 20 ]
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

        from math import sqrt, pi

        # copy configuration arguments
        pdfConfig = kwargs.copy()

        # job parameters
        dataSample = pdfConfig.pop('dataSample')
        selection  = pdfConfig.pop('selection')
        SFit       = pdfConfig.pop('SFit')
        blind      = pdfConfig.pop('blind')
        nominalPdf = pdfConfig.pop('nominalPdf')
        makePlots  = pdfConfig.pop('makePlots')

        numEvents = pdfConfig.pop('numEvents')

        nTupleName           = pdfConfig.pop('nTupleName')
        nTupleFile           = pdfConfig.pop('nTupleFile')
        nominalDataSet       = pdfConfig.pop('nominalDataSet')
        angEffMomentsFile    = pdfConfig.pop('angEffMomentsFile')
        timeEffHistFile      = pdfConfig.pop('timeEffHistFile')
        timeEffHistUBName    = pdfConfig.pop('timeEffHistUBName')
        timeEffHistExclBName = pdfConfig.pop('timeEffHistExclBName')

        fitOpts = pdfConfig.pop('fitOptions')

        angleNames   = pdfConfig.pop('angleNames')
        numTimeBins  = pdfConfig.pop('numTimeBins')
        numAngleBins = pdfConfig.pop('numAngleBins')

        # PDF parameters
        parameters = pdfConfig.pop('parameters')
        tagCatsOS  = pdfConfig.pop('tagCatsOS')
        tagCatsSS  = pdfConfig.pop('tagCatsSS')

        # PDF options
        transAngles       = pdfConfig.pop('transversityAngles')
        angRanges         = pdfConfig.pop('angularRanges')
        sigMassModel      = pdfConfig.pop('sigMassModel')
        bkgMassModel      = pdfConfig.pop('bkgMassModel')
        bkgAnglePdf       = pdfConfig.pop('bkgAnglePdf')
        sigTaggingPdf     = pdfConfig.pop('sigTaggingPdf')
        bkgTaggingPdf     = pdfConfig.pop('bkgTaggingPdf')
        multiplyByTagPdf  = pdfConfig.pop('multiplyByTagPdf')
        multiplyByTimeEff = pdfConfig.pop('multiplyByTimeEff')
        timeEffType       = pdfConfig.pop('timeEffType')
        multiplyByAngEff  = pdfConfig.pop('multiplyByAngEff')
        paramKKMass       = pdfConfig.pop('parameterizeKKMass')
        numBMassBins      = pdfConfig.pop('numBMassBins')
        ambiguityPars     = pdfConfig.pop('ambiguityParameters')
        SWeightsType      = pdfConfig.pop('SWeightsType')
        KKMassBinBounds   = pdfConfig.pop('KKMassBinBounds')
        SWaveAmpVals      = pdfConfig.pop('SWaveAmplitudeValues')
        CSPValues         = pdfConfig.pop('CSPValues')
        lifetimeRange     = pdfConfig.pop('lifetimeRange')

        assert lifetimeRange[0] > -5. and lifetimeRange[1] < 50. and lifetimeRange[1] > lifetimeRange[0]

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

        timeResType     = pdfConfig.pop('timeResType')
        numTimeResBins  = pdfConfig.pop('numTimeResBins')
        constrTResScale = pdfConfig.pop('constrainTimeResScale')

        sigFrac = pdfConfig.pop('signalFraction')
        massRangeBackground = pdfConfig.pop('massRangeBackground')

        amplitudeParam = pdfConfig.pop('amplitudeParam')
        ASParam        = pdfConfig.pop('ASParam')
        AparParam      = pdfConfig.pop('AparParam')

        if not angRanges : angRanges = dict( cpsi = [ ( -1., +1. ) ], ctheta = [ ( -1., +1. ) ], phi = [ ( -pi, +pi ) ] )

        if not paramKKMass :
            if not KKMassBinBounds : KKMassBinBounds = [ 1020. - 12., 1020. + 12. ]
        elif ambiguityPars :
            if nominalPdf or ( amplitudeParam == 'phasesSWaveFrac' and ASParam == 'deltaPerp' )\
                    or ( amplitudeParam == 'bank' and ASParam != 'ReIm' ) :
                for phaseIter, phase in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][phaseIter] = ( pi - phase[0], phase[1] )
            elif amplitudeParam == 'phases' and ASParam in [ 'Mag2ReIm', 'ReIm' ] :
                for ImIter, Im in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][ImIter] = ( -Im[0], Im[1] )

        constrainDeltaM = pdfConfig.pop('constrainDeltaM')

        lambdaCPParam = pdfConfig.pop('lambdaCPParam')


        ###################################################################################################################################
        ## imports and function definitions ##
        ######################################

        # python garbage collector
        import gc

        # ROOT imports
        from ROOT import TFile

        # RooFit infinity
        from ROOT import RooNumber
        RooInf = RooNumber.infinity()

        # RooObject wrappers
        from RooFitWrappers import RooObject, ConstVar, RealVar, Category, FormulaVar
        ws = RooObject().ws()

        if makePlots :
            # import plotting tools
            from P2VVGeneralUtils import plot
            from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

        # function for finding a splitted parameter in a simultaneous PDF
        from P2VVGeneralUtils import getSplitPar


        ###################################################################################################################################
        ## create variables (except for tagging category) ##
        ####################################################

        # angular functions
        if not nominalPdf and transAngles :
            helAngs = False
            from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs

        else :
            helAngs = True
            from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs

        self._angleFuncs = AngleFuncs(  cpsi   = dict(Name = 'helcosthetaK' if helAngs else 'trcospsi',   MinMax = angRanges['cpsi'  ][0])
                                      , ctheta = dict(Name = 'helcosthetaL' if helAngs else 'trcostheta', MinMax = angRanges['ctheta'][0])
                                      , phi    = dict(Name = 'helphi'       if helAngs else 'trphi',      MinMax = angRanges['phi'   ][0])
                                     )

        # variables in PDF (except for tagging category)
        time = RealVar( 'time', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = lifetimeRange
                       , Ranges = dict( Bulk = ( None, 5. ) ), nBins = numTimeBins )
        timeRes = RealVar(  'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, Value = 0.10, MinMax = (0.0001, 0.12) # > 0.0075
                          , nBins = numTimeResBins )
        timeRes.setBins( numTimeResBins, 'cache' )

        cpsi   = self._angleFuncs.angles['cpsi']
        ctheta = self._angleFuncs.angles['ctheta']
        phi    = self._angleFuncs.angles['phi']

        for ran in angRanges['cpsi'  ][ 1 : ] : cpsi.setRange(   ran[0], ran[ 1 : 3 ] )
        for ran in angRanges['ctheta'][ 1 : ] : ctheta.setRange( ran[0], ran[ 1 : 3 ] )
        for ran in angRanges['phi'   ][ 1 : ] : phi.setRange(    ran[0], ran[ 1 : 3 ] )

        if nominalPdf or not self._iTagZeroTrick :
            iTagOS = Category( 'iTagOS', Title = 'Initial state flavour tag opposite side', Observable = True, States = iTagStates )
            iTagSS = Category( 'iTagSS', Title = 'Initial state flavour tag same side',     Observable = True, States = iTagStates )
        estWTagComb = RealVar( 'tagomega',    Title = 'Estimated wrong tag probability OS/SS combination', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ) )
        estWTagOS   = RealVar( 'tagomega_os', Title = 'Estimated wrong tag probability opposite side', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ) )
        estWTagSS   = RealVar( 'tagomega_ss', Title = 'Estimated wrong tag probability same side', Observable = True
                              , Value = 0.25, MinMax = ( 0., 0.50001 ) )

        BMass = RealVar( 'mass',  Title = 'm(J/#psiKK)', Unit = 'MeV', Observable = True
                        , Value = 5368., MinMax = ( 5200., 5550. ), nBins = numBMassBins[0] + numBMassBins[1] + numBMassBins[2]
                        ,  Ranges = dict(  LeftSideBand  = ( 5200., 5320. )
                                         , Signal        = ( 5320., 5420. )
                                         , RightSideBand = ( 5420., 5550. )
                                        )
                       )

        angles = [ cpsi, ctheta, phi ]
        obsSetP2VV = [ time ] + angles
        if not self._iTagZeroTrick and ( nominalPdf or multiplyByTagPdf or not condTagging ) :
            obsSetP2VV.append(iTagOS)
            if SSTagging : obsSetP2VV.append(iTagSS)
        if not SFit :
            obsSetP2VV.append(BMass)

        # ntuple variables
        mumuMass = RealVar( 'mdau1', Title = 'm(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. )
                           , nBins =  51 )
        KKMass = RealVar( 'mdau2', Title = 'm(KK)', Unit = 'MeV', Observable = True, MinMax = ( KKMassBinBounds[0], KKMassBinBounds[-1] )
                         , nBins =  51 )
        #if paramKKMass == 'functions' : obsSetP2VV.append(KKMass)

        tagDecisionComb = Category( 'tagdecision',    Title = 'Tag decision OS/SS combination', Observable = True, States = iTagStatesDecision )
        tagDecisionOS   = Category( 'tagdecision_os', Title = 'Tag decision opposite side',     Observable = True, States = iTagStatesDecision )
        tagDecisionSS   = Category( 'tagdecision_ss', Title = 'Tag decision same side',         Observable = True, States = iTagStatesDecision )
        tagCatOS = Category( 'tagcat_os',   Title = 'Tagging category opposite side', Observable = True
                            , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                           )

        hlt1ExclBName = 'hlt1_excl_biased_dec'
        #hlt1ExclBName = 'hlt1_excl_biased'
        #hlt1ExclBName = 'triggerDecisionBiasedExcl'
        #hlt1ExclBName = 'trigger_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
        hlt1BName     = 'hlt1_biased'
        hlt1UBName    = 'hlt1_unbiased_dec'
        #hlt1UBName    = 'hlt1_unbiased'
        #hlt1UBName    = 'triggerDecisionUnbiased'
        #hlt1UBName    = 'trigger_Hlt1DiMuon_Hlt2DiMuonDetached'
        hlt2BName     = 'hlt2_biased'
        hlt2UBName    = 'hlt2_unbiased'

        sel       = Category( 'sel',         Observable = True, States = { 'sel'   : 1, 'notSel'   : 0 } )
        selA      = Category( 'selA',        Observable = True, States = { 'sel'   : 1, 'notSel'   : 0 } )
        selB      = Category( 'selB',        Observable = True, States = { 'sel'   : 1, 'notSel'   : 0 } )
        hlt1ExclB = Category( hlt1ExclBName, Observable = True, States = { 'exclB' : 1, 'notExclB' : 0 } )
        hlt1B     = Category( hlt1BName,     Observable = True, States = { 'B'     : 1, 'notB'     : 0 } )
        hlt1UB    = Category( hlt1UBName,    Observable = True, States = { 'UB'    : 1, 'notUB'    : 0 } )
        hlt2B     = Category( hlt2BName,     Observable = True, States = { 'B'     : 1, 'notB'     : 0 } )
        hlt2UB    = Category( hlt2UBName,    Observable = True, States = { 'UB'    : 1, 'notUB'    : 0 } )

        muPlusTrackChi2 = RealVar( 'muplus_track_chi2ndof',  Title = 'mu+ track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
        muMinTrackChi2  = RealVar( 'muminus_track_chi2ndof', Title = 'mu- track chi^2/#dof', Observable = True, MinMax = ( 0., 4. ) )
        KPlusTrackChi2  = RealVar( 'Kplus_track_chi2ndof',   Title = 'K+ track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )
        KMinTrackChi2   = RealVar( 'Kminus_track_chi2ndof',  Title = 'K- track chi^2/#dof',  Observable = True, MinMax = ( 0., 4. ) )

        observables = dict(  time             = time
                           , cpsi             = angles[0]
                           , ctheta           = angles[1]
                           , phi              = angles[2]
                           , iTagOS           = iTagOS if nominalPdf or not self._iTagZeroTrick else tagDecisionOS
                           , iTagSS           = iTagSS if nominalPdf or not self._iTagZeroTrick else tagDecisionSS
                           , tagDecisionComb  = tagDecisionComb
                           , tagDecisionOS    = tagDecisionOS
                           , estWTagComb      = estWTagComb
                           , estWTagOS        = estWTagOS
                           , estWTagSS        = estWTagSS
                           , tagCatOS         = tagCatOS
                           , BMass            = BMass
                           , mumuMass         = mumuMass
                           , KKMass           = KKMass
                           , timeRes          = timeRes
                           , sel              = sel
                           , selA             = selA
                           , selB             = selB
                           , hlt1ExclB        = hlt1ExclB
                           , hlt1B            = hlt1B
                           , hlt1UB           = hlt1UB
                           , hlt2B            = hlt2B
                           , hlt2UB           = hlt2UB
                           , muPlusTrackChi2  = muPlusTrackChi2
                           , muMinTrackChi2   = muMinTrackChi2
                           , KPlusTrackChi2   = KPlusTrackChi2
                           , KMinTrackChi2    = KMinTrackChi2
                          )

        obsSetNTuple = [ time ] + angles +  [ BMass, KKMass, timeRes, tagDecisionOS, tagDecisionSS, estWTagOS, estWTagSS ]
        if makePlots : obsSetNTuple += [ mumuMass ]
        if selection in [ 'timeEffFit', 'paper2012' ] or timeEffType in [ 'fit', 'paper2012' ] : obsSetNTuple += [ hlt1ExclB ]
        if selection == 'timeEffFit'                  or timeEffType == 'fit'                  : obsSetNTuple += [ hlt2B, hlt2UB ]
        if nominalDataSet : obsSetNTuple += [  sel, selA, selB, muPlusTrackChi2, muMinTrackChi2, KPlusTrackChi2, KMinTrackChi2
                                             , tagDecisionComb, estWTagComb ]


        ###################################################################################################################################
        ## read data ##
        ###############

        self._dataSets = { }
        ntupleCuts = ''
        if nTupleFile :
            if dataSample and type(dataSample) == tuple :
                if dataSample[0] and dataSample[1] : dataSampleCuts = 'runNumber > %d && runNumber < %d' % dataSample
                elif dataSample[0]                 : dataSampleCuts = 'runNumber > %d' % dataSample[0]
                elif dataSample[1]                 : dataSampleCuts = 'runNumber < %d' % dataSample[1]
                else                               : dataSampleCuts = 'runNumber > 0'

                if len(dataSample) > 2 and dataSample[2] : ntupleCuts += dataSample[2] + ' && '

            elif dataSample == 'Summer2011' :
                dataSampleCuts = 'runNumber > 87219 && runNumber < 94386'

            elif dataSample and type(dataSample) == str :
                dataSampleCuts = dataSample

            else :
                dataSampleCuts = 'runNumber > 0'
            trackChi2Cuts = ' && '.join( '%s < %f' % ( trackChi2, trackChi2.getMax() ) for trackChi2 in\
                                                      [ muPlusTrackChi2, muMinTrackChi2, KPlusTrackChi2, KMinTrackChi2 ] )

            if selection == 'HLT1Unbiased' :
                ntupleCuts += '%s && %s==1 && %s==1 && %s==1 && %s'\
                       % ( dataSampleCuts, sel, hlt1UB, hlt2B, trackChi2Cuts )
            elif selection == 'HLT1ExclBiased' :
                ntupleCuts += '%s && %s==1 && %s==1 && %s==1 && %s'\
                       % ( dataSampleCuts, sel, hlt1ExclB, hlt2B, trackChi2Cuts )
            elif selection == 'paper2012' :
                ntupleCuts += '%s && %s==1 && (%s==1 || %s==1) && %s==1 && %s'\
                       % ( dataSampleCuts, sel, hlt1B, hlt1UB, hlt2B, trackChi2Cuts )
            elif selection == 'timeEffFit' :
                ntupleCuts += '%s && %s==1 && (%s==1 || %s==1) && (%s==1 || %s==1) && %s'\
                       % ( dataSampleCuts, sel, hlt1B, hlt1UB, hlt2B, hlt2UB, trackChi2Cuts )
            else :
                raise ValueError( 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: unknown selection: "%s"' % selection )

            from P2VVGeneralUtils import readData
            self._dataSets['data'] = readData(  filePath = nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSetNTuple
                                              , Rename = 'JpsiphiData', ntupleCuts = ntupleCuts )

        else :
            self._dataSets['data'] = None

        if paramKKMass or makePlots :
            # create KK mass binning
            from array import array
            KKMassBinsArray = array( 'd', KKMassBinBounds )

            from ROOT import RooBinning
            self._KKMassBinning = RooBinning( len(KKMassBinBounds) - 1, KKMassBinsArray, 'KKMassBinning' )
            KKMass.setBinning( self._KKMassBinning, 'KKMassBinning' )

            # add KK mass split category to data
            from RooFitWrappers import BinningCategory
            self._KKMassCat = BinningCategory( 'KKMassCat'
                                              , Observable = KKMass
                                              , Binning = self._KKMassBinning
                                              , Fundamental = True
                                              , Data = [ self._dataSets['data'] ]
                                              , CatTypeName = 'bin'
                                             )


        ###################################################################################################################################
        ## initialize PDF component objects ##
        ######################################

        nEvents     = self._dataSets['data'].sumEntries() if self._dataSets['data'] else numEvents
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
        sigMassArgs = dict( Name = 'sig_m', mass = BMass )
        if sigMassModel.startswith('box') :
            from P2VVParameterizations.MassPDFs import Box_Signal_Mass as SignalBMass
            if sigMassModel.endswith('Fixed') :
                boxWidth = 0.5 * ( BMass.getMin('RightSideBand') - BMass.getMax('LeftSideBand') )
                boxMean  = BMass.getMax('LeftSideBand') + boxWidth
                sigMassArgs['m_sig_mean']  = dict( Value = boxMean,  Constant = True )
                sigMassArgs['m_sig_width'] = dict( Value = boxWidth, Constant = True, MinMax = ( 0.1, 2. * boxWidth ) )

        elif sigMassModel.startswith('Gauss') :
            from P2VVParameterizations.MassPDFs import Gauss_Signal_Mass as SignalBMass

        elif sigMassModel.startswith('DoubleGauss') :
            from P2VVParameterizations.MassPDFs import DoubleGauss_Signal_Mass as SignalBMass
            if sigMassModel.endswith('Diag') :
                sigMassArgs['TransformWidthPars'] = dict(  m_sig_frac    = ( +0.033129, -0.008339, -0.007473 )
                                                         , m_sig_sigma_1 = ( +0.115025, -0.067412, +0.000953 )
                                                         , m_sig_sigma_2 = ( +0.756560, +0.010614, +0.000182 )
                                                        )

        elif sigMassModel.startswith('CB') :
            from P2VVParameterizations.MassPDFs import CB_Signal_Mass as SignalBMass

        elif sigMassModel.startswith('DoubleCB') :
            from P2VVParameterizations.MassPDFs import DoubleCB_Signal_Mass as SignalBMass

        else :
            from P2VVParameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass

        bkgMassArgs = dict( Name = 'bkg_m', mass = BMass )
        if bkgMassModel.startswith('linear') :
            from P2VVParameterizations.MassPDFs import Linear_Background_Mass as BackgroundBMass
            if bkgMassModel.endswith('Constant') : bkgMassArgs['Constant'] = True
        else :
            from P2VVParameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass

        self._signalBMass     = SignalBMass(     **sigMassArgs )
        self._backgroundBMass = BackgroundBMass( **bkgMassArgs )

        from RooFitWrappers import buildPdf
        if SFit :
            self._sigMassComps += self._signalBMass.pdf()
            self._bkgMassComps += self._backgroundBMass.pdf()
            self._massPdf = buildPdf( [ self._sigMassComps, self._bkgMassComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )

        else :
            self._signalComps     += self._signalBMass.pdf()
            self._backgroundComps += self._backgroundBMass.pdf()
            self._massPdf = buildPdf( [ self._signalComps, self._backgroundComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )

        if self._dataSets['data'] :
            # determine mass parameters with a fit
            print 120 * '='
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: fitting with mass PDF'
            self._massFitResult = self._massPdf.fitTo( self._dataSets['data'], Save = True, **fitOpts )

            from P2VVImports import parNames
            self._massFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
            self._massFitResult.covarianceMatrix().Print()
            self._massFitResult.correlationMatrix().Print()

            if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit'] or paramKKMass == 'simultaneous' ) :
                # categories for splitting the PDF
                splitCats = [ [ ] ]
                splitCats[0] += [ hlt1ExclB ] if selection == 'paper2012' else [ hlt1ExclB, hlt2B ] if selection == 'timeEffFit' else [ ]
                splitCats[0] += [ self._KKMassCat ] if paramKKMass == 'simultaneous' else [ ]

                # get mass parameters that are split
                splitParams = [ [ par for par in self._massPdf.Parameters() if par.getAttribute('Yield') ] ]
                if 'FreeBkg' in SWeightsType and paramKKMass == 'simultaneous' :
                    splitCats.append( [ self._KKMassCat ] )
                    splitParams.append( [ par for par in self._backgroundBMass.parameters() if not par.isConstant() ] )

                # build simultaneous mass PDF
                from RooFitWrappers import SimultaneousPdf
                self._sWeightMassPdf = SimultaneousPdf(  self._massPdf.GetName() + '_simul'
                                                       , MasterPdf       = self._massPdf
                                                       , SplitCategories = splitCats
                                                       , SplitParameters = splitParams
                                                      )

                # set yields for categories
                splitCat     = self._sWeightMassPdf.indexCat()
                splitCatIter = splitCat.typeIterator()
                splitCatState = splitCatIter.Next()
                splitCatPars = self._sWeightMassPdf.getVariables()
                while splitCatState :
                    sigYield = getSplitPar( 'N_sigMass' if SFit else 'N_signal', splitCatState.GetName(), splitCatPars )
                    bkgYield = getSplitPar( 'N_bkgMass' if SFit else 'N_bkg',    splitCatState.GetName(), splitCatPars )

                    if splitCat.isFundamental() :
                        selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
                    else :
                        splitCat.setLabel( splitCatState.GetName() )
                        selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
                    nEv    = self._dataSets['data'].sumEntries()
                    nEvBin = self._dataSets['data'].sumEntries(selStr)

                    sigYield.setVal( sigYield.getVal() * nEvBin / nEv )
                    sigYield.setError( sqrt( sigYield.getVal() ) )
                    sigYield.setMin(0.)
                    sigYield.setMax(nEvBin)
                    bkgYield.setVal( bkgYield.getVal() * nEvBin / nEv )
                    bkgYield.setError( sqrt( bkgYield.getVal() ) )
                    bkgYield.setMin(0.)
                    bkgYield.setMax(nEvBin)

                    splitCatState = splitCatIter.Next()

                if SWeightsType.endswith( 'Fixed' ) :
                    # fix mass shape parameters for fit
                    fixedMassPars = [ par for par in self._sWeightMassPdf.Parameters()\
                                      if not ( par.getAttribute('Yield') or par.isConstant() or 'Category' in par.ClassName() ) ]
                    #from P2VVImports import parValues
                    #parVals = {  'm_bkg_exp_bin0'  : -2.1814e-03
                    #           , 'm_bkg_exp_bin1'  : -4.6151e-03
                    #           , 'm_bkg_exp_bin2'  : -1.4166e-03
                    #           , 'm_bkg_exp_bin3'  : -2.5203e-03
                    #           , 'm_bkg_exp_bin4'  : -1.3963e-03
                    #           , 'm_bkg_exp_bin5'  : -2.0610e-03
                    #           , 'm_sig_frac'      :  7.4479e-01
                    #           , 'm_sig_mean'      :  5.36809e+03
                    #           , 'm_sig_sigma_1'   :  6.1690e+00
                    #           , 'm_sig_sigma_sf'  :  2.0769e+00
                    #          }
                    for par in fixedMassPars :
                        #par.setVal( parValues[ par.GetName() ][0]\
                        #            + ( 3. * parValues[ par.GetName() ][1] if par.GetName().startswith('m_bkg_exp') else 0. ) )
                        #par.setError( parValues[ par.GetName() ][1] )
                        #par.setVal( parVals[ par.GetName() ] )
                        par.setConstant(True)

                # hack to do fit with fixed shape parameters first
                fixedShapeFit = False
                if fixedShapeFit :
                    from P2VVImports import parNames, parValues
                    fixedMassPars = [ par for par in self._sWeightMassPdf.Parameters()\
                                      if not ( par.getAttribute('Yield') or par.isConstant() or 'Category' in par.ClassName() ) ]
                    #parValues = {  'm_sig_mean'      : (  5368.236,    5.47e-02, -1.       )
                    #             , 'm_sig_frac'      : (  8.0283e-01,  2.73e-02, -1.       )
                    #             , 'm_sig_sigma_1'   : (  6.2728,      1.19e-01, -1.       )
                    #             , 'm_sig_sigma_sf'  : (  2.2479,      1.24e-01, -1.       )
                    #             , 'm_bkg_exp'       : ( -1.6249e-03,  9.67e-05, -1.       )
                    #            }

                    #fixedMassParVals = { }
                    for par in fixedMassPars :
                        par.setConstant(True)
                        #fixedMassParVals[par.GetName()] = ( par.getVal(), par.getError() )
                        par.setVal( parValues[par.GetName()][0] )
                        par.setError( parValues[par.GetName()][1] )

                    massNLL = self._sWeightMassPdf.createNLL( self._dataSets['data'] )
                    self._simMassFitResult = self._sWeightMassPdf.fitTo( self._dataSets['data'], Save = True, **fitOpts )
                    self._simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
                    massNLLValNom = massNLL.getVal()
                    for par in fixedMassPars :
                        par.setConstant(False)
                        #par.setVal( fixedMassParVals[par.GetName()][0] )
                        #par.setError( fixedMassParVals[par.GetName()][1] )

                # determine mass parameters in each sub-sample with a fit
                print 120 * '='
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: fitting with simultaneous mass PDF'
                self._simMassFitResult = self._sWeightMassPdf.fitTo( self._dataSets['data'], Save = True, **fitOpts )

                from P2VVImports import parValues
                self._simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
                self._simMassFitResult.covarianceMatrix().Print()
                self._simMassFitResult.correlationMatrix().Print()

                # hack to do fit with fixed shape parameters first
                if fixedShapeFit :
                    massNLLValThis = massNLL.getVal()
                    print '  mass NLL values:'
                    print '  nominal: %f;  this: %f; 2*DeltaNLL = %f'\
                          % ( massNLLValNom, massNLLValThis, 2. * ( massNLLValNom - massNLLValThis ) )

                if SWeightsType.endswith( 'Fixed' ) :
                    # free parameters that were fixed for mass fit
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: constant parameters in mass fit:'
                    for par in fixedMassPars :
                        par.Print()
                        par.setConstant(False)

            else :
                splitCatPars = None
                self._sWeightMassPdf = self._massPdf


        ###################################################################################################################################
        ## compute S-weights and create signal and background data sets ##
        ##################################################################

        if self._dataSets['data'] :
            print 120 * '='
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: computing S-weights'

            # create sWeigthed data sets
            from P2VVGeneralUtils import SData, splot
            self._SData = SData( Pdf = self._sWeightMassPdf, Data = self._dataSets['data'], Name = 'massSData' )
            self._dataSets['SWeightData']    = self._SData.data()
            self._dataSets['sigSWeightData'] = self._SData.data( 'sigMass' if SFit else 'signal' )
            self._dataSets['bkgSWeightData'] = self._SData.data( 'bkgMass' if SFit else 'bkg'    )

            # print signal/background info to screen
            splitCats = [  self._dataSets['data'].get().find( hlt1ExclB.GetName() )
                         , self._dataSets['data'].get().find( hlt2B.GetName() )
                        ]
            if hasattr( self, '_KKMassCat' ) and self._KKMassCat.numTypes() > 1 :
                splitCats.append( self._dataSets['data'].get().find( self._KKMassCat.GetName() ) )
            splitCats = [ cat for cat in splitCats if cat ]
            self._dataSets['sigSWeightData'].Print()
            self._dataSets['bkgSWeightData'].Print()
            nEv    = self._dataSets['data'].sumEntries()
            nSigEv = self._dataSets['sigSWeightData'].sumEntries()
            nBkgEv = self._dataSets['bkgSWeightData'].sumEntries()
            S_B    = nSigEv / nBkgEv
            S_SB   = nSigEv / nEv
            signif = nSigEv / sqrt(nEv)
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: number of events:'
            print '                  | total | signal   | backgr.  | S/B    | S/(S+B) | S/sqrt(S+B)'
            print '    --------------|-------|----------|----------|--------|---------|------------'
            print '            total | %5.0f | %8.2f | %8.2f | %6.4f | %6.4f  | %7.3f' % ( nEv, nSigEv, nBkgEv, S_B, S_SB, signif )
            print '    --------------|-------|----------|----------|--------|---------|------------'

            if splitCats :
                iters = { }
                inds  = { }
                labs  = { }
                for cat in splitCats :
                    iters[cat] = 0
                    inds[cat]  = [ ]
                    labs[cat]  = [ ]
                    catIter = cat.typeIterator()
                    catState = catIter.Next()
                    while catState :
                        inds[cat].append( catState.getVal() )
                        labs[cat].append( catState.GetName() )

                        cut    = '!(%s-%d)' % ( cat.GetName(), catState.getVal() )
                        nEv    = self._dataSets['data'].sumEntries(cut)
                        nSigEv = self._dataSets['sigSWeightData'].sumEntries(cut)
                        nBkgEv = self._dataSets['bkgSWeightData'].sumEntries(cut)
                        S_B    = nSigEv / nBkgEv
                        S_SB   = nSigEv / nEv
                        signif = nSigEv / sqrt(nEv)
                        print '    %13s | %5.0f | %8.2f | %8.2f | %6.4f | %6.4f  | %7.3f'\
                              % ( catState.GetName(), nEv, nSigEv, nBkgEv, S_B, S_SB, signif )

                        catState = catIter.Next()

                    print '    --------------|-------|----------|----------|--------|---------|------------'

                if len(splitCats) > 1 :
                    cont = True
                    while cont :
                        cut    = '&&'.join( '!(%s-%d)' % ( cat.GetName(), inds[cat][ iters[cat] ] ) for cat in splitCats )
                        nEv    = self._dataSets['data'].sumEntries(cut)
                        nSigEv = self._dataSets['sigSWeightData'].sumEntries(cut)
                        nBkgEv = self._dataSets['bkgSWeightData'].sumEntries(cut)
                        S_B    = nSigEv / nBkgEv
                        S_SB   = nSigEv / nEv
                        signif = nSigEv / sqrt(nEv)
                        print '    %13s | %5.0f | %8.2f | %8.2f | %6.4f | %6.4f  | %7.3f'\
                               % ( ';'.join( labs[cat][ iters[cat] ] for cat in splitCats ), nEv, nSigEv, nBkgEv, S_B, S_SB, signif )

                        iters[ splitCats[-1] ] += 1
                        for catIt in range( len(splitCats) ) :
                            if iters[ splitCats[ -catIt - 1 ] ] >= splitCats[ -catIt - 1 ].numTypes() :
                                if catIt == len(splitCats) - 1 :
                                    cont = False
                                else :
                                    iters[ splitCats[ -catIt - 1 ] ] = 0
                                    iters[ splitCats[ -catIt - 2 ] ] +=1
                            else :
                                continue

                    print '    --------------|-------|----------|----------|--------|---------|--------'
                print

                if splitCatPars :
                    # print yields and fractions with error from S/(S+B) fraction only (no Poisson error for number of events!)
                    print '                  | total |        signal       |     background      |        S/B        |       S/(S+B)     |   S/sqrt(S+B)   '
                    print '    --------------|-------|---------------------|---------------------|-------------------|-------------------|-----------------'
                    for cat in iters.iterkeys() : iters[cat] = 0
                    cont = True
                    while cont :
                        stateName = ';'.join( labs[cat][ iters[cat] ] for cat in splitCats )
                        sigYield = getSplitPar( 'N_sigMass' if SFit else 'N_signal', '{%s}' % stateName, splitCatPars )
                        bkgYield = getSplitPar( 'N_bkgMass' if SFit else 'N_bkg',    '{%s}' % stateName, splitCatPars )

                        nSigEv = sigYield.getVal()
                        nBkgEv = bkgYield.getVal()
                        nEv    = nSigEv + nBkgEv
                        S_SB   = nSigEv / nEv
                        S_B    = nSigEv / nBkgEv
                        signif = nSigEv / sqrt(nEv)

                        nSigErr     = sigYield.getError()
                        nSigErrCorr = sqrt( nSigErr**2 - nSigEv**2 / nEv )
                        S_SBErr     = nSigErrCorr / nEv
                        S_BErr      = S_SBErr / ( 1 - S_SB )**2
                        signifErr   = S_SBErr * sqrt(nEv)
                        print '    %13s | %5.0f | %8.2f +/- %6.2f | %8.2f +/- %6.2f | %6.4f +/- %6.4f | %6.4f +/- %6.4f | %6.2f +/- %4.2f'\
                               % ( ';'.join( labs[cat][ iters[cat] ] for cat in splitCats ), nEv, nSigEv, nSigErrCorr\
                                  , nBkgEv, nSigErrCorr, S_B, S_BErr, S_SB, S_SBErr, signif, signifErr )

                        iters[ splitCats[-1] ] += 1
                        for catIt in range( len(splitCats) ) :
                            if iters[ splitCats[ -catIt - 1 ] ] >= splitCats[ -catIt - 1 ].numTypes() :
                                if catIt == len(splitCats) - 1 :
                                    cont = False
                                else :
                                    iters[ splitCats[ -catIt - 1 ] ] = 0
                                    iters[ splitCats[ -catIt - 2 ] ] +=1
                            else :
                                continue
                    print '    --------------|-------|---------------------|---------------------|-------------------|-------------------|-----------------'

            # create signal and background data sets with side band ranges
            self._dataSets['sigRangeData'] = self._dataSets['data'].reduce( CutRange = 'Signal'       )
            self._dataSets['bkgRangeData'] = self._dataSets['data'].reduce( CutRange = 'LeftSideBand' )
            self._dataSets['bkgRangeData'].append( self._dataSets['data'].reduce( CutRange = 'RightSideBand' ) )

            if makePlots :
                if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit']\
                      or paramKKMass == 'simultaneous' ) :
                    # create projection data set
                    indexCat = self._sWeightMassPdf.indexCat()
                    if indexCat.isFundamental() :
                        projWDataSet = [ indexCat ]
                    else :
                        projWDataSet = [ cat for cat in indexCat.getObservables( self._dataSets['data'] ) ]

                    projWData = dict( ProjWData = ( self._dataSets['data'].reduce( ArgSet = projWDataSet ), False ) )
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: projection data set for mumuKK mass plots:'
                    projWData['ProjWData'][0].Print()

                else :
                    projWData = dict()

                # plot mumuKK mass distributions
                self._massCanvs = [  TCanvas( 'massCanvLog',   'B mass logarithmic scale' )
                                   , TCanvas( 'massCanvSig',   'B mass signal range'      )
                                   , TCanvas( 'massCanvLeft',  'B mass left side band'    )
                                   , TCanvas( 'massCanvRight', 'B mass right side band'   )
                                  ]
                for ( pad, frameRange, nBins, plotTitle, logy, scale )\
                      in zip(  self._massCanvs
                             , [ '', 'Signal', 'LeftSideBand', 'RightSideBand' ]
                             , [ sum(numBMassBins) ] + numBMassBins
                             , [  BMass.GetTitle()
                                , BMass.GetTitle() + ' mass fit - signal'
                                , BMass.GetTitle() + ' mass fit - left side band'
                                , BMass.GetTitle() + ' mass fit - right side band'
                               ]
                             , [ True, False, False, False ]
                             , [ ( 8.e1, 1.e4 ), ( None, None ), ( None, None ), ( None, None ) ] # [ ( 8.e1, 1.e4 ), ( 0., 4500. ), ( 0., 660. ), ( 0., 640. ) ]
                            ) :
                    plot(  pad, BMass, self._dataSets['data'], self._sWeightMassPdf, logy = logy, yScale = scale, plotResidHist = True
                         , normalize = True, symmetrize = True
                         , frameOpts  = dict( Range = frameRange, Bins = nBins, Title = plotTitle )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.6, LineWidth = 2 )
                         , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 3 )
                         , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed, LineWidth = 3 )
                                         , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed, LineWidth = 3 )
                                        }
                        )

                if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit']\
                      or paramKKMass == 'simultaneous' ) :
                    # get simultaneous PDFs
                    indexCatIter  = indexCat.typeIterator()
                    indexCatState = indexCatIter.Next()
                    bins = [ ]
                    pdfs = [ ]
                    while indexCatState :
                        indexCat.setIndex( indexCatState.getVal() )
                        bins.append( [ ( indexCat.GetName(), indexCatState.getVal(), indexCatState.GetName() ) ] )
                        if indexCat.isFundamental() :
                            bins[-1].append( bins[-1][0] )
                        else :
                            for cat in indexCat.getObservables( self._dataSets['data'] ) :
                                bins[-1].append( ( cat.GetName(), cat.getIndex(), cat.getLabel() ) )

                        pdfs.append( self._sWeightMassPdf.getPdf( indexCatState.GetName() ) )
                        indexCatState = indexCatIter.Next()

                    # plot mumuKK mass distributions in KK mass bins
                    if   len(bins) <= 4 : nPads = ( 2, 2 )
                    elif len(bins) <= 6 : nPads = ( 3, 2 )
                    elif len(bins) <= 9 : nPads = ( 3, 3 )
                    else :                nPads = ( 4, 3 )
                    self._massCanvs.append( TCanvas( 'massCanvSigBins', 'B mass signal range bins' ) )
                    for ( pad, pdf, plotTitle, dataCuts, norm )\
                            in zip(  self._massCanvs[-1].pads( nPads[0], nPads[1] )
                                   , pdfs
                                   , [ BMass.GetTitle() + ' bin %d - signal' % bin[0][1] for bin in bins ]
                                   , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                                   , [ self._dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                                       / self._dataSets['data'].sumEntries() for bin in bins ]
                                  ) :
                        plot(  pad, BMass, self._dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                             , frameOpts  = dict( Range = 'Signal', Bins = numBMassBins[0], Title = plotTitle )
                             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                             , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                             , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                                             , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                                            }
                            )

                    self._massCanvs.append( TCanvas( 'massCanvLeftBins', 'B mass left side band bins' ) )
                    for ( pad, pdf, plotTitle, dataCuts, norm )\
                            in zip(  self._massCanvs[-1].pads( nPads[0], nPads[1] )
                                   , pdfs
                                   , [ BMass.GetTitle() + ' bin %d - left side band' % bin[0][1] for bin in bins ]
                                   , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                                   , [ self._dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                                       / self._dataSets['data'].sumEntries() for bin in bins ]
                                  ) :
                        plot(  pad, BMass, self._dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                             , frameOpts  = dict( Range = 'LeftSideBand', Bins = numBMassBins[1], Title = plotTitle )
                             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                             , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                             , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                                             , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                                            }
                            )

                    self._massCanvs.append( TCanvas( 'massCanvRightBins', 'B mass right side band bins' ) )
                    for ( pad, pdf, plotTitle, dataCuts, norm )\
                            in zip(  self._massCanvs[-1].pads( nPads[0], nPads[1] )
                                   , pdfs
                                   , [ BMass.GetTitle() + ' bin %d - right side band' % bin[0][1] for bin in bins ]
                                   , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                                   , [ self._dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                                       / self._dataSets['data'].sumEntries() for bin in bins ]
                                  ) :
                        plot(  pad, BMass, self._dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                             , frameOpts  = dict( Range = 'RightSideBand', Bins = numBMassBins[2], Title = plotTitle )
                             , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                             , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                             , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = kDashed )
                                             , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = kDashed )
                                            }
                            )

        else :
            self._dataSets['sigSWeightData'] = None
            self._dataSets['bkgSWeightData'] = None
            self._dataSets['sigRangeData']   = None
            self._dataSets['bkgRangeData']   = None


        ###################################################################################################################################
        ## build tagging categories ##
        ##############################

        if nominalPdf or not self._iTagZeroTrick :
            print 120 * '='
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building tagging categories'

            # build tagging categories opposite side
            from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
            if not nominalPdf and contEstWTag :
                self._tagCatsOS = TaggingCategories(  tagCat = 'tagCatP2VVOS', DataSet = self._dataSets['sigSWeightData']
                                                    , estWTag = estWTagOS, wTagP0Constraint = 'constrain' if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = 'constrain' if nominalPdf else constrainTagging )
                self._tagCatsSS = TaggingCategories(  tagCat = 'tagCatP2VVSS', DataSet = self._dataSets['sigSWeightData']
                                                    , estWTag = estWTagSS, SameSide = True
                                                    , wTagP0Constraint = 'constrain' if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = 'constrain' if nominalPdf else constrainTagging )
            else :
                self._tagCatsOS = TaggingCategories(  tagCat = 'tagCatP2VVOS', DataSet = self._dataSets['sigSWeightData']
                                                    , estWTagName = estWTagOS.GetName()
                                                    , TagCats = tagCatsOS, NumSigmaTagBins = 1.
                                                    , wTagP0Constraint = 'constrain' if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = 'constrain' if nominalPdf else constrainTagging )
                self._tagCatsSS = TaggingCategories(  tagCat = 'tagCatP2VVSS', DataSet = self._dataSets['sigSWeightData']
                                                    , SameSide = True
                                                    , estWTagName = estWTagSS.GetName()
                                                    , TagCats = tagCatsSS, NumSigmaTagBins = 1.
                                                    , wTagP0Constraint = 'constrain' if nominalPdf else constrainTagging
                                                    , wTagP1Constraint = 'constrain' if nominalPdf else constrainTagging )

            tagCatP2VVOS = self._tagCatsOS['tagCat']
            tagCatP2VVOS.setIndex(1)
            observables[tagCatP2VVOS.GetName()] = tagCatP2VVOS

            tagCatP2VVSS = self._tagCatsSS['tagCat']
            tagCatP2VVSS.setIndex(1)
            observables[tagCatP2VVSS.GetName()] = tagCatP2VVSS

            if not self._iTagZeroTrick and ( nominalPdf or multiplyByTagPdf or not condTagging ) :
                obsSetP2VV.append(tagCatP2VVOS)
                if SSTagging : obsSetP2VV.append(tagCatP2VVSS)

            if nominalPdf or condTagging :
                self._tagCatsOS.addConditional(tagCatP2VVOS)
                self._tagCatsOS.addConditional(iTagOS)
                self._tagCatsSS.addConditional(tagCatP2VVSS)
                self._tagCatsSS.addConditional(iTagSS)

            # add tagging categories to data sets
            from P2VVGeneralUtils import addTaggingObservables
            for dataKey, data in self._dataSets.iteritems() :
                if data and not data.get().find( iTagOS.GetName() ) :
                    addTaggingObservables( data, iTagOS.GetName(), tagCatP2VVOS.GetName(), tagDecisionOS.GetName(), estWTagOS.GetName()
                                          , self._tagCatsOS['tagCats'] )
                    addTaggingObservables( data, iTagSS.GetName(), tagCatP2VVSS.GetName(), tagDecisionSS.GetName(), estWTagSS.GetName()
                                          , self._tagCatsSS['tagCats'] )

            from P2VVParameterizations.FlavourTagging import Combined_TaggingCategories as CombTaggingCategories
            self._tagCatsComb = CombTaggingCategories( self._tagCatsOS, self._tagCatsSS )

            if self._dataSets['sigSWeightData'] :
                # print tagging categories distribution for signal and background
                from RooFitWrappers import ArgSet
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in opposite side tagging category for signal:'
                self._dataSets['sigSWeightData'].table( ArgSet( 'sigOSTagSet', [ tagCatP2VVOS, iTagOS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in opposite side tagging category for background:'
                self._dataSets['bkgSWeightData'].table( ArgSet( 'bkgOSTagSet', [ tagCatP2VVOS, iTagOS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in same side tagging category for signal:'
                self._dataSets['sigSWeightData'].table( ArgSet( 'sigSSTagSet', [ tagCatP2VVSS, iTagSS ] ) ).Print('v')
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: distribution in same side tagging category for background:'
                self._dataSets['bkgSWeightData'].table( ArgSet( 'bkgSSTagSet', [ tagCatP2VVSS, iTagSS ] ) ).Print('v')


        ###################################################################################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################

        print 120 * '='
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building PDFs'

        # transversity amplitudes
        commonArgs = dict( AmbiguityParameters = ambiguityPars )
        if paramKKMass == 'functions' :
            commonArgs[ 'KKMass' ]        = KKMass
            commonArgs[ 'KKMassBinning' ] = self._KKMassBinning
        if nominalPdf or not ASParam.startswith('Mag2ReIm') :
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
        self._lifetimeParams = LifetimeParams( dGamma = dGammaVar, dMConstraint = 'constrain' if nominalPdf else constrainDeltaM )
        if ambiguityPars : self._lifetimeParams['dGamma'].setVal( -self._lifetimeParams['dGamma'].getVal() )

        if nominalPdf or timeResType.startswith('event') :
            from P2VVParameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
            timeResArgs = dict(  time = time
                               , timeResSigma = timeRes
                               , timeResSFConstraint = 'constrain' if nominalPdf else constrTResScale
                               , Cache = multiplyByTimeEff not in [ 'all', 'signal', 'background' ]  # make sure we do not 'double cache'
                              )
            if not nominalPdf and 'nomean' in timeResType.lower() :
                timeResArgs['timeResMean']   = ConstVar( Name = 'timeResMean',   Value = 0. )
                timeResArgs['timeResMeanSF'] = ConstVar( Name = 'timeResMeanSF', Value = 1. )
            elif not nominalPdf and 'constmean' in timeResType.lower() :
                timeResArgs['timeResMean']   = dict( Value = -0.01, Error = 0.005 )
                timeResArgs['timeResMeanSF'] = ConstVar( Name = 'timeResMeanSF', Value = 1. )
                timeResArgs['timeResMeanConstraint'] = constrTResScale
            elif not nominalPdf and 'stoffset' in timeResType.lower():
                timeResArgs['timeResMeanConstraint'] = 'constrain'
                timeResArgs['timeResSigmaSF'] = dict(Name = 'timeResSigmaSF', Value = 1.30, Error = 0.06, MinMax = ( 0.1, 5. ) )
                ## timeResArgs['timeResSigmaSF'] = ConstVar( Name = 'timeResMeanSF', Value = 1.30 )
                timeResArgs['timeResSigmaOffset'] = dict(Name = 'timeResSigmaOffset', Value = 0.0065, Error = 0.001, MinMax = ( 0.00001, 0.1 ) )
                ## timeResArgs['timeResSigmaOffset'] = ConstVar(Name = 'timeResSigmaOffset', Value = 0.012)
                timeResArgs['timeResSFOffset'] = True
            else :
                timeResArgs['timeResMeanConstraint'] = 'constrain' if nominalPdf else constrTResScale

            self._timeResModel = TimeResolution( **timeResArgs )

        elif timeResType == '3Gauss' :
            from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution( time = time
                                                , timeResSFConstraint = 'constrain' if nominalPdf else constrTResScale )

        else :
            from P2VVParameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution(  time          = time
                                                , timeResMu     = dict( Value = 0.,    Constant = True )
                                                , timeResSigma  = dict( Value = 0.045, Constant = True )
                                                , PerEventError = False
                                                , Cache = multiplyByTimeEff not in [ 'all', 'signal', 'background' ]
                                               )

        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: decay time resolution model:'
        self._timeResModel['model'].Print()
        for par in self._timeResModel.parameters() : par.Print()

        if makePlots :
            # plot event-by-event estimated time resolution
            self._timeResCanv = TCanvas( 'timeResCanv', 'Decay time resolution' )
            for ( pad, data, plotTitle )\
                  in zip(  self._timeResCanv.pads( 2, 2 )
                         , [ self._dataSets['sigSWeightData'], self._dataSets['bkgSWeightData'] ]
                         , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                        ) :
                plot(  pad, timeRes, data, None
                     , frameOpts  = dict( Title = timeRes.GetTitle() + plotTitle )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 , MarkerColor = kBlue, LineColor = kBlue )
                    )


        # CP violation parameters
        if lambdaCPParam == 'ReIm' : 
            from P2VVParameterizations.CPVParams import LambdaCarth_CPParam as CPParam
            ReLambdaCPVar = dict( Name = 'ReLambdaCP' )
            ImLambdaCPVar = dict( Name = 'ImLambdaCP' )
            if blind: ReLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsGoofyMoriond2012', 0.1 )
            if blind: ImLambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsPlutoMoriond2012', 0.1 )
            self._lambdaCP = CPParam( ReLambdaCP = ReLambdaCPVar, ImLambdaCP = ImLambdaCPVar )

        elif lambdaCPParam == 'lambPhiRel_CPVDecay' :
            from P2VVParameterizations.CPVParams import LambdaAbsArgRel_CPVDecay_CPParam as CPParam
            phiCPVar = dict( Name = 'phiCP_m' )
            if blind: phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
            self._lambdaCP = CPParam( phiCP_m = phiCPVar, AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes )
            if ambiguityPars : self._lambdaCP['phiCP_m'].setVal( pi - self._lambdaCP['phiCP_m'].getVal() )

        elif lambdaCPParam == 'lambPhi_CPVDecay' :
            from P2VVParameterizations.CPVParams import LambdaAbsArg_CPVDecay_CPParam as CPParam
            self._lambdaCP = CPParam( AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes )

            #rhoCP = RealVar( 'rhoCP', Title = 'CPV in decay param. |rho|', Value = 1., Error = 0.04, MinMax = ( 0.,      5.      ) )
            #phiCP = RealVar( 'phiCP', Title = 'CPV in decay param. phi',   Value = 0., Error = 0.1,  MinMax = ( -RooInf, +RooInf ) )
            #self._lambdaCP = CPParam(  AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes
            #                         , rhoCP_A0 = rhoCP, rhoCP_Apar = rhoCP, rhoCP_Aperp = rhoCP, rhoCP_AS = rhoCP
            #                         , phiCP_A0 = phiCP, phiCP_Apar = phiCP, phiCP_Aperp = phiCP, phiCP_AS = phiCP
            #                        )

        else :
            if lambdaCPParam == 'lambPhi' :
                from P2VVParameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
            else :
                from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam

            phiCPVar = dict( Name = 'phiCP' )
            lambdaCPVar = dict( Name = 'lambdaCP' )
            if blind:
                phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
                lambdaCPVar['Blind'] = ( 'UnblindUniform', 'BsCrumbleMoriond2012', 0.3 )
            self._lambdaCP = CPParam( phiCP = phiCPVar, lambdaCP = lambdaCPVar )
            if ambiguityPars :
                self._lambdaCP['phiCP'].setVal( pi - self._lambdaCP['phiCP'].getVal() )

        # coefficients for time functions
        from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, self._amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

        # tagging parameters
        if not nominalPdf and self._iTagZeroTrick :
            from P2VVParameterizations.FlavourTagging import LinearEstWTag_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams(  estWTag = estWTagOS, p0 = dict( Name = 'wTagP0' ), p1 = dict( Name = 'wTagP1' )
                                                , p0Constraint = 'constrain' if nominalPdf else constrainTagging
                                                , p1Constraint = 'constrain' if nominalPdf else constrainTagging )

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
            if SSTagging : tagCatsDictSS = self._tagCatsSS.tagCatsDict()

            if not nominalPdf and sigTaggingPdf.startswith('tagUntag') :
                # assume products of asymmetries are small and B-Bbar asymmetries are equal for all tagged categories
                if tagCatP2VVOS.numTypes() > 2 or ( SSTagging and tagCatP2VVSS.numTypes() > 2 ) :
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: tagging in signal PDF:\n'\
                        + '    * assuming B-Bbar asymmetries are equal for all tagged categories'

                # provide the same asymmetry for all tagged categories
                asymVal = 0. #asymVal = -self._lambdaCP['C'].getVal()
                asymErr = ( 10. / sqrt( self._dataSets['sigSWeightData'].sumEntries() ) ) if self._dataSets['sigSWeightData'] else 0.1
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
                        tagCatsDict[ 'tagCatCoef0_%d'  % catOS ] = tagCatsDictOS[ 'tagCatCoef%d'  % catOS ]
                        tagCatsDict[ 'ATagEff0_%d'     % catOS ] = tagCatsDictOS[ 'ATagEff%d'     % catOS ]
                        tagCatsDict[ 'tagDilution0_%s' % catOS ] = tagCatsDictOS[ 'tagDilution%d' % catOS ]
                        tagCatsDict[ 'ADilWTag0_%s'    % catOS ] = tagCatsDictOS[ 'ADilWTag%d'    % catOS ]

                    for catSS in range( 1, tagCatsDictSS['NumTagCats'] ) :
                        tagCatsDict[ 'tagCatCoef1_%d'  % catSS ] = tagCatsDictSS[ 'tagCatCoef%d'  % catSS ]
                        tagCatsDict[ 'ATagEff1_%d'     % catSS ] = tagCatsDictSS[ 'ATagEff%d'     % catSS ]
                        tagCatsDict[ 'tagDilution1_%s' % catSS ] = tagCatsDictSS[ 'tagDilution%d' % catSS ]
                        tagCatsDict[ 'ADilWTag1_%s'    % catSS ] = tagCatsDictSS[ 'ADilWTag%d'    % catSS ]

                # add production asymmetry and normalization asymmetry to tagging categories dictionary
                tagCatsDict['AProd'] = 0.
                tagCatsDict['ANorm'] = 0. #-self._lambdaCP['C'].getVal()

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
        self._BTagDecay = sigPdf


        ###################################################################################################################################
        ## time acceptance ##
        #####################

        if multiplyByTimeEff in [ 'all', 'signal', 'background' ] :
            self._timeResModelOriginal = self._timeResModel
            if timeEffType == 'fit' and selection == 'timeEffFit' :
                hists = {  hlt1ExclB : {  'exclB'    : { 'histogram' : 'hlt1_shape', 'average' : ( 6.285e-01, 1.633e-02 ) }
                                        , 'notExclB' : { 'bins'      : time.getRange(), 'heights' : [0.5]                 }
                                       }
                         , hlt2B     : { 'B'         : { 'histogram' : 'hlt2_shape', 'average' : ( 6.3290e-01, 1.65e-02 ) } }
                         , hlt2UB    : { 'UB'        : { 'bins'      : time.getRange(), 'heights' : [0.5]                 } }
                        }

                from P2VVParameterizations.TimeAcceptance import Paper2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance( time = time, Input = timeEffHistFile, Histograms = hists
                                                    , Data = self._dataSets['data'], Fit = True, Original = sigPdf
                                                    , ResolutionModel = self._timeResModel )

            elif timeEffType == 'paper2012' and selection == 'paper2012' :
                hists = { hlt1ExclB : {  'exclB'    : { 'histogram' : timeEffHistExclBName }
                                       , 'notExclB' : { 'histogram' : timeEffHistUBName    }
                                      }
                        }

                from P2VVParameterizations.TimeAcceptance import Paper2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance( time = time, Input = timeEffHistFile, Histograms = hists
                                                    , Data = self._dataSets['data'], Fit = False, Original = sigPdf
                                                    , ResolutionModel = self._timeResModel, BinHeightMinMax = ( -RooInf, RooInf ) )

            elif timeEffType in [ 'HLT1Unbiased', 'HLT1ExclBiased' ] or ( timeEffType == 'paper2012' and selection == 'paper2012' ) :
                from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance(  time = time
                                                    , Input = timeEffHistFile
                                                    , Histogram = timeEffHistExclBName if timeEffType == 'HLT1ExclBiased'\
                                                                  else timeEffHistUBName
                                                    , ResolutionModel = self._timeResModel
                                                   )

            else:
                raise ValueError( 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: unknown time efficiency type: "%s" (with "%s" selection)'\
                                 % ( timeEffType, selection ) )

        if multiplyByTimeEff in [ 'all', 'signal' ] :
            # multiply signal PDF with time acceptance
            print 'P2VV - INFO:  Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with lifetime efficiency function'
            args['resolutionModel'] = self._timeResModel['model']
            args['ConditionalObservables'] = list( set( args['ConditionalObservables'] + self._timeResModel.conditionalObservables() ) )
            args['ExternalConstraints'] = list( set( args['ExternalConstraints'] + self._timeResModel.externalConstraints() ) )
            sigPdfTimeAcc = BTagDecay( 'sig_t_angles_timeEff', **args )
        else :
            sigPdfTimeAcc = sigPdf


        ###################################################################################################################################
        ## angular acceptance ##
        ########################

        if multiplyByAngEff :
            # multiply signal PDF with angular efficiency
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with angular efficiency moments from file "%s"'\
                  % angEffMomentsFile

            from P2VVGeneralUtils import RealMomentsBuilder
            moments = RealMomentsBuilder()
            if multiplyByAngEff == 'weights' :
                angMomInds = [ ( 0, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 0, 0 ), ( 0, 2, 1 ), ( 0, 2, -1 ), ( 0, 2, -2 )
                              , ( 1, 0, 0 ), ( 1, 2, 1 ), ( 1, 2, -1 ) ]
            elif multiplyByAngEff == 'basis012' :
                angMomInds = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                              for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            elif multiplyByAngEff == 'basis0123' :
                angMomInds = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(4) for YIndex0 in range(4)\
                              for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
            elif multiplyByAngEff == 'basisSig3' :
                angMomInds = [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if nominalPdf or not transAngles\
                              else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
            elif multiplyByAngEff == 'basisSig6' :
                angMomInds = [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 2, 2, 0 ), ( 1, 1, 0 ), ( 1, 2, 0 ) ]\
                             if nominalPdf or not transAngles else\
                             [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 2, 0 ), ( 2, 2, 2 )\
                              , ( 1, 1, 1 ), ( 1, 2, 0 ), ( 1, 2, 2 ) ]
            else :
                raise RuntimeError('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: not a valid angular efficiency configuration: %s'\
                                   % multiplyByAngEff)

            moments.appendPYList( self._angleFuncs.angles, angMomInds )
            moments.read(angEffMomentsFile)
            moments.Print()
            sigPdfTimeAcc = moments * sigPdfTimeAcc

        self._BTagDecayAngEff = sigPdfTimeAcc

        # add PDF to signal components
        self._signalComps += sigPdfTimeAcc
        if nominalPdf or condTagging : self._sig_t_angles = sigPdfTimeAcc
        else :                         self._sig_t_angles_tagCat_iTag = sigPdfTimeAcc


        ###################################################################################################################################
        ## plot mumu mass ##
        ####################

        if makePlots :
            self._mumuMassCanv = TCanvas( 'mumuMassCanv', 'mumu Mass' )
            for ( pad, data, plotTitle )\
                  in zip(  self._mumuMassCanv.pads( 2, 2 )
                         , [ self._dataSets['sigSWeightData'], self._dataSets['bkgSWeightData'] ]
                         , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                        ) :
                plot(  pad, mumuMass, data, None
                     , frameOpts  = dict( Title = mumuMass.GetTitle() + plotTitle )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 , MarkerColor = kBlue, LineColor = kBlue )
                    )


        ###################################################################################################################################
        ## build KK mass PDFs ##
        ########################

        if paramKKMass or makePlots :
            # build the signal and background KK mass PDFs
            from P2VVParameterizations.MassPDFs import Binned_MassPdf
            self._signalKKMass = Binned_MassPdf( 'sig_mKK', KKMass, Binning = self._KKMassBinning
                                                , Data = self._dataSets['sigSWeightData'] )
            self._backgroundKKMass = Binned_MassPdf( 'bkg_mKK', KKMass, Binning = self._KKMassBinning
                                                    , Data = self._dataSets['bkgSWeightData'] )
            #if paramKKMass == 'functions' :
            #    self._signalComps += self._signalKKMass.pdf()
            #    if not SFit: self._backgroundComps += self._backgroundKKMass.pdf()

            if makePlots :
                self._KKMassCanv = TCanvas( 'KKMassCanv', 'KK Mass' )
                for ( pad, data, pdf, plotTitle, scale )\
                      in zip(  self._KKMassCanv.pads( 2, 2 )
                             , [ self._dataSets['sigSWeightData'], self._dataSets['bkgSWeightData'] ]
                             , [ self._signalKKMass.pdf(), self._backgroundKKMass.pdf() ]
                             , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                             , [ ( 5., 1.e4 ), ( 1.2e2, 1.2e3 ) ]
                            ) :
                    plot(  pad, KKMass, data, pdf, logy = True, yScale = scale
                         , frameOpts  = dict( Title = KKMass.GetTitle() + plotTitle )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )#, MarkerColor = kBlue, LineColor = kBlue )
                        )


        ###################################################################################################################################
        ## build signal tagging PDF ##
        ##############################

        if ( nominalPdf and not SFit ) or ( condTagging and multiplyByTagPdf ) :
            if not nominalPdf and self._iTagZeroTrick :
                # no implementation for signal tagging PDF with tag = { B, Bbar, Untagged } (yet)
                pass

            else :
                # tagCat = { Untagged, TagCat1, TagCat2, ... }, tag = { B, Bbar }
                sigTaggingData = self._dataSets['sigRangeData'] if massRangeBackground else self._dataSets['sigSWeightData']
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
                        if tagCatP2VVOS.numTypes() > 2 or ( SSTagging and tagCatP2VVSS.numTypes() > 2 ) :
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
                    self._sigTaggingPdf = TaggingPdf(  'tagCat_iTag'
                                                     , tagCatP2VVOS, tagCatP2VVSS if SSTagging else None
                                                     , iTagOS, iTagSS if SSTagging else None
                                                     , NamePrefix    = 'sig'
                                                     , TagCatCoefs   = [ sigPdf.tagCatCoefs(i) for i in range( tagCatP2VVOS.numTypes() ) ]\
                                                                       if SSTagging else [ sigPdf.tagCatCoefs(0) ]
                                                     , TaggedCatName = 'TagCat' if tagCatP2VVOS.numTypes() > 2 else 'Tagged'
                                                     , Data          = sigTaggingData
                                                     , RelativeCoefs = True if sigTaggingPdf.endswith('Relative') else False
                                                     , Dilutions     = [ sigPdf.dilutions(False), sigPdf.dilutions(True) ]\
                                                                       if sigTaggingPdf.startswith('tagUntag') and SSTagging else None
                                                    )

                    self._sig_tagCat_iTag = self._sigTaggingPdf.pdf()
                    self._signalComps += self._sig_tagCat_iTag


        ###################################################################################################################################
        ## build PDFs for estimated wrong-tag probabilities ##
        ######################################################

        if self._dataSets['data'] and self._dataSets['sigSWeightData'] and self._dataSets['bkgSWeightData'] and makePlots :
            # build PDFs for estimated wrong-tag probabilities
            from RooFitWrappers import HistPdf
            self._estWTagDataOS = self._dataSets['sigSWeightData'].reduce( '%s > 0' % tagCatP2VVOS.GetName() )
            if tagCatP2VVSS :
                self._estWTagDataSS = self._dataSets['sigSWeightData'].reduce( '%s > 0' % tagCatP2VVSS.GetName() )
            else :
                self._estWTagDataSS = self._dataSets['sigSWeightData']

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
            untagFracOS    = self._dataSets['data'].table(tagCatP2VVOS).getFrac('Untagged')
            untagFracSigOS = self._dataSets['sigSWeightData'].table(tagCatP2VVOS).getFrac('Untagged')
            untagFracBkgOS = self._dataSets['bkgSWeightData'].table(tagCatP2VVOS).getFrac('Untagged')
            if tagCatP2VVSS :
                untagFracSS    = self._dataSets['data'].table(tagCatP2VVSS).getFrac('Untagged')
                untagFracSigSS = self._dataSets['sigSWeightData'].table(tagCatP2VVSS).getFrac('Untagged')
                untagFracBkgSS = self._dataSets['bkgSWeightData'].table(tagCatP2VVSS).getFrac('Untagged')
            else :
                untagFracSS    = 1.
                untagFracSigSS = 1.
                untagFracBkgSS = 1.

            # plot estimated wrong-tag probabilities for signal and for background
            self._estWTagCanvOS = TCanvas( 'estWTagCanvOS', 'Estimated wrong-tag probability OS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvOS.pads( 1, 1 ) if SFit else self._estWTagCanvOS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ self._dataSets['sigSWeightData'] if SFit else self._dataSets['data']
                            , self._dataSets['sigSWeightData'], self._dataSets['bkgSWeightData'] ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigOS if SFit else untagFracOS ), 1. - untagFracSigOS, 1. - untagFracBkgOS ]
                        ) :
                plot(  pad, estWTagOS, data, self._sig_bkg_estWTagOS
                     , frameOpts  = dict( Bins = nBins, Title = estWTagOS.GetTitle() + plotTitle, Range = ( 0., 0.499999 ) )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                )
                     , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm                           )
                    )

            self._estWTagCanvSS = TCanvas( 'estWTagCanvSS', 'Estimated wrong-tag probability SS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvSS.pads( 1, 1 ) if SFit else self._estWTagCanvSS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ self._dataSets['sigSWeightData'] if SFit else self._dataSets['data']
                            , self._dataSets['sigSWeightData'], self._dataSets['bkgSWeightData'] ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigSS if SFit else untagFracSS ), 1. - untagFracSigSS, 1. - untagFracBkgSS ]
                        ) :
                plot(  pad, estWTagSS, data, self._sig_bkg_estWTagSS
                     , frameOpts  = dict( Bins = nBins, Title = estWTagSS.GetTitle() + plotTitle, Range = ( 0., 0.499999 ) )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4                                                )
                     , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm                           )
                    )


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
        self._backgroundTime = BackgroundTime(  Name = 'bkg_t', time = time, resolutionModel = self._timeResModel['model']
                                              , Efficiency = timeAcceptance if multiplyByTimeEff in [ 'all', 'background' ] else None
                                             )
        self._bkg_t = self._backgroundTime.pdf()

        if not SFit : self._backgroundComps += self._bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit or makePlots :
            bkgAngleData = self._dataSets['bkgRangeData'] if massRangeBackground else self._dataSets['bkgSWeightData']
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
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a %s angular PDF for background'\
                      % ( 'hybrid' if bkgAnglePdf == 'hybrid' else 'basis functions' if bkgAnglePdf == 'basis' else 'binned' )

                if nominalPdf or bkgAnglePdf != 'basis' :
                    # create a binned PDF for background angular shape
                    from array import array
                    baselineBin = -1
                    if not nominalPdf and transAngles :
                        nBins = [ 5, 7, 9 ]
                        cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
                        cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
                        phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )

                    elif not nominalPdf and bkgAnglePdf == 'hybrid' :
                        nBins = [ 20, 40, 20 ]
                        cpsBinBounds = array( 'd', [ -1.,                                        1.  ] )
                        cthBinBounds = array( 'd', [ -1., -0.95, -0.90, -0.85, 0.85, 0.90, 0.95, 1.  ] )
                        phiBinBounds = array( 'd', [ -pi,                                        pi  ] )
                        baselineBin = 3

                    else :
                        nBins = [ 10, 24, 5 ]
                        cpsBinBounds = array( 'd', [ -1. + 2. / 10.     * float(i) for i in range(11) ] )
                        cthBinBounds = array( 'd', [ -1. + 2. / 24.     * float(i) for i in range(25) ] )
                        #cthBinBounds = array( 'd', [ -1., -0.95, -0.90, -0.85, -0.6, -0.2, 0.2, 0.6, 0.85, 0.90, 0.95, 1. ] )
                        phiBinBounds = array( 'd', [ -pi + 2. * pi / 5. * float(i) for i in range(6)  ] )

                        #cpsBinBounds = array( 'd', [ -1., -0.5, 0., 0.5, 1. ] )
                        #cthBinBounds = array( 'd', [ -1., -0.5, 0., 0.5, 1. ] )
                        #phiBinBounds = array( 'd', [ -pi, -0.5 * pi, 0., 0.5 * pi, pi ] )

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
                                                  , Constant = False if not nominalPdf and bkgAnglePdf == 'hybrid' and bin1 != baselineBin\
                                                               else True
                                                 )\
                                          if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                                          for bin2 in range( phiNumBins ) for bin1 in range( cthNumBins ) for bin0 in range( cpsNumBins )
                                        ]

                    from RooFitWrappers import ComplementCoef
                    self._bkgAngCoefs[0] = ComplementCoef(  Name         = 'bkg_angBin_0_0_0'
                                                          , Coefficients = self._bkgAngCoefs[ 1 : ]
                                                         )

                    # create a BinnedPdf
                    from RooFitWrappers import BinnedPdf
                    self._bkg_angBins = BinnedPdf(  Name = 'bkg_angBins' if not nominalPdf and bkgAnglePdf == 'hybrid' else 'bkg_angles'
                                                  , Observables = angles
                                                  , Binnings = [ cpsBins, cthBins, phiBins ]
                                                  , Coefficients = self._bkgAngCoefs
                                                  , BinIntegralCoefs = True
                                                 )
                    self._bkg_angBins.setForceUnitIntegral(True)

                    sumWeights = 0.
                    sumBinWeights = ( cpsNumBins * cthNumBins * phiNumBins - 1 ) * [ 0. ]
                    spikesFrac = 1.
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

                        if not nominalPdf and bkgAnglePdf == 'hybrid' :
                            baseBinVal = sumBinWeights[ baselineBin - 1 ] / sumWeights / cthBins.binWidth(baselineBin)
                            spikesFrac -= baseBinVal * 2.

                        # set bin coefficient values
                        #binCoefVals = [ 0.10, 0.07, 0., 0.07, 0.10, 0.37 ]
                        #binCoefErrs = [ 0.03, 0.03, 0., 0.03, 0.03, 0.03 ]
                        for bin, ( coef, weight ) in enumerate( zip( self._bkgAngCoefs[ 1 : ], sumBinWeights ) ) :
                            binInt = weight / sumWeights
                            assert binInt >= 0.,\
                                   'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: background angular PDF coefficient \"%s\" has negative value: %f'\
                                   % (coef.GetName(), binInt)
                            if not nominalPdf and bkgAnglePdf == 'hybrid' :
                                coef.setVal( ( ( binInt - baseBinVal * cthBins.binWidth(bin + 1) ) / spikesFrac )\
                                             if bin != baselineBin - 1 else 0. )
                                #coef.setVal( binCoefVals[bin] )
                                #coef.setError( binCoefErrs[bin] )
                            else :
                                coef.setVal(binInt)

                    if not nominalPdf and bkgAnglePdf == 'hybrid' :
                        self._bkgAngBinsFrac = RealVar(  'bkgAngBinsFrac'
                                                       , Title  = 'Binned PDF fraction in angular background'
                                                       , Value  = spikesFrac
                                                       , MinMax = ( 0., 1. )
                                                      )

                if not nominalPdf and bkgAnglePdf in [ 'basis', 'hybrid' ] :
                    if bkgAnglePdf == 'basis' : nBins = [ 20, 20, 20 ]
                    # create an angular PDF with basis functions
                    angPDFIndices = [ ( 2, 0, 0 ) ]
                    #angPDFIndices = [ ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 0, 0 ), ( 2, 2, 0 ), ( 2, 2, 2 ) ]
                    #angPDFIndices = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                    #                  for YIndex1 in range( -YIndex0, YIndex0 + 1 )\
                    #                  if PIndex != 0 or YIndex0 != 0 or YIndex1 != 0 ]

                    from P2VVParameterizations.AngularPDFs import AngleBasis_AngularPdfTerms as angularPdfTerms
                    cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
                    angPdfTerms = angularPdfTerms(  Angles = self._angleFuncs.angles
                                                  , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                             , {  'Name'    : 'bkg_angCoef_%s_%s_%s'\
                                                                              % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                                , 'Value'   : 0.
                                                                , 'Error'   : 0.01
                                                                , 'MinMax'  : ( -0.2, 0.2 )
                                                                , 'Indices' : inds
                                                               }
                                                            ) for inds in angPDFIndices
                                                          )
                                                 )
                    self._bkg_angFuncs = angPdfTerms.buildSumPdf('bkg_angFuncs')

                # build total angular PDF
                if not nominalPdf and bkgAnglePdf == 'hybrid' :
                    from RooFitWrappers import SumPdf
                    self._bkg_angles = SumPdf(  Name   = 'bkg_angles'
                                              , PDFs   = [ self._bkg_angBins, self._bkg_angFuncs ]
                                              , Yields = { self._bkg_angBins.GetName() : self._bkgAngBinsFrac }
                                             )

                elif not nominalPdf and bkgAnglePdf == 'basis' :
                    self._bkg_angles = self._bkg_angFuncs

                else :
                    self._bkg_angles = self._bkg_angBins

                if bkgAngleData :
                    # fit background angular distribution
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: fitting background angular distribution'
                    print '  initial parameter values:'
                    for par in self._bkg_angles.getParameters(bkgAngleData) : par.Print()
                    self._bkg_angles.fitTo( bkgAngleData, SumW2Error = False, Save = False, **fitOpts )

            if not SFit : self._backgroundComps += self._bkg_angles

            if self._dataSets['bkgRangeData'] and self._dataSets['bkgSWeightData'] and makePlots :
                # plot background angles with S-weights
                self._bkgAnglesSWeightCanv = TCanvas( 'bkgAnglesSWeightCanv', 'Background Decay Angles with S-Weights' )
                for ( pad, obs, data, bins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesSWeightCanv.pads( 3, 2 )
                             , angles
                             , 3 * ( self._dataSets['bkgSWeightData'], )
                             , nBins
                             , [ angle.GetTitle() for angle in angles ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
                        )

                # plot background angles with side bands
                self._bkgAnglesSideBandCanv = TCanvas( 'bkgAnglesSideBandCanv', 'Background Decay Angles with Side Bands' )
                for ( pad, obs, data, bins, plotTitle, xTitle, scale )\
                      in zip(  self._bkgAnglesSideBandCanv.pads( 3, 2 )
                             , angles
                             , 3 * ( self._dataSets['bkgRangeData'], )
                             , nBins
                             , [ angle.GetTitle() for angle in angles ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                             , ( ( None, None ), ( None, None ), ( None, None ) ) # ( ( 0., 740. * 9585. / 12005. ), ( 0., 580. * 9585. / 12005. ), ( 0., 760. * 9585. / 12005. ) ) # ( ( 0., 740. ), ( 0., 580. ), ( 0., 760. ) )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle, yScale = scale
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
                        )

                # plot 2-D angular distributions
                from P2VVGeneralUtils import _P2VVPlotStash
                for angle0, angle1, data, canv, padNr in [  ( 1, 0, self._dataSets['bkgSWeightData'], self._bkgAnglesSWeightCanv,  4 )
                                                          , ( 2, 0, self._dataSets['bkgSWeightData'], self._bkgAnglesSWeightCanv,  5 )
                                                          , ( 1, 2, self._dataSets['bkgSWeightData'], self._bkgAnglesSWeightCanv,  6 )
                                                          , ( 1, 0, self._dataSets['bkgRangeData'],   self._bkgAnglesSideBandCanv, 4 )
                                                          , ( 2, 0, self._dataSets['bkgRangeData'],   self._bkgAnglesSideBandCanv, 5 )
                                                          , ( 1, 2, self._dataSets['bkgRangeData'],   self._bkgAnglesSideBandCanv, 6 )
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

        if not SFit and ( nominalPdf or multiplyByTagPdf or not condTagging ) :
            if not nominalPdf and self._iTagZeroTrick :
                # no implementation for background tagging PDF with tag = { B, Bbar, Untagged } (yet)
                pass

            else :
                # tagCat = { Untagged, TagCat1, TagCat2, ... }, tag = { B, Bbar }
                bkgTaggingData = self._dataSets['bkgRangeData'] if massRangeBackground else self._dataSets['bkgSWeightData']
                if not nominalPdf and bkgTaggingData and bkgTaggingPdf == 'histPdf' :
                    # create histogram from background data and use the (fixed) bin coefficients for the PDF
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: creating background tagging PDFs from %s data'\
                          % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )
                    from RooFitWrappers import HistPdf
                    self._bkg_tagCat_iTag = HistPdf( Name = 'bkg_tagCat_iTag', Observables = [tagCatP2VVOS, iTagOS], Data = bkgTaggingData )
                    self._backgroundComps += self._bkg_tagCat_iTag

                else :
                    # use a PDF with variable bin coefficients
                    if nominalPdf or bkgTaggingPdf.startswith('tagUntag')\
                            or ( tagCatP2VVOS.numTypes() == 2 and ( not SSTagging or tagCatP2VVSS.numTypes() == 2 ) ) :
                        # couple background tagging category coefficients to signal tagging category coefficients
                        # and assume B-Bbar asymmetry is equal for all tagged categories
                        if tagCatP2VVOS.numTypes() > 2  or ( SSTagging and tagCatP2VVSS.numTypes() > 2 ) :
                            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: background tagging PDF:\n'\
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
                                + '    * assuming equal tagging category coefficients\n'\
                                + '    * assuming absence of B-Bbar asymmetries'

                        from P2VVParameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

                    # build PDF
                    self._bkgTaggingPdf = TaggingPdf(  'tagCat_iTag'
                                                     , tagCatP2VVOS, tagCatP2VVSS if SSTagging else None
                                                     , iTagOS, iTagSS if SSTagging else None
                                                     , NamePrefix    = 'bkg'
                                                     , TaggedCatName = 'TagCat' if tagCatP2VVOS.numTypes() > 2 else 'Tagged'
                                                     , Data          = bkgTaggingData
                                                     , RelativeCoefs = True if bkgTaggingPdf.endswith('Relative') else False
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
        ## split PDF for different data samples ##
        ##########################################

        # first garbage collect
        gc.collect()

        # check number of KK mass bin parameters
        if paramKKMass in [ 'simultaneous', 'functions' ] :
            assert len(SWaveAmpVals[0]) == len(SWaveAmpVals[1]) == len(CSPValues) == len(KKMassBinBounds) - 1,\
                   'P2VV - ERROR: wrong number of KK mass bin parameters specified'
            print 'P2VV - INFO: KK mass bins: %s' % ' - '.join( '%.1f' % binEdge for binEdge in KKMassBinBounds )
        else :
            assert len(SWaveAmpVals[0]) == len(SWaveAmpVals[1]) == 0 and len(CSPValues) == 1,\
                   'P2VV - ERROR: only one S-P-wave coupling factor and no S-wave amplitude values should be specified'

        if ASParam != 'Mag2ReIm' :
            print 'P2VV - INFO: using S-P-wave coupling factors:',
            for iter, fac in enumerate(CSPValues) : print '%d: %.4f%s' % ( iter, fac, '' if iter == len(CSPValues) - 1 else ',' ),
            print

        if ( not SFit and selection in ['paper2012', 'timeEffFit'] ) or paramKKMass == 'simultaneous' :
            # categories for splitting the PDF
            splitCats = [ [ ] ]
            if not SFit :
                splitCats[0] += [ hlt1ExclB ] if selection == 'paper2012' else [ hlt1ExclB, hlt2B ] if selection == 'timeEffFit' else [ ]
            splitCats[0] += [ self._KKMassCat ] if paramKKMass == 'simultaneous' else [ ]

            # specify parameters that are different in simultaneous categories
            splitParams = [ [ ] ]
            if not SFit :
                splitParams[0].append( self._signalComps.getYield() )
                splitParams[0].append( self._backgroundComps.getYield() )
            if not SFit and selection in ['paper2012', 'timeEffFit'] :
                splitCats.append( [hlt1ExclB] if selection == 'paper2012' else [hlt1ExclB, hlt2B] if selection == 'timeEffFit' else [] )
                splitParams.append( [ par for par in self._backgroundTime.parameters() if not par.isConstant() ] )
            #if selection == 'paper2012' :
            #    if SFit :
            #        splitCats.append( [hlt1ExclB] )
            #        splitParams.append( [ par for par in self._timeResModel['binHeights'] ] )
            #    else :
            #        splitParams[-1] += [ par for par in self._timeResModel['binHeights'] ]
            if paramKKMass == 'simultaneous' :
                splitCats.append( [ self._KKMassCat ] )
                splitParams.append( [ ] )
                for amp in self._amplitudes.parameters() :
                    if any( name in amp.GetName() for name in [ 'AS', 'A_S', 'fS', 'f_S', 'C_SP' ] ) : splitParams[-1].append(amp)
                if not SFit :
                    splitParams[-1] += [ par for par in self._backgroundBMass.parameters() if not par.isConstant() ]

            # build simultaneous PDF
            from RooFitWrappers import SimultaneousPdf
            self._simulPdf = SimultaneousPdf(  self._fullPdf.GetName() + '_simul'
                                             , MasterPdf       = self._fullPdf
                                             , SplitCategories = splitCats
                                             , SplitParameters = splitParams
                                            )

            massPars = self._sWeightMassPdf.getVariables()
            splitCatPars = self._simulPdf.getVariables()
            if not SFit :
                # set values for splitted yields
                splitCatIter = self._simulPdf.indexCat().typeIterator()
                splitCatState = splitCatIter.Next()
                while splitCatState :
                    sigYield = getSplitPar( 'N_sigMass' if SFit else 'N_signal', splitCatState.GetName(), splitCatPars )
                    bkgYield = getSplitPar( 'N_bkgMass' if SFit else 'N_bkg',    splitCatState.GetName(), splitCatPars )

                    if not sigYield in massPars or not bkgYield in massPars :
                        if splitCat.isFundamental() :
                            selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
                        else :
                            splitCat.setLabel( splitCatState.GetName() )
                            selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
                        nEv    = self._dataSets['data'].sumEntries()
                        nEvBin = self._dataSets['data'].sumEntries(selStr)

                        sigYield.setVal( sigYield.getVal() * nEvBin / nEv )
                        sigYield.setError( sqrt( sigYield.getVal() ) )
                        sigYield.setMin(0.)
                        sigYield.setMax(nEvBin)
                        bkgYield.setVal( bkgYield.getVal() * nEvBin / nEv )
                        bkgYield.setError( sqrt( bkgYield.getVal() ) )
                        bkgYield.setMin(0.)
                        bkgYield.setMax(nEvBin)

                    splitCatState = splitCatIter.Next()

            if not SFit and selection in ['paper2012', 'timeEffFit'] :
                # set values for background parameters in different trigger samples
                splitCatIter = hlt1ExclB.typeIterator()
                splitCatState = splitCatIter.Next()
                while splitCatState :
                    if splitCatState.getVal() != 0 :
                        # get names of background time parameters (FIXME: isn't there a simpler way to do this?)
                        mlFracName = self._backgroundTime.pdf().coefList().at(0).GetName()
                        mlTauName = ''
                        mlPdfPars = self._backgroundTime.pdf().pdfList().at(0).getVariables()
                        for par in self._backgroundTime.parameters() :
                            mlTau = mlPdfPars.find( par.GetName() )
                            if mlTau : mlTauName = mlTau.GetName()
                        assert mlTauName

                        # "remove" medium-lived lifetime PDF for HLT1 exclusively biased events
                        mlFrac = getSplitPar( mlFracName, splitCatState.GetName(), splitCatPars )
                        mlTau  = getSplitPar( mlTauName,  splitCatState.GetName(), splitCatPars )
                        mlFrac.setVal(0.)
                        mlTau.setVal(0.)
                        mlFrac.setConstant(True)
                        mlTau.setConstant(True)

                    splitCatState = splitCatIter.Next()

            #if selection == 'paper2012' :
            #    accFile = TFile.Open(timeEffHistFile)
            #    accHist = accFile.Get(timeEffHistExclBName)
            #    exclBName = hlt1ExclB.lookupType(1).GetName()
            #    for binIt in range( accHist.GetNbinsX() ) :
            #        parName = self._timeResModel['binHeights'][binIt]
            #        binHeight = getSplitPar( parName, exclBName, splitCatPars )
            #        binHeight.setVal( accHist.GetBinContent( binIt + 1 ) )
            #    accFile.Close()

            if paramKKMass == 'simultaneous' :
                # set values for splitted amplitudes
                splitCatIter = self._KKMassCat.typeIterator()
                splitCatState = splitCatIter.Next()
                while splitCatState :
                    if ASParam != 'Mag2ReIm' :
                        # S-P-wave coupling factors
                        C_SP = getSplitPar( 'C_SP', splitCatState.GetName(), splitCatPars )
                        C_SP.setVal( CSPValues[ splitCatState.getVal() ] )
                        C_SP.setConstant(True)

                    if nominalPdf or ( amplitudeParam == 'bank' and ASParam != 'ReIm' )\
                            or ( amplitudeParam == 'phasesSWaveFrac' and ASParam == 'deltaPerp' ) :
                        # amplitude parameterization with delta_S-delta_perp
                        if amplitudeParam == 'phasesSWaveFrac' :
                            f_S = getSplitPar( 'f_S', splitCatState.GetName(), splitCatPars )
                        else :
                            ASOddMag2 = getSplitPar( 'ASOddMag2', splitCatState.GetName(), splitCatPars )
                        ASOddPhase = getSplitPar( 'ASOddPhase', splitCatState.GetName(), splitCatPars )

                        if amplitudeParam == 'phasesSWaveFrac' :
                            f_S.setVal(   SWaveAmpVals[0][ splitCatState.getVal() ][0] )
                            f_S.setError( SWaveAmpVals[0][ splitCatState.getVal() ][1] )
                        else :
                            ASOddMag2.setVal(   SWaveAmpVals[0][ splitCatState.getVal() ][0] )
                            ASOddMag2.setError( SWaveAmpVals[0][ splitCatState.getVal() ][1] )
                            if ASOddMag2.getMax() < 5. : ASOddMag2.setMax(5.)
                        ASOddPhase.setVal(   SWaveAmpVals[1][ splitCatState.getVal() ][0] )
                        ASOddPhase.setError( SWaveAmpVals[1][ splitCatState.getVal() ][1] )

                    elif amplitudeParam == 'phases' and ASParam in [ 'ReIm', 'Mag2ReIm' ] :
                        # amplitude parameterization with Re(A_S) and Im(A_S)
                        if ASParam == 'Mag2ReIm' :
                            ASMag2 = getSplitPar( 'ASMag2', splitCatState.GetName(), splitCatPars )
                        ReAS = getSplitPar( 'ReAS', splitCatState.GetName(), splitCatPars )
                        ImAS = getSplitPar( 'ImAS', splitCatState.GetName(), splitCatPars )

                        if ASParam == 'Mag2ReIm' :
                            ASMag2.setVal( SWaveAmpVals[0][ splitCatState.getVal() ][0]**2\
                                          + SWaveAmpVals[1][ splitCatState.getVal() ][0]**2 )
                        ReAS.setVal(   SWaveAmpVals[0][ splitCatState.getVal() ][0] )
                        ReAS.setError( SWaveAmpVals[0][ splitCatState.getVal() ][1] )
                        ImAS.setVal(   SWaveAmpVals[1][ splitCatState.getVal() ][0] )
                        ImAS.setError( SWaveAmpVals[1][ splitCatState.getVal() ][1] )

                    splitCatState = splitCatIter.Next()

        assert not pdfConfig, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: superfluous arguments found: %s' % pdfConfig
        PdfBuilder.__init__( self, self._simulPdf if paramKKMass == 'simultaneous' else self._fullPdf, observables, { } )
        print 120 * '='

        # garbage collect
        gc.collect()

    def __getitem__( self, kw ) :
        return self._dataSets[kw] if kw in [ 'data', 'sigSWeightData', 'bkgSWeightData', 'sigRangeData', 'bkgRangeData' ]\
               else PdfBuilder.__getitem__( self, kw )
