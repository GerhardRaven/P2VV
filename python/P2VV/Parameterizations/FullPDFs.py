###########################################################################################################################################
## P2VV.Parameterizations.FullPDFs: Parameterizations of complete PDFs that are used in an analysis                                      ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

###########################################################################################################################################
## PDF configuration ##
#######################

# PDF configuration base class
class PdfConfiguration( dict ) :
    def __init__( self, parameters = None, **kwargs ) :
        self._parameters = { }
        if parameters != None : self.addParameters(parameters)
        self['parameters'] = self._parameters
        dict.__init__(self,**kwargs)

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
            if par.ClassName() != 'RooRealVar' : continue
            self.addParameter( par.GetName(), (  par.getVal()
                                               , par.getError()
                                               , par.getMin()
                                               , par.getMax()
                                               , False if par.isConstant() else True
                                              )
                             )

    def setParametersInPdf( self, pdf ) :
        for par in pdf.getVariables() :
            if not par.GetName() in self._parameters.keys() : continue
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
              parFloat = ( line[5] == 'True' )
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
            if ( floating and not parVals[4] ) or ( not floating and parVals[4] ) : continue

            cont += ( '  {0:<%s}   {1:<+14.8g}   {2:<13.8g}   {3:<+14.8g}   {4:<+14.8g}   {5:<}\n' % maxLenName )\
                      .format( parName, parVals[0], parVals[1], parVals[2], parVals[3], 'True' if parVals[4] else 'False' )
            numPars += 1

        cont += '# ' + '-' * (79 + maxLenName) + '\n'

        # write content to file
        parFile.write(cont)
        parFile.close()

        print 'P2VV - INFO: PdfConfiguration.writeParametersToFile: %d parameter%s written to file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )


# B_s^0 -> J/psi phi analysis of 1 fb^-1 2011 data (paper) configuration
class Bs2Jpsiphi_2011Analysis( PdfConfiguration ) :
    def __init__( self ) :
        from math import pi

        # job parameters
        self['dataSample']  = ''
        self['selection']   = 'paper2012'
        self['makePlots']   = False
        self['AllTagPlots'] = False
        self['SFit']        = True
        self['blind']       = { }

        self['numEvents'] = 54755

        self['dataSet'] = None
        self['obsDict'] = dict(  mass      = ( 'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.   )
                               , mumuMass  = ( 'mdau1',                'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.   )
                               , KKMass    = ( 'mdau2',                'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,   990.,   1050.   )
                               , KKMassCat = ( 'KKMassCat',            'KK mass category',        { }                                   )
                               , time      = ( 'time',                 'Decay time',              'ps',        1.5,    0.3,     14.     )
                               , timeRes   = ( 'sigmat',               '#sigma(t)',               'ps',        0.01,   0.0001,  0.12    )
                               , cpsi      = ( 'helcosthetaK',         'cos(#theta_{K})',         '',          0.,    -1.,     +1.      )
                               , ctheta    = ( 'helcosthetaL',         'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.      )
                               , phi       = ( 'helphi',               '#phi_{h}',                'rad',       0.,    -pi,     +pi      )
                               , wTagOS    = ( 'tagomega_os',          'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001 )
                               , wTagSS    = ( 'tagomega_ss',          'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001 )
                               , iTagOS    = ( 'iTagOS',               'OS flavour tag',          { 'B' : +1, 'Bbar' : -1 }             )
                               , iTagSS    = ( 'iTagSS',               'SS flavour tag',          { 'B' : +1, 'Bbar' : -1 }             )
                               , tagCatOS  = ( 'tagCatP2VVOS',         'OS flavour tag',          { 'Untagged' : 0, 'Tagged' : 1 }      )
                               , tagCatSS  = ( 'tagCatP2VVSS',         'SS flavour tag',          { 'Untagged' : 0, 'Tagged' : 1 }      )
                               , hlt1ExclB = ( 'hlt1_excl_biased_dec', 'HLT1 excl. B.',           { 'exclB' : 1, 'notExclB' : 0 }       )
                               , hlt2B     = ( 'hlt2_biased',          'HLT2 B.',                 { 'B'     : 1, 'notB'     : 0 }       )
                               , hlt2UB    = ( 'hlt2_unbiased',        'HLT2 UB.',                { 'UB'    : 1, 'notUB'    : 0 }       )
                              )

        self['timeEffHistFile']      = 'data/Bs_HltPropertimeAcceptance_Data-20120816.root'
        self['timeEffHistUBName']    =\
            'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
        self['timeEffHistExclBName'] =\
            'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
        self['angEffMomentsFile']    = 'data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

        # fit options
        self['fitOptions'] = dict( NumCPU = 1, Optimize = 1, Timer = True, Minimizer = 'Minuit2' )

        # PDF parameters
        self['tagCatsOS'] = [ ]
        self['tagCatsSS'] = [ ]

        # PDF options
        self['transversityAngles']   = False
        self['angularRanges']        = dict( cpsi = [ ], ctheta = [ ], phi = [ ] )
        self['sigMassModel']         = ''
        self['bkgMassModel']         = ''
        self['bkgAnglePdf']          = 'hybrid'
        self['sigTaggingPdf']        = 'tagUntag'
        self['bkgTaggingPdf']        = 'tagUntagRelative'
        self['multiplyByTagPdf']     = False
        self['multiplyByTimeEff']    = 'signal'
        self['timeEffType']          = 'paper2012'
        self['timeEffParameters']    = dict( Spline= False,  SmoothSpline = 0.8 )
        self['multiplyByAngEff']     = 'weights'
        self['parameterizeKKMass']   = 'simultaneous'
        self['ambiguityParameters']  = False
        self['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
        self['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.067, 0.029), (0.008, 0.011), (0.016, 0.011), (0.055, 0.026), (0.17,  0.04) ]
                                        , [ (1.3,  0.7 ), (0.77,  0.28 ), (0.50,  0.47 ), (-0.51, 0.25 ), (-0.46, 0.21 ), (-0.65, 0.20) ] )
        self['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

        self['sameSideTagging']    = True
        self['conditionalTagging'] = True
        self['continuousEstWTag']  = True
        self['numEstWTagBins']     = 50
        self['constrainTagging']   = 'constrain'

        self['timeResType']           = 'eventNoMean'
        self['numTimeResBins']        = 40
        self['constrainTimeResScale'] = 'constrain'

        self['signalFraction'] = 0.504

        self['amplitudeParam'] = 'phasesSWaveFrac'
        self['ASParam']        = 'deltaPerp'
        self['AparParam']      = 'phase'

        self['constrainDeltaM'] = 'constrain'

        self['lambdaCPParam'] = 'lambPhi'

        self['angleNames'] = (  ( 'helcosthetaK', 'cos#kern[0.1]{#theta_{K}}'   )
                              , ( 'helcosthetaL', 'cos#kern[0.1]{#theta_{#mu}}' )
                              , ( 'helphi',       '#varphi_{h} [rad]'           )
                             )

        self['numBMassBins'] = [ 70, 40, 20, 20, 20 ]
        self['numTimeBins']  = 30
        self['numAngleBins'] = ( 10, 24, 5 )

        # initialize PdfConfiguration object
        PdfConfiguration.__init__( self )


###########################################################################################################################################
## PDF builders ##
##################

# PDF builder base class
class PdfBuilder ( object ) :
    def __init__( self, pdf, observables, parameters ) :
        self._pdf         = pdf
        self._observables = observables
        self._parameters  = parameters

    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf(self)         : return self._pdf
    def observables(self) : return self._observables
    def parameters(self)  : return self._parameters


# B_s^0 -> J/psi phi PDF builder
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
        dataSample  = pdfConfig.pop('dataSample')    # '' / 'Summer2011'
        selection   = pdfConfig.pop('selection')     # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'timeEffFit'
        SFit        = pdfConfig.pop('SFit')
        blind       = pdfConfig.pop('blind')
        makePlots   = pdfConfig.pop('makePlots')
        AllTagPlots = pdfConfig.pop('AllTagPlots')

        numEvents = pdfConfig.pop('numEvents')

        dataSet = pdfConfig.pop('dataSet')
        obsDict = pdfConfig.pop('obsDict')

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
        sigMassModel      = pdfConfig.pop('sigMassModel')           # '' / 'LP2011' / 'Gauss' / 'DoubleGauss' / 'DoubleGaussDiag' / 'CB' / 'DoubleCB' / 'box' / 'boxFixed'
        bkgMassModel      = pdfConfig.pop('bkgMassModel')           # '' / 'LP2011' / 'linear'
        bkgAnglePdf       = pdfConfig.pop('bkgAnglePdf')            # '' / 'histPdf' / 'binned' / 'basis' / 'hybrid'
        sigTaggingPdf     = pdfConfig.pop('sigTaggingPdf')          # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        bkgTaggingPdf     = pdfConfig.pop('bkgTaggingPdf')          # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        multiplyByTagPdf  = pdfConfig.pop('multiplyByTagPdf')
        multiplyByTimeEff = pdfConfig.pop('multiplyByTimeEff')      # '' / 'all' / 'signal'
        timeEffType       = pdfConfig.pop('timeEffType')            # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'fit'
        timeEffParameters = pdfConfig.pop('timeEffParameters')
        multiplyByAngEff  = pdfConfig.pop('multiplyByAngEff')       # '' / 'weights' / 'basis012' / 'basis012Plus' / 'basis012Thetal' / 'basis0123' / 'basis01234' / 'basisSig3' / 'basisSig4'
        paramKKMass       = pdfConfig.pop('parameterizeKKMass')     # '' / 'parameters' / 'simultaneous'
        numBMassBins      = pdfConfig.pop('numBMassBins')
        ambiguityPars     = pdfConfig.pop('ambiguityParameters')
        KKMassBinBounds   = pdfConfig.pop('KKMassBinBounds')
        SWaveAmpVals      = pdfConfig.pop('SWaveAmplitudeValues')
        CSPValues         = pdfConfig.pop('CSPValues')

        if not SWaveAmpVals : SWaveAmpVals = ( [ ], [ ] )

        SSTagging        = pdfConfig.pop('sameSideTagging')
        condTagging      = pdfConfig.pop('conditionalTagging')
        contEstWTag      = pdfConfig.pop('continuousEstWTag')
        numEstWTagBins   = pdfConfig.pop('numEstWTagBins')
        constrainTagging = pdfConfig.pop('constrainTagging')      # '' / 'constrain' / 'fixed'
        condTagging = True if contEstWTag else condTagging

        timeResType     = pdfConfig.pop('timeResType')              # '' / 'event' / 'eventNoMean' / 'eventConstMean' / '3Gauss'
        numTimeResBins  = pdfConfig.pop('numTimeResBins')
        constrTResScale = pdfConfig.pop('constrainTimeResScale')    # '' / 'constrain' / 'fixed'

        sigFrac = pdfConfig.pop('signalFraction')

        amplitudeParam = pdfConfig.pop('amplitudeParam')    # 'phases' / 'phasesSWaveFrac' / 'ReIm' / 'bank'
        ASParam        = pdfConfig.pop('ASParam')           # 'delta0' / 'deltaPerp' / 'ReIm' / 'Mag2ReIm' / 'Mag2ReImPerp'
        AparParam      = pdfConfig.pop('AparParam')         # 'phase' / 'ReIm' / 'Mag2ReIm' / 'cos' / 'real'

        if not angRanges : angRanges = dict( cpsi = [ ], ctheta = [ ], phi = [ ] )
        if 'cpsi'   not in angRanges : angRanges['cpsi']   = [ ]
        if 'ctheta' not in angRanges : angRanges['ctheta'] = [ ]
        if 'phi'    not in angRanges : angRanges['phi']    = [ ]

        if paramKKMass and ambiguityPars :
            if ( amplitudeParam == 'phasesSWaveFrac' and ASParam == 'deltaPerp' ) or ( amplitudeParam == 'bank' and ASParam != 'ReIm' ) :
                for phaseIter, phase in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][phaseIter] = ( pi - phase[0], phase[1] )
            elif amplitudeParam == 'phases' and ASParam in [ 'Mag2ReIm', 'ReIm' ] :
                for ImIter, Im in enumerate( SWaveAmpVals[1] ) : SWaveAmpVals[1][ImIter] = ( -Im[0], Im[1] )

        constrainDeltaM = pdfConfig.pop('constrainDeltaM')    # '' / 'constrain' / 'fixed'

        lambdaCPParam = pdfConfig.pop('lambdaCPParam')    # 'ReIm' / 'lambSqPhi' / 'lambPhi' / 'lambPhi_CPVDecay' / 'lambPhiRel_CPVDecay'

        # check KK mass bin parameters
        assert obsDict['KKMass'][4] == KKMassBinBounds[0] and obsDict['KKMass'][5] == KKMassBinBounds[-1]\
               , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: KK mass range in "KKMassBinBounds" is not the same as in "obsDict"'

        if paramKKMass == 'parameters' or 'simultaneous' :
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
        from P2VV.RooFitWrappers import RooObject, ConstVar, RealVar, Category, FormulaVar
        ws = RooObject().ws()

        if makePlots :
            # import plotting tools
            from P2VV.GeneralUtils import plot
            from ROOT import TCanvas, kBlue, kRed, kGreen, kSolid, kDashed, kFullCircle, TPaveText
            from P2VV.GeneralUtils import _P2VVPlotStash

            LHCbLabel = TPaveText( 0.24, 0.81, 0.37, 0.89, 'BRNDC' )
            LHCbLabel.AddText('LHCb')
            LHCbLabel.SetFillColor(0)
            LHCbLabel.SetTextAlign(12)
            LHCbLabel.SetTextSize(0.072)
            LHCbLabel.SetBorderSize(0)
            _P2VVPlotStash.append(LHCbLabel)


        ###################################################################################################################################
        ## create variables ##
        ######################

        # get observables
        observables = { }
        for name in obsDict.keys():
            if dataSet :
                if   ws.var( obsDict[name][0] ) : observables[name] = RealVar(  obsDict[name][0] )
                elif ws.cat( obsDict[name][0] ) : observables[name] = Category( obsDict[name][0] )
                else : raise RuntimeError( 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: variable "%s" not in workspace' % name )
            else :
                if name.startswith('tagCat') :
                    observables[name] = None
                    continue

                if type( obsDict[name][2] ) != dict :
                    observables[name] = RealVar( obsDict[name][0], Title = obsDict[name][1], Unit = obsDict[name][2]
                                                , Value = obsDict[name][3], MinMax = ( obsDict[name][4], obsDict[name][5] )
                                               )
                else :
                    states = obsDict[name][2] if name != 'KKMassCat' else\
                             dict( [ ( 'bin%d' % it, it ) for it in range( len(KKMassBinBounds) - 1 ) ] )
                    observables[name] = Category( obsDict[name][0], Title = obsDict[name][1], States = states )

        # set observable properties
        for obs in observables.itervalues() : obs.setObservable(True)

        observables['time'].setRange( 'Bulk', ( None, 5. ) )
        observables['time'].setBins(numTimeBins)
        observables['timeRes'].setBins(numTimeResBins)
        observables['timeRes'].setBins( numTimeResBins, 'cache' )
        observables['mass'].setBins(numBMassBins[0])
        observables['mass'].setRanges( dict(  LeftSideBand  = ( 5200., 5320. )
                                            , Signal        = ( 5320., 5420. )
                                            , RightSideBand = ( 5420., 5550. )
                                            , PeakBkg       = ( 5390., 5440. )
                                           )
                                     )

        for ran in angRanges['cpsi'  ] : observables['cpsi'].setRange(   ran[0], ran[ 1 : 3 ] )
        for ran in angRanges['ctheta'] : observables['ctheta'].setRange( ran[0], ran[ 1 : 3 ] )
        for ran in angRanges['phi'   ] : observables['phi'].setRange(    ran[0], ran[ 1 : 3 ] )

        if dataSet and 'KKMassCat' in obsDict:
            # get KK mass binning
            assert observables['KKMass'].hasBinning('KKMassBinning')\
                   , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: KK mass observable does not have a binning named "KKMassBinning"'
            self._KKMassBinning = observables['KKMass'].getBinning('KKMassBinning')

            assert self._KKMassBinning.numBins() == observables['KKMassCat'].numTypes() == len(KKMassBinBounds) - 1\
                   , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: different numbers of bins in KK mass category and/or binning and specified bin bounds'
            for it in range( self._KKMassBinning.numBins() ) :
                assert self._KKMassBinning.binLow(it) == KKMassBinBounds[it]\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: different numbers of bins in KK mass binning and specified bin bounds'

        else :
            # create KK mass binning
            from array import array
            KKMassBinsArray = array( 'd', KKMassBinBounds )

            from ROOT import RooBinning
            self._KKMassBinning = RooBinning( len(KKMassBinBounds) - 1, KKMassBinsArray, 'KKMassBinning' )
            observables['KKMass'].setBinning( self._KKMassBinning, 'KKMassBinning' )

        # PDF observables set
        obsSetP2VV = [ observables[name] for name in [ 'time', 'cpsi', 'ctheta', 'phi' ] ]
        if multiplyByTagPdf or not condTagging :
            obsSetP2VV.append( observables['iTagOS'] )
            if SSTagging : obsSetP2VV.append( observables['iTagSS'] )
        if not SFit :
            obsSetP2VV.append( observables['mass'] )


        ###################################################################################################################################
        ## build tagging categories ##
        ##############################

        print 120 * '='
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building tagging categories'

        if dataSet :
            # get category bins
            assert observables['wTagOS'].hasBinning('tagCats'),\
                   'P2VV - ERROR:  Bs2Jpsiphi_PdfBuilder: binning "tagCats" not found for OS estimated wrong-tag probability'
            assert observables['wTagSS'].hasBinning('tagCats'),\
                   'P2VV - ERROR:  Bs2Jpsiphi_PdfBuilder: binning "tagCats" not found for SS estimated wrong-tag probability'
            etaBinsOS = observables['wTagOS'].getBinning('tagCats')
            etaBinsSS = observables['wTagSS'].getBinning('tagCats')

            # get bin parameters from data for OS and SS
            tagBins = [ ]
            for wTag, cat, bins, isSS in zip(  ( obsDict['wTagOS'][0],                        obsDict['wTagSS'][0]                        )
                                             , ( observables['tagCatOS'],                     observables['tagCatSS']                     )
                                             , ( observables['wTagOS'].getBinning('tagCats'), observables['wTagSS'].getBinning('tagCats') )
                                             , ( False,                                       True                                        )
                                            ) :
                for it in range( bins.numBins() ) :
                    assert cat.isValidIndex(it), 'P2VV - ERROR: no bin %d found for tagging category "%s" ' % ( it, cat.GetName() )
                binPars = [ ( cat.lookupType(it).GetName(), it, bins.binHigh( bins.numBins()-it-1 ) ) for it in range( bins.numBins() ) ]

                from P2VV.Parameterizations.FlavourTagging import getTagCatParamsFromData as getTagParams
                tagBins.append( getTagParams( dataSet, estWTagName = wTag, tagCats = binPars, numSigmas = 1., SameSide = isSS ) )

            tagCatsOS = tagBins[0]
            tagCatsSS = tagBins[1]

        else :
            assert tagCatsOS and tagCatsSS, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no tagging category binnings found'

        # build tagging categories
        from P2VV.Parameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
        self._tagCatsOS = TaggingCategories(  tagCat = observables['tagCatOS'] if observables['tagCatOS'] else obsDict['tagCatOS'][0]
                                            , estWTag = observables['wTagOS'] if contEstWTag else None
                                            , estWTagName = obsDict['wTagOS'][0], TagCats = tagCatsOS, SameSide = False
                                            , wTagP0Constraint = constrainTagging, wTagP1Constraint = constrainTagging
                                           )
        self._tagCatsSS = TaggingCategories(  tagCat = observables['tagCatSS'] if observables['tagCatSS'] else obsDict['tagCatSS'][0]
                                            , estWTag = observables['wTagSS'] if contEstWTag else None
                                            , estWTagName = obsDict['wTagSS'][0], TagCats = tagCatsSS, SameSide = True
                                            , wTagP0Constraint = constrainTagging, wTagP1Constraint = constrainTagging
                                           )

        observables['tagCatOS'] = self._tagCatsOS['tagCat']
        observables['tagCatSS'] = self._tagCatsSS['tagCat']

        if multiplyByTagPdf or not condTagging :
            obsSetP2VV.append( observables['tagCatOS'] )
            if SSTagging : obsSetP2VV.append( observables['tagCatSS'] )

        if condTagging :
            self._tagCatsOS.addConditional( observables['tagCatOS'] )
            self._tagCatsOS.addConditional( observables['iTagOS']   )
            self._tagCatsSS.addConditional( observables['tagCatSS'] )
            self._tagCatsSS.addConditional( observables['iTagSS']   )


        ###################################################################################################################################
        ## initialize PDF component objects ##
        ######################################

        nEvents     = dataSet.sumEntries() if dataSet else numEvents
        nSignal     = nEvents * sigFrac
        nBackground = nEvents * ( 1. - sigFrac )

        from P2VV.RooFitWrappers import Component
        self._signalComps     = Component( 'signal', [ ], Yield = ( nSignal,     0., nEvents ) )
        self._backgroundComps = Component( 'bkg'   , [ ], Yield = ( nBackground, 0., nEvents ) )


        ###################################################################################################################################
        ## build B mass PDFs ##
        #######################

        if not SFit :
            # build the signal and background mass PDFs
            sigMassArgs = dict( Name = 'sig_m', mass = observables['mass'] )
            if sigMassModel.startswith('box') :
                from P2VV.Parameterizations.MassPDFs import Box_Signal_Mass as SignalBMass
                if sigMassModel.endswith('Fixed') :
                    boxWidth = 0.5 * ( observables['mass'].getMin('RightSideBand') - observables['mass'].getMax('LeftSideBand') )
                    boxMean  = observables['mass'].getMax('LeftSideBand') + boxWidth
                    sigMassArgs['m_sig_mean']  = dict( Value = boxMean,  Constant = True )
                    sigMassArgs['m_sig_width'] = dict( Value = boxWidth, Constant = True, MinMax = ( 0.1, 2. * boxWidth ) )

            elif sigMassModel.startswith('Gauss') :
                from P2VV.Parameterizations.MassPDFs import Gauss_Signal_Mass as SignalBMass

            elif sigMassModel.startswith('DoubleGauss') :
                from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as SignalBMass
                if sigMassModel.endswith('Diag') :
                    sigMassArgs['TransformWidthPars'] = dict(  m_sig_frac    = ( +0.033129, -0.008339, -0.007473 )
                                                             , m_sig_sigma_1 = ( +0.115025, -0.067412, +0.000953 )
                                                             , m_sig_sigma_2 = ( +0.756560, +0.010614, +0.000182 )
                                                            )

            elif sigMassModel.startswith('CB') :
                from P2VV.Parameterizations.MassPDFs import CB_Signal_Mass as SignalBMass

            elif sigMassModel.startswith('DoubleCB') :
                from P2VV.Parameterizations.MassPDFs import DoubleCB_Signal_Mass as SignalBMass

            else :
                from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass

            bkgMassArgs = dict( Name = 'bkg_m', mass = observables['mass'] )
            if bkgMassModel.startswith('linear') :
                from P2VV.Parameterizations.MassPDFs import Linear_Background_Mass as BackgroundBMass
                if bkgMassModel.endswith('Constant') : bkgMassArgs['Constant'] = True
            else :
                from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass

            self._signalBMass     = SignalBMass(     **sigMassArgs )
            self._backgroundBMass = BackgroundBMass( **bkgMassArgs )

            from P2VV.RooFitWrappers import buildPdf
            self._signalComps     += self._signalBMass.pdf()
            self._backgroundComps += self._backgroundBMass.pdf()


        ###################################################################################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################

        print 120 * '='
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building PDFs'

        # angular functions
        if transAngles : from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
        else :           from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles     as AngleFuncs
        self._angleFuncs = AngleFuncs( cpsi = observables['cpsi'], ctheta = observables['ctheta'], phi = observables['phi'] )

        # transversity amplitudes
        commonArgs = dict( AmbiguityParameters = ambiguityPars )
        if paramKKMass == 'parameters' :
            commonArgs[ 'KKMassCategory' ] = observables['KKMassCat']
        if not ASParam.startswith('Mag2ReIm') :
            if len(CSPValues) == 1 :
                commonArgs['C_SP'] = CSPValues[0]
            elif paramKKMass == 'parameters' :
                for it, val in enumerate(CSPValues) : commonArgs[ 'C_SP_bin%d' %it ] = val

        if amplitudeParam == 'phasesSWaveFrac' :
            from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = ASParam, AparParameterization = AparParam, **commonArgs )

        elif amplitudeParam == 'phases' :
            from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = ASParam, **commonArgs )

        elif amplitudeParam == 'bank' :
            from P2VV.Parameterizations.DecayAmplitudes import JpsiVBank_AmplitudeSet as Amplitudes
            self._amplitudes = Amplitudes( ASParameterization = ASParam, AparParameterization = AparParam, **commonArgs )

        else :
            raise RuntimeError('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no valid amplitude parameterization specified')

        # B decay time
        from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        dGammaVar = dict( Name = 'dGamma' )
        if blind and 'dGamma' in blind :
            if blind['dGamma'] : dGammaVar['Blind'] = blind['dGamma']
            else :               dGammaVar['Blind'] = ( 'UnblindUniform', 'BsDGs2013EPS', 0.02 )
        self._lifetimeParams = LifetimeParams( dGamma = dGammaVar, dMConstraint = constrainDeltaM )
        if ambiguityPars : self._lifetimeParams['dGamma'].setVal( -self._lifetimeParams['dGamma'].getVal() )

        if timeResType.startswith('event') :
            from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
            timeResArgs = dict(  time = observables['time']
                               , timeResSigma = observables['timeRes']
                               , Cache = multiplyByTimeEff not in [ 'all', 'signal', 'background' ]  # make sure we do not 'double cache'
                              )
            if 'doublegauss' in timeResType.lower():
                constant = 'constant' in timeResType.lower()
                timeResArgs['timeResComb'] = dict(Name = 'timeResComb', Value = 1.4918, Error = 4.08e-03, MinMax = ( 0.1, 5. ), Constant = constant)
                timeResArgs['timeResSigmaSF2'] = dict( Name = 'timeResSigmaSF2', Value = 6.0074, Error = 1.89e-01, MinMax = (1, 10), Constant = constant)
                timeResArgs['timeResSigmaFrac2'] = dict( Name = 'timeResSigmaFrac2', Value = 1.5818e-02, Error = 1.07e-03, MinMax = (0.001, 0.999), Constant = constant)
                timeResArgs['Covariance'] = { ('timeResComb',       'timeResComb'       ) :  1.663e-05,
                                              ('timeResComb',       'timeResSigmaFrac2' ) :  1.322e-06,
                                              ('timeResComb',       'timeResSigmaSF2'   ) :  0.0001297,
                                              ('timeResSigmaFrac2', 'timeResSigmaFrac2' ) :  1.146e-06,
                                              ('timeResSigmaFrac2', 'timeResSigmaSF2'   ) : -0.0001486,
                                              ('timeResSigmaSF2',   'timeResSigmaSF2'   ) :  0.03556   }
                timeResArgs['nGauss'] = 2
                if 'constmean' in timeResType.lower() :
                    timeResArgs['timeResMean'] = dict(Value = -4.0735e-03, Error = 1.33e-04)
                    timeResArgs['timeResMeanConstraint'] = 'constrain'
                elif 'fixedmean' in timeResType.lower() :
                    timeResArgs['timeResMean'] = dict(Value = -4.0735e-03, Error = 1.33e-04)
                    timeResArgs['timeResMeanConstraint'] = 'fixed'
                else:
                    timeResArgs['timeResMean'] = dict(Value = 0, Error = 0)
                    timeResArgs['timeResMeanConstraint'] = 'fixed'
            elif 'nomean' in timeResType.lower() :
                timeResArgs['timeResMean']   = ConstVar( Name = 'timeResMean',   Value = 0. )
                timeResArgs['timeResMeanSF'] = ConstVar( Name = 'timeResMeanSF', Value = 1. )
                timeResArgs['timeResSFConstraint'] = constrTResScale
            elif 'constmean' in timeResType.lower() :
                timeResArgs['timeResMean']   = dict( Value = -0.01, Error = 0.005 )
                timeResArgs['timeResMeanSF'] = ConstVar( Name = 'timeResMeanSF', Value = 1. )
                timeResArgs['timeResMeanConstraint'] = constrTResScale
                timeResArgs['timeResSFConstraint'] = constrTResScale
            elif 'stlinear' in timeResType.lower():
                timeResArgs['timeResMeanConstraint'] = 'constrain'
                timeResArgs['timeResSigmaSF'] = dict(Name = 'timeResSigmaSF', Value = 1.253, Error = 0.014, MinMax = (0.1, 5 ), Constant = True)
                timeResArgs['timeResSigmaOffset'] = dict(Name = 'timeResSigmaOffset', Value = 0.0153,
                                                         Error = 0.00011, Constant = True)
                covariance = {('timeResSigmaOffset', 'timeResSigmaOffset'): 1.301e-08,
                              ('timeResSigmaOffset', 'timeResSigmaSF'): 5.545e-07,
                              ('timeResSigmaSF', 'timeResSigmaSF'): 0.0002012}
                timeResArgs['Covariance'] = covariance
                timeResArgs['timeResSFModel'] = 'linear'
            elif 'stquad' in timeResType.lower():
                timeResArgs['timeResMeanConstraint'] = 'constrain'
                timeResArgs['timeResSigmaOffset'] = dict( Name = 'timeResSigmaOffset', Value = 0.0159, Error = 0.000148, MinMax = (0.001, 0.1))
                timeResArgs['timeResSigmaSF'] = dict( Name = 'timeResSigmaSF', Value = 1.245, Error = 0.0143, MinMax = ( 0.1, 5. ))
                timeResArgs['timeResSigmaSF2'] = dict( Name = 'timeResSigmaSF2', Value = -8.812, Error = 1.507, MinMax = (-11, -1))
                covariance = {('timeResSigmaOffset', 'timeResSigmaOffset'): 2.178e-08,
                                ('timeResSigmaOffset', 'timeResSigmaSF'): 4.389e-07,
                                ('timeResSigmaOffset', 'timeResSigmaSF2'): -0.000141,
                                ('timeResSigmaSF', 'timeResSigmaSF'): 0.0002041,
                                ('timeResSigmaSF', 'timeResSigmaSF2'): 0.001858,
                                ('timeResSigmaSF2', 'timeResSigmaSF2'): 2.271}
                timeResArgs['Covariance'] = covariance
                timeResArgs['timeResSFModel'] = 'quadratic'
            else :
                timeResArgs['timeResMean'] = dict( Value = -4.0735e-03, Error = 1.33e-04 )
                timeResArgs['timeResMeanConstraint'] = constrTResScale
                timeResArgs['timeResSFConstraint'] = constrTResScale

            self._timeResModel = TimeResolution( **timeResArgs )

        elif timeResType == '3Gauss' :
            from P2VV.Parameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution( time = observables['time'], timeResSFConstraint = constrTResScale)

        else :
            from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
            self._timeResModel = TimeResolution(  time          = observables['time']
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
                         , [ sigData, cbkgData ]
                         , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                        ) :
                plot(  pad, observables['timeRes'], data, None
                     , frameOpts  = dict( Title = observables['timeRes'].GetTitle() + plotTitle )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 , MarkerColor = kBlue, LineColor = kBlue )
                    )


        # CP violation parameters
        if lambdaCPParam == 'ReIm' : 
            from P2VV.Parameterizations.CPVParams import LambdaCarth_CPParam as CPParam
            ReLambdaCPVar = dict( Name = 'ReLambdaCP' )
            ImLambdaCPVar = dict( Name = 'ImLambdaCP' )
            self._lambdaCP = CPParam( ReLambdaCP = ReLambdaCPVar, ImLambdaCP = ImLambdaCPVar )

        elif lambdaCPParam == 'lambPhiRel_CPVDecay' :
            from P2VV.Parameterizations.CPVParams import LambdaAbsArgRel_CPVDecay_CPParam as CPParam
            phiCPVar = dict( Name = 'phiCP_m' )
            self._lambdaCP = CPParam( phiCP_m = phiCPVar, AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes )
            if ambiguityPars : self._lambdaCP['phiCP_m'].setVal( pi - self._lambdaCP['phiCP_m'].getVal() )

        elif lambdaCPParam.startswith('lambPhi_CPVDecay') :
            from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPVDecay_CPParam as CPParam
            if lambdaCPParam.endswith('PSWaves') :
                rhoCP_P = RealVar( 'rhoCP_P', Title = 'CPV in decay param. |rho|', Value = 1., Error = 0.04, MinMax = ( 0., 5. ) )
                rhoCP_S = RealVar( 'rhoCP_S', Title = 'CPV in decay param. |rho|', Value = 1., Error = 0.04, MinMax = ( 0., 5. ) )
                phiCP_P = RealVar( 'phiCP_P', Title = 'CPV in decay param. phi',   Value = 0., Error = 0.1,  MinMax = (-RooInf, +RooInf) )
                phiCP_S = RealVar( 'phiCP_S', Title = 'CPV in decay param. phi',   Value = 0., Error = 0.1,  MinMax = (-RooInf, +RooInf) )
                self._lambdaCP = CPParam(  AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes
                                         , rhoCP_A0 = rhoCP_P, rhoCP_Apar = rhoCP_P, rhoCP_Aperp = rhoCP_P, rhoCP_AS = rhoCP_S
                                         , phiCP_A0 = phiCP_P, phiCP_Apar = phiCP_P, phiCP_Aperp = phiCP_P, phiCP_AS = phiCP_S
                                        )

            else :
                self._lambdaCP = CPParam( AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self._amplitudes )

        else :
            if lambdaCPParam == 'lambPhi' :
                from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
            else :
                from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam

            phiCPVar = dict( Name = 'phiCP' )
            lambdaCPVar = dict( Name = 'lambdaCP' )
            if blind and 'phiCP' in blind :
                if blind['phiCP'] : phiCPVar['Blind'] = blind['phiCP']
                else              : phiCPVar['Blind'] = ( 'UnblindUniform', 'BsPhis2013EPS',  0.2 )
            if blind and 'lambdaCP' in blind :
                if blind['lambdaCP'] : phiCPVar['Blind'] = blind['lambdaCP']
                else                 : phiCPVar['Blind'] = ( 'UnblindUniform', 'BsLambdas2013EPS', 0.1 )
            self._lambdaCP = CPParam( phiCP = phiCPVar, lambdaCP = lambdaCPVar )
            if ambiguityPars :
                self._lambdaCP['phiCP'].setVal( pi - self._lambdaCP['phiCP'].getVal() )

        # coefficients for time functions
        from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, self._amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] )

        # tagging parameters
        if sigTaggingPdf == 'histPdf' :
            raise RuntimeError('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: A histogram tagging PDF can only be used when tagging observables are conditional')

        else :
            # get tagging category parameters dictionary/dictionaries
            tagCatsDictOS = self._tagCatsOS.tagCatsDict()
            if SSTagging : tagCatsDictSS = self._tagCatsSS.tagCatsDict()

            if sigTaggingPdf.startswith('tagUntag') :
                # assume products of asymmetries are small and B-Bbar asymmetries are equal for all tagged categories
                if observables['tagCatOS'].numTypes() > 2 or ( SSTagging and observables['tagCatSS'].numTypes() > 2 ) :
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: tagging in signal PDF:\n'\
                        + '    * assuming B-Bbar asymmetries are equal for all tagged categories'

                # provide the same asymmetry for all tagged categories
                asymVal = 0. #asymVal = -self._lambdaCP['C'].getVal()
                asymErr = ( 10. / sqrt( dataSet.sumEntries() ) ) if dataSet else 0.1
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
                                       , Conditionals = tagCatsDictOS['Conditionals'] | tagCatsDictSS['Conditionals']
                                       , Constraints  = tagCatsDictOS['Constraints']  | tagCatsDictSS['Constraints']
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

            from P2VV.Parameterizations.FlavourTagging import CatDilutionsCoefAsyms_TaggingParams as TaggingParams
            self._taggingParams = TaggingParams( **tagCatsDict )

            if condTagging :
                # don't float tagging category coefficients if PDF is conditional on tagging observables
                for coefList in self._taggingParams['singleTagCatCoefs'] :
                    for coef in coefList :
                        if coef.isFundamental() : coef.setConstant(True)

            if not SSTagging :
                args = dict(  tagCat      = observables['tagCatOS']
                            , iTag        = observables['iTagOS']
                            , dilutions   = self._taggingParams['dilutions']
                            , ADilWTags   = self._taggingParams['ADilWTags']
                            , avgCEvens   = self._taggingParams['avgCEvens']
                            , avgCOdds    = self._taggingParams['avgCOdds']
                            , tagCatCoefs = self._taggingParams['tagCatCoefs']
                           )
            else :
                args = dict(  tagCat0     = observables['tagCatOS']
                            , tagCat1     = observables['tagCatSS']
                            , iTag0       = observables['iTagOS']
                            , iTag1       = observables['iTagSS']
                            , dilutions0  = self._taggingParams['dilutions'][0]
                            , dilutions1  = self._taggingParams['dilutions'][1]
                            , ADilWTags0  = self._taggingParams['ADilWTags'][0]
                            , ADilWTags1  = self._taggingParams['ADilWTags'][1]
                            , avgCEvens   = self._taggingParams['avgCEvens']
                            , avgCOdds    = self._taggingParams['avgCOdds']
                            , tagCatCoefs = self._taggingParams['tagCatCoefs']
                           )

        args.update(  time                   = observables['time']
                    , tau                    = self._lifetimeParams['MeanLifetime']
                    , dGamma                 = self._lifetimeParams['dGamma']
                    , dm                     = self._lifetimeParams['dM']
                    , coshCoef               = timeBasisCoefs['cosh']
                    , sinhCoef               = timeBasisCoefs['sinh']
                    , cosCoef                = timeBasisCoefs['cos']
                    , sinCoef                = timeBasisCoefs['sin']
                    , resolutionModel        = self._timeResModel['model']
                    , ConditionalObservables = self._amplitudes.conditionalObservables()\
                                               | self._timeResModel.conditionalObservables()\
                                               | self._taggingParams.conditionalObservables()
                    , ExternalConstraints    = self._lifetimeParams.externalConstraints()\
                                               | self._timeResModel.externalConstraints()\
                                               | self._taggingParams.externalConstraints()
                   )

        # build signal PDF
        from P2VV.RooFitWrappers import BTagDecay
        if condTagging : sigPdf = BTagDecay( 'sig_t_angles',             **args )
        else :           sigPdf = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )
        self._BTagDecay = sigPdf


        ###################################################################################################################################
        ## time acceptance ##
        #####################

        if multiplyByTimeEff in [ 'all', 'signal', 'background' ] :
            hlt1ExclB = observables['hlt1ExclB']
            hlt2B     = observables['hlt2B']
            hlt2UB    = observables['hlt2UB']
            self._timeResModelOriginal = self._timeResModel
            if timeEffType == 'fit' and selection == 'timeEffFit' :
                hists = {  hlt1ExclB : {  'exclB'    : { 'histogram' : 'hlt1_shape', 'average' : ( 6.285e-01, 1.633e-02 ) }
                                        , 'notExclB' : { 'bins'      : observables['time'].getRange(), 'heights' : [0.5]  }
                                       }
                         , hlt2B     : { 'B'         : { 'histogram' : 'hlt2_shape', 'average' : ( 6.3290e-01, 1.65e-02 ) } }
                         , hlt2UB    : { 'UB'        : { 'bins'      : observables['time'].getRange(), 'heights' : [0.5]  } }
                        }

                from P2VV.Parameterizations.TimeAcceptance import Paper2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance( time = observables['time'], Input = timeEffHistFile, Histograms = hists
                                                    , Data = dataSet, Fit = True, Original = sigPdf
                                                    , ResolutionModel = self._timeResModel
                                                    , **timeEffParameters )

            elif timeEffType == 'paper2012' and selection == 'paper2012' :
                hists = { hlt1ExclB : {  'exclB'    : { 'histogram' : timeEffHistExclBName }
                                       , 'notExclB' : { 'histogram' : timeEffHistUBName    }
                                      }
                        }

                from P2VV.Parameterizations.TimeAcceptance import Paper2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance( time = observables['time'], Input = timeEffHistFile, Histograms = hists
                                                    , Data = dataSet, Fit = False, Original = sigPdf
                                                    , ResolutionModel = self._timeResModel, BinHeightMinMax = ( -RooInf, RooInf )
                                                    , **timeEffParameters )

            elif timeEffType in [ 'HLT1Unbiased', 'HLT1ExclBiased' ] or ( timeEffType == 'paper2012' and selection == 'paper2012' ) :
                from P2VV.Parameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
                self._timeResModel = TimeAcceptance(  time = observables['time']
                                                    , Input = timeEffHistFile
                                                    , Histogram = timeEffHistExclBName if timeEffType == 'HLT1ExclBiased'\
                                                                  else timeEffHistUBName
                                                    , ResolutionModel = self._timeResModel
                                                    , **timeEffParameters )
            else:
                raise ValueError( 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: unknown time efficiency type: "%s" (with "%s" selection)'\
                                 % ( timeEffType, selection ) )

        if multiplyByTimeEff in [ 'all', 'signal' ] :
            # multiply signal PDF with time acceptance
            print 'P2VV - INFO:  Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with lifetime efficiency function'
            args.update( resolutionModel= self._timeResModel['model']
                       , ConditionalObservables =  args['ConditionalObservables'] | self._timeResModel.conditionalObservables() 
                       , ExternalConstraints =  args['ExternalConstraints'] | self._timeResModel.externalConstraints()  
                       )
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

            from P2VV.GeneralUtils import RealMomentsBuilder,angularMomentIndices
            moments = RealMomentsBuilder()
            moments.appendPYList( self._angleFuncs.angles, angularMomentIndices(multiplyByAngEff,self._angleFuncs ) )
            moments.read(angEffMomentsFile)
            moments.Print()
            sigPdfTimeAcc = moments * sigPdfTimeAcc

        self._BTagDecayAngEff = sigPdfTimeAcc

        # add PDF to signal components
        self._signalComps += sigPdfTimeAcc
        if condTagging : self._sig_t_angles = sigPdfTimeAcc
        else           : self._sig_t_angles_tagCat_iTag = sigPdfTimeAcc


        ###################################################################################################################################
        ## plot mumu mass ##
        ####################

        if makePlots :
            self._mumuMassCanv = TCanvas( 'mumuMassCanv', 'mumu Mass' )
            for ( pad, data, plotTitle )\
                  in zip(  self._mumuMassCanv.pads( 2, 2 )
                         , [ sigData, cbkgData ]
                         , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                        ) :
                plot(  pad, observables['mumuMass'], data, None
                     , frameOpts  = dict( Title = observables['mumuMass'].GetTitle() + plotTitle )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 , MarkerColor = kBlue, LineColor = kBlue )
                    )


        ###################################################################################################################################
        ## build KK mass PDFs ##
        ########################

        if makePlots :
            # build the signal and background KK mass PDFs
            from P2VV.Parameterizations.MassPDFs import Binned_MassPdf
            self._signalKKMass = Binned_MassPdf( 'sig_mKK', observables['KKMass'], Binning = self._KKMassBinning
                                                , Data = sigData )
            self._backgroundKKMass = Binned_MassPdf( 'bkg_mKK', observables['KKMass'], Binning = self._KKMassBinning
                                                    , Data = cbkgData )

            self._KKMassCanv = TCanvas( 'KKMassCanv', 'KK Mass' )
            for ( pad, data, pdf, plotTitle, scale )\
                  in zip(  self._KKMassCanv.pads( 2, 2 )
                         , [ sigData, cbkgData ]
                         , [ self._signalKKMass.pdf(), self._backgroundKKMass.pdf() ]
                         , [ ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ ( 5., 1.e4 ), ( 1.2e2, 1.2e3 ) ]
                        ) :
                plot(  pad, observables['KKMass'], data, pdf, logy = True, yScale = scale
                     , frameOpts  = dict( Title = observables['KKMass'].GetTitle() + plotTitle )
                     , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )#, MarkerColor = kBlue, LineColor = kBlue )
                    )


        ###################################################################################################################################
        ## build signal tagging PDF ##
        ##############################

        if condTagging and multiplyByTagPdf :
            if sigData and sigTaggingPdf == 'histPdf' :
                # create histogram from signal data and use the (fixed) bin coefficients for the PDF
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: creating signal tagging PDFs from signal data'
                from P2VV.RooFitWrappers import HistPdf
                if not SSTagging :
                    self._sig_tagCat_iTag = HistPdf(  Name = 'sig_tagCat_iTag'
                                                    , Observables = [ observables['tagCatOS'], observables['iTagOS'] ]
                                                    , Data = sigData
                                                   )
                else :
                    self._sig_tagCat_iTag = HistPdf(  Name = 'sig_tagCat_iTag'
                                                    , Observables = [ observables['tagCatOS'], observables['tagCatSS']
                                                                     , observables['iTagOS'], observables['iTagSS'] ]
                                                    , Data = sigData
                                                   )
                self._signalComps += self._sig_tagCat_iTag

            else :
                # use a PDF with variable bin coefficients
                if sigTaggingPdf.startswith('tagUntag') or ( observables['tagCatOS'].numTypes() == 2\
                                                            and ( not SSTagging or observables['tagCatSS'].numTypes() == 2 ) ) :
                    # assume B-Bbar asymmetry is equal for all tagged categories
                    if observables['tagCatOS'].numTypes() > 2 or ( SSTagging and observables['tagCatSS'].numTypes() > 2 ) :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: signal tagging PDF:\n'\
                            + '    * assuming B-Bbar asymmetries are equal for all tagged categories'
                    else :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a binned signal tagging PDF'

                    from P2VV.Parameterizations.FlavourTagging import TagUntag_BinnedTaggingPdf as TaggingPdf

                else :
                    # create independent tagging bin coefficients
                    if sigData :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining signal tagging coefficients from signal data'
                    else :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: WARNING:'\
                            + ' no signal data available to determine signal tagging coefficients:\n'\
                            + '    * assuming absence of B-Bbar asymmetries'

                    from P2VV.Parameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

                # build PDF
                self._sigTaggingPdf = TaggingPdf(  'tagCat_iTag'
                                                 , observables['tagCatOS'], observables['tagCatSS'] if SSTagging else None
                                                 , observables['iTagOS'], observables['iTagSS'] if SSTagging else None
                                                 , NamePrefix    = 'sig'
                                                 , TagCatCoefs   = [ sigPdf.tagCatCoefs(it)\
                                                                     for it in range( observables['tagCatOS'].numTypes() ) ]\
                                                                     if SSTagging else [ sigPdf.tagCatCoefs(0) ]
                                                 , TaggedCatName = 'TagCat' if observables['tagCatOS'].numTypes() > 2 else 'Tagged'
                                                 , Data          = sigData
                                                 , RelativeCoefs = True if sigTaggingPdf.endswith('Relative') else False
                                                 , Dilutions     = [ sigPdf.dilutions(False), sigPdf.dilutions(True) ]\
                                                                   if sigTaggingPdf.startswith('tagUntag') and SSTagging else None
                                                )

                self._sig_tagCat_iTag = self._sigTaggingPdf.pdf()
                self._signalComps += self._sig_tagCat_iTag


        ###################################################################################################################################
        ## build PDFs for estimated wrong-tag probabilities ##
        ######################################################

        if makePlots and fullData and sigData and cbkgData :
            # build PDFs for estimated wrong-tag probabilities
            from P2VV.RooFitWrappers import HistPdf
            self._estWTagDataOS = sigData.reduce( '%s > 0' % observables['tagCatOS'].GetName() )
            self._estWTagDataSS = sigData.reduce( '%s > 0' % observables['tagCatSS'].GetName() )

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
            untagFracOS    = fullData.table(observables['tagCatOS']).getFrac('Untagged')
            untagFracSigOS = sigData.table(observables['tagCatOS']).getFrac('Untagged')
            untagFracBkgOS = cbkgData.table(observables['tagCatOS']).getFrac('Untagged')
            untagFracSS    = fullData.table(observables['tagCatSS']).getFrac('Untagged')
            untagFracSigSS = sigData.table(observables['tagCatSS']).getFrac('Untagged')
            untagFracBkgSS = cbkgData.table(observables['tagCatSS']).getFrac('Untagged')

            # plot estimated wrong-tag probabilities for signal and for background
            self._estWTagCanvOS = TCanvas( 'estWTagCanvOS', 'Estimated wrong-tag probability OS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvOS.pads( 1, 1 ) if SFit else self._estWTagCanvOS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ sigData if SFit else fullData
                            , sigData, cbkgData ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigOS if SFit else untagFracOS ), 1. - untagFracSigOS, 1. - untagFracBkgOS ]
                        ) :
                pad.SetLeftMargin(0.18)
                pad.SetRightMargin(0.05)
                pad.SetBottomMargin(0.18)
                pad.SetTopMargin(0.05)

                #plot(  pad, estWTagOS, data, self._sig_bkg_estWTagOS
                plot(  pad, estWTagOS,  data, yScale = ( 0., None ), xTitleOffset = 1.10, yTitleOffset = 1.15
                     , xTitle     = '#eta^{OS}'
                     , yTitle     = 'Candidates / %.2f' % ( 0.499999 / float(nBins) )
                     , frameOpts  = dict( Bins = nBins, Title = estWTagOS.GetTitle() + plotTitle, Range = ( 0., 0.499999 )
                                         , Name = estWTagOS.GetName() )
                     , dataOpts   = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, LineWidth = 3 )
                     #, pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
                    )

            self._estWTagCanvSS = TCanvas( 'estWTagCanvSS', 'Estimated wrong-tag probability SS' )
            for ( pad, data, nBins, plotTitle, norm )\
                  in zip(  self._estWTagCanvSS.pads( 1, 1 ) if SFit else self._estWTagCanvSS.pads( 2, 2, lambda pad : pad != 2 )
                         , [ sigData if SFit else fullData, sigData, cbkgData ]
                         , 3 * [ numEstWTagBins ]
                         , [ '', ' - signal (B mass S-weights)', ' - background (B mass S-weights)' ]
                         , [ 1. - ( untagFracSigSS if SFit else untagFracSS ), 1. - untagFracSigSS, 1. - untagFracBkgSS ]
                        ) :
                pad.SetLeftMargin(0.18)
                pad.SetRightMargin(0.05)
                pad.SetBottomMargin(0.18)
                pad.SetTopMargin(0.05)

                #plot(  pad, estWTagSS, data, self._sig_bkg_estWTagSS
                plot(  pad, estWTagSS, data, yScale = ( 0., None ), xTitleOffset = 1.10, yTitleOffset = 1.15
                     , xTitle     = '#eta^{SSK}'
                     , yTitle     = 'Candidates / %.2f' % ( 0.499999 / float(nBins) )
                     , frameOpts  = dict( Bins = nBins, Title = estWTagSS.GetTitle() + plotTitle, Range = ( 0., 0.499999 )
                                         , Name = estWTagSS.GetName() )
                     , dataOpts   = dict( MarkerStyle = kFullCircle, MarkerSize = 0.7, LineWidth = 3 )
                     #, pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
                    )

            #TODO: Make it availble in the case of cFit
            if AllTagPlots:
                self._estWTagCanvSS_B    = TCanvas('estWTagCanvSS_B'   , 'Est. wrong-tag probability SS, B'   )
                self._estWTagCanvSS_Bbar = TCanvas( 'estWTagCanvSS_Bbar','Est. wrong-tag probability SS, Bbar')
                self._estWTagCanvOS_B    = TCanvas( 'estWTagCanvOS_B'   ,'Est. wrong-tag probability OS, B'   )
                self._estWTagCanvOS_Bbar = TCanvas( 'estWTagCanv0S_Bbar','Est. wrong-tag probability OS, Bbar')

                self._estWTagCanvSSOS_B     = TCanvas( 'estWTagCanvSSOS_B'    , 'Est. wrong-tag probability SS+OS, B'    )
                self._estWTagCanvSSOS_Bbar  = TCanvas( 'estWTagCanvSS0S_Bbar' , 'Est. wrong-tag probability SS+OS, Bbar' )
                self._estWTagCanvSSOS_BbarB = TCanvas( 'estWTagCanvSS0S_BbarB', 'Est. wrong-tag probability SS+OS, BbarB')

                tagCutSS_B    = 'tagdecision_ss==tagdecision_ss::B'
                tagCutSS_Bbar = 'tagdecision_ss==tagdecision_ss::Bbar'
                tagCutOS_B    = 'tagdecision_os==tagdecision_os::B'
                tagCutOS_Bbar = 'tagdecision_os==tagdecision_os::Bbar'

                tagCutComb_B     = 'tagdecision==tagdecision::B'
                tagCutComb_Bbar  = 'tagdecision==tagdecision::Bbar'
                tagCutComb_BbarB = tagCutComb_B + '|' + tagCutComb_Bbar

                for ( pad, data, nBins, BorBbar, titleX, obs )\
                        in zip(  [ self._estWTagCanvSS_B  ,   self._estWTagCanvSS_Bbar,
                                   self._estWTagCanvOS_B  ,   self._estWTagCanvOS_Bbar,
                                   self._estWTagCanvSSOS_B,   self._estWTagCanvSSOS_Bbar,
                                   self._estWTagCanvSSOS_BbarB ]

                                 , [ sigData.reduce(tagCutSS_B      ),
                                     sigData.reduce(tagCutSS_Bbar   ),
                                     sigData.reduce(tagCutOS_B      ),
                                     sigData.reduce(tagCutOS_Bbar   ),
                                     sigData.reduce(tagCutComb_B    ),
                                     sigData.reduce(tagCutComb_Bbar ),
                                     sigData.reduce(tagCutComb_BbarB) ]

                                 , 7 * [ numEstWTagBins ]
                                 , 3 * [ 'B', 'Bbar'  ] +     [ 'BbarB'    ]
                                 , 2 * [ '#eta^{SS}'  ] + 2 * [ '#eta^{OS}'] + 3 * [ '#eta^{SS+OS}']
                                 , 2 * [  estWTagSS   ] + 2 * [  estWTagOS ] + 3 * [  estWTagComb  ]
                                 ) :
                    plot(  pad, obs , data
                         , xTitle    = titleX
                         , yScale    = [0, None]
                         , frameOpts = dict( Bins = nBins, Title = obs.GetTitle() + BorBbar, Range = ( 0., 0.499999 )
                                            , Name = obs.GetName() + BorBbar  )
                         , dataOpts  = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts   = dict( LineColor = kBlue, LineWidth = 3  )
                        )


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
        self._backgroundTime = BackgroundTime(  Name = 'bkg_t', time = observables['time'], resolutionModel = self._timeResModel['model']
                                              , Efficiency = timeAcceptance if multiplyByTimeEff in [ 'all', 'background' ] else None
                                             )
        self._bkg_t = self._backgroundTime.pdf()

        if not SFit : self._backgroundComps += self._bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit or makePlots :
            if cbkgData and bkgAnglePdf == 'histPdf' :
                # create a histogram from background data and use it for background angular PDF
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining angular shape of background from background data'\

                from P2VV.RooFitWrappers import HistPdf
                nBins = numAngleBins
                self._bkg_angles = HistPdf(  Name = 'bkg_angles'
                                           , Observables = [ observables['cpsi'], observables['ctheta'], observables['phi'] ]
                                           , Binning = { cpsi : nBins[0], ctheta : nBins[1], phi : nBins[2] }
                                           , Data = cbkgData
                                          )

            else :
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a %s angular PDF for background'\
                      % ( 'hybrid' if bkgAnglePdf == 'hybrid' else 'basis functions' if bkgAnglePdf == 'basis' else 'binned' )

                if bkgAnglePdf != 'basis' :
                    # create a binned PDF for background angular shape
                    from array import array
                    baselineBin = -1
                    if transAngles :
                        nBins = [ 5, 7, 9 ]
                        cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
                        cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
                        phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )

                    elif bkgAnglePdf == 'hybrid' :
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
                                                  , Constant = False if bkgAnglePdf == 'hybrid' and bin1 != baselineBin else True
                                                 )\
                                          if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                                          for bin2 in range( phiNumBins ) for bin1 in range( cthNumBins ) for bin0 in range( cpsNumBins )
                                        ]

                    from P2VV.RooFitWrappers import ComplementCoef
                    self._bkgAngCoefs[0] = ComplementCoef(  Name         = 'bkg_angBin_0_0_0'
                                                          , Coefficients = self._bkgAngCoefs[ 1 : ]
                                                         )

                    # create a BinnedPdf
                    from P2VV.RooFitWrappers import BinnedPdf
                    self._bkg_angBins = BinnedPdf(  Name = 'bkg_angBins' if bkgAnglePdf == 'hybrid' else 'bkg_angles'
                                                  , Observables = [ observables['cpsi'], observables['ctheta'], observables['phi'] ]
                                                  , Binnings = [ cpsBins, cthBins, phiBins ]
                                                  , Coefficients = self._bkgAngCoefs
                                                  , BinIntegralCoefs = True
                                                 )
                    self._bkg_angBins.setForceUnitIntegral(True)

                    sumWeights = 0.
                    sumBinWeights = ( cpsNumBins * cthNumBins * phiNumBins - 1 ) * [ 0. ]
                    spikesFrac = 1.
                    if cbkgData :
                        # determine bin coefficient values
                        angleInitVals = [ observables[name].getVal() for name in [ 'cpsi', 'ctheta', 'phi' ] ]
                        for obsSet in cbkgData :
                            for name in [ 'cpsi', 'ctheta', 'phi' ] :
                                observables[name].setVal( obsSet.getRealValue( observables[name].GetName() ) )
                            sumWeights += cbkgData.weight()
                            cpsBin = cpsi.getBin('bkg_cpsBins')
                            cthBin = ctheta.getBin('bkg_cthBins')
                            phiBin = phi.getBin('bkg_phiBins')
                            bin = cpsBin + cthBin * cpsNumBins + phiBin * cpsNumBins * cthNumBins - 1
                            if bin >= 0 : sumBinWeights[ bin ] += cbkgData.weight()
                        for name, val in zip( [ 'cpsi', 'ctheta', 'phi' ], angleInitVals ) : observables[name].setVal(val)

                        if bkgAnglePdf == 'hybrid' :
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
                            if bkgAnglePdf == 'hybrid' :
                                coef.setVal( ( ( binInt - baseBinVal * cthBins.binWidth(bin + 1) ) / spikesFrac )\
                                             if bin != baselineBin - 1 else 0. )
                                #coef.setVal( binCoefVals[bin] )
                                #coef.setError( binCoefErrs[bin] )
                            else :
                                coef.setVal(binInt)

                    if bkgAnglePdf == 'hybrid' :
                        self._bkgAngBinsFrac = RealVar(  'bkgAngBinsFrac'
                                                       , Title  = 'Binned PDF fraction in angular background'
                                                       , Value  = spikesFrac
                                                       , MinMax = ( 0., 1. )
                                                      )

                if bkgAnglePdf in [ 'basis', 'hybrid' ] :
                    if bkgAnglePdf == 'basis' : nBins = [ 20, 20, 20 ]
                    # create an angular PDF with basis functions
                    angPDFIndices = [ ( 2, 0, 0 ), ( 0, 2, 0 ) ]
                    #angPDFIndices = [ ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 0, 0 ), ( 2, 2, 0 ), ( 2, 2, 2 ) ]
                    #angPDFIndices = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
                    #                  for YIndex1 in range( -YIndex0, YIndex0 + 1 )\
                    #                  if PIndex != 0 or YIndex0 != 0 or YIndex1 != 0 ]

                    from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms as angularPdfTerms
                    cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
                    angPdfTerms = angularPdfTerms(  Angles = self._angleFuncs.angles
                                                  , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                             , {  'Name'    : 'bkg_angCoef_%s_%s_%s'\
                                                                              % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                                , 'Value'   : 0.
                                                                , 'Error'   : 0.01
                                                                , 'MinMax'  : ( -0.4, 0.4 )
                                                                , 'Indices' : inds
                                                               }
                                                            ) for inds in angPDFIndices
                                                          )
                                                 )
                    self._bkg_angFuncs = angPdfTerms.buildSumPdf('bkg_angFuncs')

                # build total angular PDF
                if bkgAnglePdf == 'hybrid' :
                    from P2VV.RooFitWrappers import SumPdf
                    self._bkg_angles = SumPdf(  Name   = 'bkg_angles'
                                              , PDFs   = [ self._bkg_angBins, self._bkg_angFuncs ]
                                              , Yields = { self._bkg_angBins.GetName() : self._bkgAngBinsFrac }
                                             )

                elif bkgAnglePdf == 'basis' :
                    self._bkg_angles = self._bkg_angFuncs

                else :
                    self._bkg_angles = self._bkg_angBins

                if cbkgData :
                    # fit background angular distribution
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: fitting background angular distribution'
                    print '  initial parameter values:'
                    for par in self._bkg_angles.getParameters(cbkgData) : par.Print()
                    self._bkg_angles.fitTo( cbkgData, SumW2Error = False, Save = False, **fitOpts )

            if not SFit : self._backgroundComps += self._bkg_angles

            if cbkgData and cbkgData and makePlots :
                # plot background angles with S-weights
                self._bkgAnglesSWeightCanv = TCanvas( 'bkgAnglesSWeightCanv', 'Background Decay Angles with S-Weights' )
                for ( pad, obs, data, bins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesSWeightCanv.pads( 3, 2 )
                             , [ observables[name] for name in [ 'cpsi', 'ctheta', 'phi' ] ]
                             , 3 * ( cbkgData, )
                             , nBins
                             , [ observables[name].GetTitle() for name in [ 'cpsi', 'ctheta', 'phi' ] ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle, yTitleOffset = 1.5
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
                        )

                # plot background angles with side bands
                self._bkgAnglesSideBandCanv = TCanvas( 'bkgAnglesSideBandCanv', 'Background Decay Angles with Side Bands' )
                for ( pad, obs, data, bins, plotTitle, xTitle, scale )\
                      in zip(  self._bkgAnglesSideBandCanv.pads( 3, 2 )
                             , [ observables[name] for name in [ 'cpsi', 'ctheta', 'phi' ] ]
                             , 3 * ( cbkgData, )
                             , nBins
                             , [ observables[name].GetTitle() for name in [ 'cpsi', 'ctheta', 'phi' ] ]
                             , ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                             , ( ( None, None ), ( None, None ), ( None, None ) ) # ( ( 0., 740. * 9585. / 12005. ), ( 0., 580. * 9585. / 12005. ), ( 0., 760. * 9585. / 12005. ) ) # ( ( 0., 740. ), ( 0., 580. ), ( 0., 760. ) )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle, yTitleOffset = 1.5, yScale = scale
                         , frameOpts  = dict( Bins = bins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3  )
                        )

                # plot 2-D angular distributions
                #for angle0, angle1, data, canv, padNr in [  ( 1, 0, cbkgData, self._bkgAnglesSWeightCanv,  4 )
                #                                          , ( 2, 0, cbkgData, self._bkgAnglesSWeightCanv,  5 )
                #                                          , ( 1, 2, cbkgData, self._bkgAnglesSWeightCanv,  6 )
                #                                          , ( 1, 0, cbkgData,   self._bkgAnglesSideBandCanv, 4 )
                #                                          , ( 2, 0, cbkgData,   self._bkgAnglesSideBandCanv, 5 )
                #                                          , ( 1, 2, cbkgData,   self._bkgAnglesSideBandCanv, 6 )
                #                                         ] :
                #    bkgAngHist = data.createHistogram( angles[angle0]._var, angles[angle1]._var, nBins[angle0], nBins[angle1] )
                #    _P2VVPlotStash.append(bkgAngHist)
                #    bkgAngHist.SetStats(False)
                #    bkgAngHist.SetTitle( '%s vs. %s' % ( angleNames[angle0][1], angleNames[angle1][1] ) )
                #    bkgAngHist.SetMinimum(0.)
                #    bkgAngHist.GetXaxis().SetTitle( angleNames[angle0][1] )
                #    bkgAngHist.GetYaxis().SetTitle( angleNames[angle1][1] )
                #    bkgAngHist.GetXaxis().SetLabelOffset(0.01)
                #    bkgAngHist.GetYaxis().SetLabelOffset(0.008)
                #    bkgAngHist.GetXaxis().SetTitleOffset(1.8)
                #    bkgAngHist.GetYaxis().SetTitleOffset(1.8)
                #    bkgAngPad = canv.cd(padNr)
                #    bkgAngPad.SetLeftMargin(0.08)
                #    bkgAngPad.SetRightMargin(0.05)
                #    bkgAngPad.SetBottomMargin(0.05)
                #    bkgAngPad.SetTopMargin(0.05)
                #    bkgAngHist.Draw('lego2')


        ###################################################################################################################################
        ## build background tagging PDF ##
        ##################################

        if not SFit and ( multiplyByTagPdf or not condTagging ) :
            # tagCat = { Untagged, TagCat1, TagCat2, ... }, tag = { B, Bbar }
            if cbkgData and bkgTaggingPdf == 'histPdf' :
                # create histogram from background data and use the (fixed) bin coefficients for the PDF
                print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: creating background tagging PDFs from background data'
                from P2VV.RooFitWrappers import HistPdf
                self._bkg_tagCat_iTag = HistPdf(  Name = 'bkg_tagCat_iTag'
                                                , Observables = [ observables['tagCatOS'], observables['iTagOS'] ]
                                                , Data = cbkgData
                                               )
                self._backgroundComps += self._bkg_tagCat_iTag

            else :
                # use a PDF with variable bin coefficients
                if bkgTaggingPdf.startswith('tagUntag') or ( observables['tagCatOS'].numTypes() == 2\
                                                            and ( not SSTagging or observables['tagCatSS'].numTypes() == 2 ) ) :
                    # couple background tagging category coefficients to signal tagging category coefficients
                    # and assume B-Bbar asymmetry is equal for all tagged categories
                    if observables['tagCatOS'].numTypes() > 2  or ( SSTagging and observables['tagCatSS'].numTypes() > 2 ) :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: background tagging PDF:\n'\
                            + '    * assuming B-Bbar asymmetries are equal for all tagged categories'
                    else :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building a binned background tagging PDF'

                    from P2VV.Parameterizations.FlavourTagging import TagUntag_BinnedTaggingPdf as TaggingPdf

                else :
                    # create independent tagging bin coefficients
                    if cbkgData :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: determining background tagging coefficients from background data'
                    else :
                        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: WARNING:'\
                            + ' no background data available to determine background tagging coefficients:\n'\
                            + '    * assuming equal tagging category coefficients\n'\
                            + '    * assuming absence of B-Bbar asymmetries'

                    from P2VV.Parameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

                # build PDF
                self._bkgTaggingPdf = TaggingPdf(  'tagCat_iTag'
                                                 , observables['tagCatOS'], observables['tagCatSS'] if SSTagging else None
                                                 , observables['iTagOS'], observables['iTagSS'] if SSTagging else None
                                                 , NamePrefix    = 'bkg'
                                                 , TaggedCatName = 'TagCat' if observables['tagCatOS'].numTypes() > 2 else 'Tagged'
                                                 , Data          = cbkgData
                                                 , RelativeCoefs = True if bkgTaggingPdf.endswith('Relative') else False
                                                )

                self._bkg_tagCat_iTag = self._bkgTaggingPdf.pdf()
                self._backgroundComps += self._bkg_tagCat_iTag


        ###################################################################################################################################
        ## build full PDF ##
        ####################

        from P2VV.RooFitWrappers import buildPdf
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

        if ( not SFit and selection in ['paper2012', 'timeEffFit'] ) or paramKKMass == 'simultaneous' :
            # categories for splitting the PDF
            splitCats = [ [ ] ]
            if not SFit :
                splitCats[0] += [ hlt1ExclB ] if selection == 'paper2012' else [ hlt1ExclB, hlt2B ] if selection == 'timeEffFit' else [ ]
            splitCats[0] += [ observables['KKMassCat'] ] if paramKKMass == 'simultaneous' else [ ]

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
                splitCats.append( [ observables['KKMassCat'] ] )
                splitParams.append( [ ] )
                for amp in self._amplitudes.parameters() :
                    if any( name in amp.GetName() for name in [ 'AS', 'A_S', 'fS', 'f_S', 'C_SP' ] ) : splitParams[-1].append(amp)
                if not SFit :
                    splitParams[-1] += [ par for par in self._backgroundBMass.parameters() if not par.isConstant() ]

            # build simultaneous PDF
            from P2VV.RooFitWrappers import SimultaneousPdf
            self._simulPdf = SimultaneousPdf(  self._fullPdf.GetName() + '_simul'
                                             , MasterPdf       = self._fullPdf
                                             , SplitCategories = splitCats
                                             , SplitParameters = splitParams
                                            )

            from P2VV.GeneralUtils import getSplitPar
            splitCatPars = self._simulPdf.getVariables()
            if not SFit :
                # set values for splitted yields
                splitCatIter = self._simulPdf.indexCat().typeIterator()
                splitCatState = splitCatIter.Next()
                massPars = self._sWeightMassPdf.getVariables()
                while splitCatState :
                    sigYield = getSplitPar( 'N_sigMass' if SFit else 'N_signal', splitCatState.GetName(), splitCatPars )
                    bkgYield = getSplitPar( 'N_bkgMass' if SFit else 'N_bkg',    splitCatState.GetName(), splitCatPars )

                    if not sigYield in massPars or not bkgYield in massPars :
                        if splitCat.isFundamental() :
                            selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
                        else :
                            splitCat.setLabel( splitCatState.GetName() )
                            selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
                        nEv    = fullData.sumEntries()
                        nEvBin = fullData.sumEntries(selStr)

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
                splitCatIter = observables['KKMassCat'].typeIterator()
                splitCatState = splitCatIter.Next()
                while splitCatState :
                    if ASParam != 'Mag2ReIm' :
                        # S-P-wave coupling factors
                        C_SP = getSplitPar( 'C_SP', splitCatState.GetName(), splitCatPars )
                        C_SP.setVal( CSPValues[ splitCatState.getVal() ] )
                        C_SP.setConstant(True)

                    if ( amplitudeParam == 'bank' and ASParam != 'ReIm' )\
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

        else :
            self._simulPdf = None

        assert not pdfConfig, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: superfluous arguments found: %s' % pdfConfig
        PdfBuilder.__init__( self, self._simulPdf if self._simulPdf else self._fullPdf, observables, { } )
        print 120 * '='

        # garbage collect
        gc.collect()
