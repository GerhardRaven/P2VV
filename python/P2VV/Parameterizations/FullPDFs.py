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
        dict.__init__( self, **kwargs )

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
          raise RuntimeError( 'P2VV - ERROR: PdfConfiguration.readParametersFromFile(): unable to open file \"%s\"' % filePath )

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

        print 'P2VV - INFO: PdfConfiguration.readParametersFromFile(): %d parameter%s read from file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )

    def writeParametersToFile( self, filePath = 'parameters', **kwargs ) :
        # get file path and name
        filePath = filePath.strip()
        fileName = filePath.split('/')[-1]

        # open file
        try :
            parFile = open( filePath, 'w' )
        except :
            raise RuntimeError( 'P2VV - ERROR: PdfConfiguration.writeParametersToFile(): unable to open file \"%s\"' % filePath )

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

        print 'P2VV - INFO: PdfConfiguration.writeParametersToFile(): %d parameter%s written to file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )


# B_s^0 -> J/psi phi analysis of 1 fb^-1 2011 data (paper) configuration
class Bs2Jpsiphi_2011Analysis( PdfConfiguration ) :
    def __init__( self ) :
        from math import pi

        # job parameters
        self['SFit']          = True            # fit only signal?
        self['blind']         = { }             # { 'phiCP' : ( 'UnblindUniform', 'myString', 0.2 ) }
        self['numEvents']     = 60000           # sum of event yields
        self['sigFrac']       = 0.5             # fraction of signal events
        self['parNamePrefix'] = ''              # prefix for parameter names
        self['runPeriods']    = [ 2011, 2012 ]  # run periods

        self['obsDict'] = dict(  runPeriod = ( 'runPeriod',            'run period',              { }                                   )
                               , mass      = ( 'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.   )
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
        self['readFromWS'] = False    # get observables from workspace?
        self['signalData'] = None     # data set of signal events

        # fit options
        self['fitOptions'] = dict( NumCPU = 1, Optimize = 1, Timer = True, Minimizer = 'Minuit2' )

        # PDF parameters
        self['numTimeBins']        = 30               # number of bins for the decay time observable
        self['numTimeResBins']     = 40               # number of bins for the decay-time resolution observable (used for caching)
        self['timeResType']        = 'eventNoMean'    # '' / 'event' / 'eventNoMean' / 'eventConstMean' / '3Gauss'
        self['constrainTResScale'] = ''      # '' / 'constrain' / 'fixed'
        self['timeEffType']        = 'paper2012'      # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'fit'
        self['constrainDeltaM']    = ''      # '' / 'constrain' / 'fixed'
        self['constrainBeta']      = 'noBeta'         # '' / 'constrain' / 'fixed' / 'noBeta'

        self['timeEffHistFiles']     = (  'data/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                        , 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                        , 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                       )
        self['timeEffParameters']    = { }

        self['transAngles']     = False        # use transversity angles?
        self['anglesEffType']   = 'weights'    # '' / 'weights' / 'basis012' / 'basis012Plus' / 'basis012Thetal' / 'basis0123' / 'basis01234' / 'basisSig3' / 'basisSig4'
        self['angEffMomsFiles'] = 'data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'
        self['angularRanges']   = dict( cpsi = [ ], ctheta = [ ], phi = [ ] )

        self['sigTaggingType']   = 'tagUntag'    # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        self['tagPdfType']       = ''
        self['SSTagging']        = True           # use same-side Kaon tagging?
        self['condTagging']      = True           # make tagging categories and B/Bbar tags conditional observables?
        self['contEstWTag']      = True           # use a continuous estimated wrong-tag probability instead of tagging categories?
        self['constrainTagging'] = ''    # '' / 'constrain' / 'fixed'
        self['tagCatsOS']        = [ ]            # [ ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.669, 0.0 ), ( 'Tagged',   1, 0.4999999, 0.392, 0.392, 0.0, 0.331, 0.0 ) ]
        self['tagCatsSS']        = [ ]            # [ ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.896, 0.0 ), ('Tagged',    1, 0.4999999, 0.359, 0.359, 0.0, 0.104, 0.0 ) ]

        self['amplitudeParam'] = 'phasesSWaveFrac'    # 'phases' / 'phasesSWaveFrac' / 'ReIm' / 'bank'
        self['ASParam']        = 'deltaPerp'          # 'delta0' / 'deltaPerp' / 'ReIm' / 'Mag2ReIm' / 'Mag2ReImPerp'
        self['AparParam']      = 'phase'              # 'phase' / 'ReIm' / 'Mag2ReIm' / 'cos' / 'real'
        self['ambiguityPars']  = False                # set parameters to values of the second minimum?

        self['paramKKMass']     = 'simultaneous'       # '' / 'parameters' / 'simultaneous'
        self['KKMassBinBounds'] = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ]    # KK-mass bin boundaries
        self['CSPValues']       = [ 0.959, 0.770, 0.824, 0.968 ]                      # S-P coupling factors for KK-mass bins

        self['lambdaCPParam'] = 'lambPhi'    # 'ReIm' / 'lambSqPhi' / 'lambPhi' / 'lambPhi_CPVDecay' / 'lambPhiRel_CPVDecay'

        self['splitParams'] = dict( KKMassCat = [ 'f_S', 'ASOddPhase' ] )

        self['externalConstr'] = dict(  dM             = (  17.768, 0.024  )
                                      , wTagP0OS       = (  0.392,  0.008  )
                                      , wTagP1OS       = (  1.000,  0.023  )
                                      , wTagDelP0OS    = (  0.0110, 0.0034 )
                                      , wTagDelP1OS    = (  0.000,  0.001  )
                                      , wTagP0SS       = (  0.350,  0.017  )
                                      , wTagP1SS       = (  1.00,   0.16   )
                                      , wTagDelP0SS    = ( -0.019,  0.005  )
                                      , wTagDelP1SS    = (  0.00,   0.01   )
                                      , timeResSigmaSF = (  1.45,   0.06   )
                                      , betaTimeEff    = ( -0.0083, 0.004  )
                                     )

        # initialize PdfConfiguration object
        PdfConfiguration.__init__( self )

    def __setitem__( self, key, val ) :
        if key == 'SFit' :
            if not val : self.addBkgParams()
            else       : self.rmBkgParams()

        return dict.__setitem__( self, key, val )

    def addBkgParams(self) :
        if not 'bkgAnglePdfType' in self : self['bkgAnglePdfType'] = 'hybrid'
        if not 'numAngleBins'    in self : self['numAngleBins']    = ( 10, 24, 5 )

    def rmBkgParams(self) :
        for key in [ 'bkgAnglePdfType', 'numAngleBins' ] : self.pop( key, None )


class SimulCatSettings(list) :
    def __init__( self, Name ) :
        self._name = Name
        self._cats = set()
        self._default = None

    def categories(self) :
        return self._cats

    def default(self) :
        assert self._default, 'P2VV - ERROR: SimulCatSettings.default(%s): no default settings defined' % self._name
        return self._default

    def addSettings( self, Categories, Labels, Settings ) :
        if not type(Categories) == str :
            self.append( ( dict( [ ( cat, [ lab for lab in labs ] ) for cat, labs in zip( Categories, Labels ) ] ), Settings ) )
            for cat, labs in zip( Categories, Labels ) :
                if len(labs) > 0 : self._cats.add(cat)

        else :
            self._default = Settings

    def getSettings( self, CatLabels ) :
        retSetts = self._default
        settCount = 0
        for settings in self :
            thisSett = True
            for cat, lab in CatLabels :
                if cat in settings[0] and len( settings[0][cat] ) > 0 and lab not in settings[0][cat] :
                    thisSett = False
                    break
            if thisSett :
                retSetts = settings[1]
                settCount += 1

        assert settCount < 2\
               , 'P2VV - ERROR: SimulCatSettings.getSettings(%s): multiple settings found that match criteria'\
                 % self._name
        assert retSetts\
               , 'P2VV - ERROR: SimulCatSettings.getSettings(%s): no settings found that match criteria and no default specified'\
                 % self._name

        return retSetts


###########################################################################################################################################
## PDF builders ##
##################

# function for parsing keyword arguments
def getKWArg( self, kwargs, name, default = 'noDefault' ) :
    if name in kwargs : return kwargs.pop(name)
    if 'kwargs' in self and name in self['kwargs'] : return self['kwargs'].pop(name)
    if name in self : return self[name]
    if default != 'noDefault' : return default
    raise KeyError( 'P2VV - ERROR: getKWArg(): argument "%s" not found' % name )


# PDF builder base class
class PdfBuilder ( dict ) :
    def __init__(self) :
        self['pdf']         = None
        self['observables'] = { }
        self['parameters']  = { }

    def pdf(self)         : return self['pdf']
    def observables(self) : return self['observables']
    def parameters(self)  : return self['parameters']

    def _createObservables( self, **kwargs ) :
        obsDict    = getKWArg( self, kwargs, 'obsDict' )
        readFromWS = getKWArg( self, kwargs, 'readFromWS' )

        # create observables or wrap P2VV wrappers around observables from data set
        from P2VV.RooFitWrappers import RooObject, RealVar, Category
        self['observables'] = { }
        for name in obsDict.keys() :
            if readFromWS :
                ws = RooObject().ws()
                if   ws.var( obsDict[name][0] ) : self['observables'][name] = RealVar(  obsDict[name][0] )
                elif ws.cat( obsDict[name][0] ) : self['observables'][name] = Category( obsDict[name][0] )
                else : raise RuntimeError( 'P2VV - ERROR: PdfBuilder._createObservables(): variable "%s" not in workspace' % name )
            else :
                if type( obsDict[name][2] ) != dict :
                    self['observables'][name] = RealVar( obsDict[name][0], Title = obsDict[name][1]
                                                        , Unit = obsDict[name][2], Value = obsDict[name][3]
                                                        , MinMax = ( obsDict[name][4], obsDict[name][5] )
                                                       )
                else :
                    self['observables'][name] = Category( obsDict[name][0], Title = obsDict[name][1], States = obsDict[name][2] )

        for obs in self['observables'].itervalues() : obs.setObservable(True)


# B_s^0 -> J/psi phi PDF builder
class Bs2Jpsiphi_PdfBuilder ( PdfBuilder ) :
    """builds the signal+background PDF for the measurement of phi_s in B_s -> J/psi(->mu^+ mu^-) phi(->K^+ K^-)
    """

    def __init__( self, **kwargs ) :
        # make sure keyword arguments are processed by called helper functions
        self['kwargs'] = kwargs

        # get some build parameters
        for par in [ 'SFit', 'KKMassBinBounds', 'obsDict', 'CSPValues', 'condTagging', 'contEstWTag', 'SSTagging', 'transAngles'
                    , 'numEvents', 'sigFrac', 'paramKKMass', 'amplitudeParam', 'ASParam', 'signalData', 'fitOptions', 'parNamePrefix'
                    , 'tagPdfType', 'timeEffType', 'timeEffHistFiles', 'timeEffParameters', 'anglesEffType', 'angEffMomsFiles'
                    , 'readFromWS', 'splitParams', 'externalConstr', 'runPeriods' ] :
            self[par] = getKWArg( self, { }, par )

        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix, getParNamePrefix
        setParNamePrefix( self['parNamePrefix'] )
        namePF = getParNamePrefix(True)

        self['condTagging'] = True if self['contEstWTag'] else self['condTagging']

        # create observables
        if self['runPeriods'] :
            states = dict( [ ( 'p%d' % per, per ) for per in self['runPeriods'] ] )
            self['obsDict']['runPeriod'] = tuple( [ comp if it != 2 else states for it, comp in enumerate(self['obsDict']['runPeriod']) ] )
        else :
            self['obsDict'].pop('runPeriod')

        states = dict( [ ( 'bin%d' % it, it ) for it in range( len( self['KKMassBinBounds'] ) - 1 ) ] )
        self['obsDict']['KKMassCat'] = tuple( [ comp if it != 2 else states for it, comp in enumerate( self['obsDict']['KKMassCat'] ) ] )

        self._createObservables()
        self._setObsProperties()
        self['obsSetP2VV']  = getKWArg( self, { }, 'obsSetP2VV' )
        self['observables'] = getKWArg( self, { }, 'observables' )

        # set (empty) dictionary of parameters
        self['parameters'] = { }

        # build tagging categories
        buildTaggingCategories( self, data = self['signalData'] )

        # build signal PDFs
        from P2VV.RooFitWrappers import Component
        self['pdfComps'] = [ ]
        self['pdfComps'].append( Component( namePF + 'sig', [ ], Yield = ( self['numEvents'] * self['sigFrac'], 0., self['numEvents'] ) ) )
        self['pdfComps'][0] += buildBs2JpsiphiSignalPdf(self)
        if self['tagPdfType'] or not self['condTagging'] : self['pdfComps'][0] += buildTaggingPdf( self, data = self['signalData'] )

        if not self['SFit'] : 
            # build background PDFs
            self['pdfComps'].append( Component( namePF + 'bkg', [ ]
                                    , Yield = ( self['numEvents'] * ( 1. - self['sigFrac'] ), 0., self['numEvents'] ) ) )
            self['pdfComps'][1] += buildBs2JpsiphiCombBkgTimePdf(self)
            self['pdfComps'][1] += buildBs2JpsiphiCombBkgAnglesPdf(self)
            if self['tagPdfType'] or not self['condTagging'] : self['pdfComps'][1] += buildTaggingPdf( self, data = self['signalData'] )

        # build full PDF
        from P2VV.RooFitWrappers import buildPdf
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: requesting %s PDF for observables [%s]'\
              % ( 'signal' if self['SFit'] else 'signal + background', ', '.join( str(obs) for obs in self['obsSetP2VV'] ) )
        self['fullPdf'] = buildPdf( self['pdfComps'], Observables = self['obsSetP2VV'], Name = 'Jpsiphi' )

        # build simultaneous PDF by splitting parameters
        self['pdf'] = self._createSimultaneous()

        # multiply by acceptance functions
        if self['timeEffType'] :   self._multiplyByTimeAcceptance()
        if self['anglesEffType'] : self._multiplyByAngularAcceptance()

        # create external constraints
        if self['externalConstr'] : self._createExternalConstraints()

        # collect python garbage
        import gc
        gc.collect()

        # check if no configuration arguments are left
        assert not self['kwargs'], 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: superfluous arguments found: %s' % self['kwargs']
        print 120 * '='


    def _setObsProperties( self, **kwargs ) :
        observables      = getKWArg( self, kwargs, 'observables' )
        obsDict          = getKWArg( self, kwargs, 'obsDict' )
        KKMassBinBounds  = getKWArg( self, kwargs, 'KKMassBinBounds' )
        paramKKMass      = getKWArg( self, kwargs, 'paramKKMass' )
        CSPValues        = getKWArg( self, kwargs, 'CSPValues' )
        tagPdfType       = getKWArg( self, kwargs, 'tagPdfType' )
        condTagging      = getKWArg( self, kwargs, 'condTagging' )
        SFit             = getKWArg( self, kwargs, 'SFit' )
        numTimeBins      = getKWArg( self, kwargs, 'numTimeBins' )
        numTimeResBins   = getKWArg( self, kwargs, 'numTimeResBins' )
        angularRanges    = getKWArg( self, kwargs, 'angularRanges' )

        # set observable properties
        observables['time'].setRange( 'Bulk', ( None, 5. ) )
        observables['time'].setBins(numTimeBins)
        observables['timeRes'].setBins(numTimeResBins)
        observables['timeRes'].setBins( numTimeResBins, 'cache' )
        observables['mass'].setRanges( dict(  LeftSideBand  = ( 5200., 5320. )
                                            , Signal        = ( 5320., 5420. )
                                            , RightSideBand = ( 5420., 5550. )
                                            , PeakBkg       = ( 5390., 5440. )
                                           )
                                     )

        for ran in angularRanges.get( 'cpsi',   [ ] ) : observables['cpsi'].setRange(   ran[0], ran[ 1 : 3 ] )
        for ran in angularRanges.get( 'ctheta', [ ] ) : observables['ctheta'].setRange( ran[0], ran[ 1 : 3 ] )
        for ran in angularRanges.get( 'phi',    [ ] ) : observables['phi'].setRange(    ran[0], ran[ 1 : 3 ] )

        # check KK mass bin parameters
        assert obsDict['KKMass'][4] == KKMassBinBounds[0] and obsDict['KKMass'][5] == KKMassBinBounds[-1]\
               , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: KK mass range in "KKMassBinBounds" is not the same as in "obsDict"'

        if paramKKMass in [ 'parameters', 'simultaneous' ] :
            assert len(CSPValues) == len(KKMassBinBounds) - 1,\
                   'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: wrong number of KK mass bin parameters specified'
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: KK mass bins: %s' % ' - '.join( '%.1f' % binEdge for binEdge in KKMassBinBounds )
        else :
            assert len(CSPValues) == 1,\
                   'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: only one S-P-wave coupling factor and no S-wave amplitude values should be specified'

        # get KK mass binning
        assert observables['KKMass'].hasBinning('KKMassBinning')\
               , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: KK mass observable does not have a binning named "KKMassBinning"'
        self['KKMassBinning'] = observables['KKMass'].getBinning('KKMassBinning')

        assert self['KKMassBinning'].numBins() == observables['KKMassCat'].numTypes() == len(KKMassBinBounds) - 1\
               , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: different numbers of bins in KK mass category (%s) and/or binning (%s) and specified bin bounds (%s)' \
               % ( self['KKMassBinning'].numBins(), observables['KKMassCat'].numTypes(), len(KKMassBinBounds) - 1 )
        for it in range( self['KKMassBinning'].numBins() ) :
            assert self['KKMassBinning'].binLow(it) == KKMassBinBounds[it]\
                   , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: different boundary in KK mass binning (%s) and specified bin bounds (%s)' \
                   % ( self['KKMassBinning'].binLow(it),  KKMassBinBounds[it] )

        # PDF observables set
        self['obsSetP2VV'] = [ observables[name] for name in [ 'time', 'cpsi', 'ctheta', 'phi' ] ]
        if tagPdfType or not condTagging :
            self['obsSetP2VV'].append( observables['iTagOS'] )
            if SSTagging : self['obsSetP2VV'].append( observables['iTagSS'] )
        if not SFit :
            self['obsSetP2VV'].append( observables['mass'] )


    def _createSimultaneous( self, **kwargs ) :
        observables       = getKWArg( self, kwargs, 'observables' )
        CSPValues         = getKWArg( self, kwargs, 'CSPValues' )
        contEstWTag       = getKWArg( self, kwargs, 'contEstWTag' )
        paramKKMass       = getKWArg( self, kwargs, 'paramKKMass' )
        ASParam           = getKWArg( self, kwargs, 'ASParam' )
        timeEffType       = getKWArg( self, kwargs, 'timeEffType' )
        timeEffHistFiles  = getKWArg( self, kwargs, 'timeEffHistFiles' )
        anglesEffType     = getKWArg( self, kwargs, 'anglesEffType' )
        angEffMomsFiles   = getKWArg( self, kwargs, 'angEffMomsFiles' )
        splitParams       = getKWArg( self, kwargs, 'splitParams' )
        fullPdf           = getKWArg( self, kwargs, 'fullPdf' )
        timeResModelsOrig = getKWArg( self, kwargs, 'timeResModelsOrig' )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        self['splitParsDict'] = { }

        # check split parameters
        if paramKKMass == 'simultaneous' :
            if 'KKMassCat' not in splitParams : splitParams['KKMassCat'] = [ 'C_SP' ]
            elif 'C_SP' not in splitParams['KKMassCat'] : splitParams['KKMassCat'].append('C_SP')
            for cat, pars in splitParams.iteritems() :
                if cat != 'KKMassCat' :
                    assert 'C_SP' not in pars\
                           , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: S-P coupling factors are set to be split for category "%s"' % cat

        if anglesEffType and type(angEffMomsFiles) != str :
            for cat in angEffMomsFiles.categories() :
                assert type(cat) == str\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: one of category names specified in "angEffMomsFiles" is not a string'
                if cat not in splitParams :
                    splitParams[cat] = [ 'angEffDummyCoef' ]
                elif 'angEffDummyCoef' not in splitParams[cat] :
                    splitParams[cat].append('angEffDummyCoef')

        if timeEffType and type(timeEffHistFiles) != dict :
            for cat in timeEffHistFiles.categories() :
                assert type(cat) == str\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: one of category names specified in "timeEffHistFiles" is not a string'
                assert cat in splitParams\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: PDF is not set to be split for category "%s" (required by "timeEffHistFiles")' % cat

        # split PDF for different data samples
        if splitParams :
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: splitting parameters in PDF "%s"' % fullPdf

            # get workspace and PDF variables
            ws   = fullPdf.ws()
            vars = fullPdf.getVariables()

            # get splitting categories and parameters
            splitParsDict = { }
            for cat, params in splitParams.iteritems() :
                assert ws.cat(cat), 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: category "%s" not in workspace' % cat
                for par in params :
                    par = namePF + par
                    assert ws.var(par), 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no variable "%s" in workspace' % par
                    assert vars.find(par), 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: variable "%s" not in PDF' % par
                    if ws[par] not in splitParsDict :
                        splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
                    else :
                        splitParsDict[ ws[par] ].add( ws[cat] )
            self['splitParsDict'] = splitParsDict.copy()

            # create lists of split categories and parameters
            pars = splitParsDict.keys()
            splitCats = [ ]
            splitPars = [ ]
            for par in pars :
                if par not in splitParsDict : continue
                splitPars.append( set( [par] ) )
                splitCats.append( splitParsDict.pop(par) )
                for par1 in pars :
                    if par1 not in splitParsDict : continue
                    if splitParsDict[par1] == splitCats[-1] :
                        splitPars[-1].add(par1)
                        splitParsDict.pop(par1)

            # build simultaneous PDF
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building simultaneous PDF "%s":' % ( fullPdf.GetName() + '_simul' )
            print 13 * ' ' + 'splitting categories: [ %s ]' % ' ], [ '.join(', '.join(cat.GetName() for cat in cats) for cats in splitCats)
            print 13 * ' ' + 'split parameters:     [ %s ]' % ' ], [ '.join(', '.join(par.GetName() for par in pars) for pars in splitPars)
            from P2VV.RooFitWrappers import SimultaneousPdf
            self['simulPdf'] = SimultaneousPdf(  fullPdf.GetName() + '_simul'
                                               , MasterPdf       = fullPdf
                                               , SplitCategories = splitCats
                                               , SplitParameters = splitPars
                                              )

            # set time resolution models
            splitCat      = self['simulPdf'].indexCat()
            splitCatPars  = self['simulPdf'].getVariables()
            splitCatIter  = splitCat.typeIterator()
            splitCatState = splitCatIter.Next()
            inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()
            origResParams = timeResModelsOrig['prototype']['model']._target_()._parameters
            from P2VV.Utilities.General import getSplitPar
            while splitCatState :
                splitCat.setIndex( splitCatState.getVal() )
                catPdf = self['simulPdf'].getPdf( splitCatState.GetName() )
                resModelCount = 0
                from ROOT import RooBTagDecay
                for comp in filter( lambda x : isinstance( x, RooBTagDecay ), catPdf.getComponents() ) :
                    # TODO: don't do this with RooBTagDecay, but move RooBTagDecay::resModel() to RooAbsAnaConvPdf
                    assert resModelCount < 1, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: multiple resolution models found for simultaneous category "%s"'\
                                              % splitCatState.GetName()
                    params = [ ]
                    for var in origResParams :
                        splitCats = self['splitParsDict'].get( ws[ var.GetName() ], set() )
                        if not splitCats :
                            params.append(var)
                        else :
                            catLabels = [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats if cat in splitCats ]
                            catsStr = ';'.join( lab[1] for lab in catLabels )
                            if len(catLabels) > 1 : catsStr = '{' + catsStr + '}'
                            splitVar = getSplitPar( var.GetName(), catsStr, splitCatPars )
                            assert splitVar, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: parameter "%s" is set to be constrainedi for category "%s", but it is not found in PDF'\
                                   % ( var.GetName(), catsStr )
                            from P2VV.RooFitWrappers import RealVar
                            params.append( RealVar( Name = splitVar.GetName() ) )

                    from P2VV.RooFitWrappers import ResolutionModel
                    from P2VV.Parameterizations.TimeResolution import TimeResolution
                    timeResModelsOrig[ splitCatState.GetName() ]\
                        = TimeResolution( Model = ResolutionModel( Name = comp.resolutionModel().GetName() ), Parameters = params )
                    resModelCount += 1

                assert resModelCount > 0, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no resolution model found for simultaneous category "%s"'\
                                          % splitCatState.GetName()

                splitCatState = splitCatIter.Next()

            if paramKKMass == 'simultaneous' and ASParam != 'Mag2ReIm' :
                # set values for splitted S-P coupling factors
                if ASParam != 'Mag2ReIm' :
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: using S-P-wave coupling factors:',
                    for iter, fac in enumerate( CSPValues ) :
                        print '%d: %.4f%s' % ( iter, fac, '' if iter == len( CSPValues ) - 1 else ',' ),
                    print

                splitCatIter.Reset()
                splitCatState = splitCatIter.Next()
                while splitCatState :
                    splitCat.setIndex( splitCatState.getVal() )
                    C_SP = getSplitPar( namePF + 'C_SP', observables['KKMassCat'].getLabel(), splitCatPars )
                    C_SP.setVal( CSPValues[ observables['KKMassCat'].getIndex() ] )
                    C_SP.setConstant(True)
                    splitCatState = splitCatIter.Next()

        else :
            self['simulPdf'] = None

        return self['simulPdf'] if self['simulPdf'] else self['fullPdf']


    def _multiplyByTimeAcceptance( self, **kwargs ) :
        timeEffHistFiles = getKWArg( self, kwargs, 'timeEffHistFiles' )
        pdf              = getKWArg( self, kwargs, 'pdf' )
        simulPdf         = getKWArg( self, kwargs, 'simulPdf' )
        signalData       = getKWArg( self, kwargs, 'signalData' )

        if not simulPdf :
            multiplyByTimeAcceptance( pdf, self, data = signalData, histFile = timeEffHistFiles['file']
                                     , histUBName = timeEffHistFiles['hlt1UB'], histExclBName = timeEffHistFiles['hlt1ExclB'] )
        else :
            splitCat      = simulPdf.indexCat()
            splitCatIter  = splitCat.typeIterator()
            splitCatState = splitCatIter.Next()
            inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()
            while splitCatState :
                splitCat.setIndex( splitCatState.getVal() )
                if type(timeEffHistFiles) == dict :
                    effFile = timeEffHistFiles
                else :
                    effFile = timeEffHistFiles.getSettings( [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats ] )
                catPdf = simulPdf.getPdf( splitCatState.GetName() )
                cNamePF = ( splitCatState.GetName() ).replace( '{', '' ).replace( '}', '' ).replace( ';', '_' )
                multiplyByTimeAcceptance( catPdf, self, data = signalData, histFile = effFile['file'], histUBName = effFile['hlt1UB']
                                         , histExclBName = effFile['hlt1ExclB'], coefNamePF = cNamePF, motherPdf = simulPdf
                                         , resModelKey = splitCatState.GetName() )
                splitCatState = splitCatIter.Next()


    def _multiplyByAngularAcceptance( self, **kwargs ) :
        angEffMomsFiles = getKWArg( self, kwargs, 'angEffMomsFiles' )
        pdf             = getKWArg( self, kwargs, 'pdf' )
        simulPdf        = getKWArg( self, kwargs, 'simulPdf' )

        if type(angEffMomsFiles) == str :
            multiplyByAngularAcceptance( pdf, self, angEffMomsFile = angEffMomsFiles )
        else :
            splitCat      = simulPdf.indexCat()
            splitCatIter  = splitCat.typeIterator()
            splitCatState = splitCatIter.Next()
            inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()
            while splitCatState :
                splitCat.setIndex( splitCatState.getVal() )
                effFile = angEffMomsFiles.getSettings( [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats ] )
                catPdf = simulPdf.getPdf( splitCatState.GetName() )
                cNamePF = ( splitCatState.GetName() ).replace( '{', '' ).replace( '}', '' ).replace( ';', '_' )
                multiplyByAngularAcceptance( catPdf, self, angEffMomsFile = effFile, coefNamePF = cNamePF )
                splitCatState = splitCatIter.Next()


    def _createExternalConstraints( self, **kwargs ) :
        externalConstr = getKWArg( self, kwargs, 'externalConstr' )
        splitParsDict  = getKWArg( self, kwargs, 'splitParsDict' )
        pdf            = getKWArg( self, kwargs, 'pdf' )
        simulPdf       = getKWArg( self, kwargs, 'simulPdf' )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        constraints = set()
        def buildConstraint( par, constrVals ) :
            from P2VV.RooFitWrappers import Pdf, ConstVar
            from ROOT import RooGaussian as Gaussian
            constrName = par.GetName().replace( '{', '' ).replace( '}', '' ).replace( ';', '_' )
            constraints.add( Pdf(  Name = constrName + '_constraint', Type = Gaussian
                                 , Parameters = [  par
                                                 , ConstVar( Name = constrName + '_mean',  Value = constrVals[0] )
                                                 , ConstVar( Name = constrName + '_sigma', Value = constrVals[1] )
                                                ]
                                )
                           )

        ws = pdf.ws()
        pdfVars = pdf.getVariables()
        if simulPdf :
            splitCat      = simulPdf.indexCat()
            splitCatIter  = splitCat.typeIterator()
            inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()

        for par, constrVals in externalConstr.iteritems() :
            parVar = ws[ namePF + par ]
            assert parVar, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: parameter "%s" is set to be constrained, but it is not found in workspace'\
                           % ( namePF + par )
            splitCats = splitParsDict.get( parVar, set() )
            if not splitCats :
                if type(constrVals) == SimulCatSettings : constrVals = constrVals.default()
                buildConstraint( parVar, constrVals )

            else :
                assert simulPdf\
                      , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: found splitting categories for parameter "%s", but no simultaneous PDF' % par

                splitCatIter.Reset()
                splitCatState = splitCatIter.Next()
                catsStrings = [ ]
                while splitCatState :
                    splitCat.setIndex( splitCatState.getVal() )
                    splitCatState = splitCatIter.Next()

                    catLabels = [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats if cat in splitCats ]
                    catsStr = ';'.join( lab[1] for lab in catLabels )
                    if len(catLabels) > 1 : catsStr = '{' + catsStr + '}'
                    if catsStr in catsStrings : continue

                    catsStrings.append(catsStr)
                    from P2VV.Utilities.General import getSplitPar
                    parVar = getSplitPar( namePF + par, catsStr, pdfVars )
                    assert parVar, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: parameter "%s" is set to be constrained, but it is not found in PDF'\
                                   % ( namePF + par )

                    if type(constrVals) == SimulCatSettings :
                        constrValsTuple = constrVals.getSettings(catLabels)
                    else :
                        constrValsTuple = constrVals
                    buildConstraint( parVar, constrValsTuple )

        pdf['ExternalConstraints'] = pdf['ExternalConstraints'] | constraints


###########################################################################################################################################
## Helper functions ##
######################

# function to build B_s^0 -> J/psi phi signal PDF
def buildBs2JpsiphiSignalPdf( self, **kwargs ) :
    """builds the signal PDF for the measurement of phi_s in B_s -> J/psi(->mu^+ mu^-) phi(->K^+ K^-)
    """

    # build parameters
    blind              = getKWArg( self, kwargs, 'blind' )
    obsDict            = getKWArg( self, kwargs, 'obsDict' )
    observables        = getKWArg( self, kwargs, 'observables' )
    tagCatsOS          = getKWArg( self, kwargs, 'tagCatsOS' )
    tagCatsSS          = getKWArg( self, kwargs, 'tagCatsSS' )
    transAngles        = getKWArg( self, kwargs, 'transAngles' )
    sigTaggingType     = getKWArg( self, kwargs, 'sigTaggingType' )
    ambiguityPars      = getKWArg( self, kwargs, 'ambiguityPars' )
    CSPValues          = getKWArg( self, kwargs, 'CSPValues' )
    SSTagging          = getKWArg( self, kwargs, 'SSTagging' )
    contEstWTag        = getKWArg( self, kwargs, 'contEstWTag' )
    condTagging        = getKWArg( self, kwargs, 'condTagging' )
    timeResType        = getKWArg( self, kwargs, 'timeResType' )
    constrainTResScale = getKWArg( self, kwargs, 'constrainTResScale' )
    paramKKMass        = getKWArg( self, kwargs, 'paramKKMass' )
    amplitudeParam     = getKWArg( self, kwargs, 'amplitudeParam' )
    ASParam            = getKWArg( self, kwargs, 'ASParam' )
    AparParam          = getKWArg( self, kwargs, 'AparParam' )
    constrainDeltaM    = getKWArg( self, kwargs, 'constrainDeltaM' )
    constrainBeta      = getKWArg( self, kwargs, 'constrainBeta' )
    lambdaCPParam      = getKWArg( self, kwargs, 'lambdaCPParam' )
    timeEffType        = getKWArg( self, kwargs, 'timeEffType' )
    anglesEffType      = getKWArg( self, kwargs, 'anglesEffType' )
    angEffMomsFiles    = getKWArg( self, kwargs, 'angEffMomsFiles' )
    parNamePrefix      = getKWArg( self, kwargs, 'parNamePrefix', '' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix, getParNamePrefix
    setParNamePrefix(parNamePrefix)
    namePF = getParNamePrefix(True)

    print 120 * '='
    print 'P2VV - INFO: buildBs2JpsiphiSignalPdf(): building B_s^0 -> J/psi phi signal PDF'

    # angular functions
    if transAngles : from P2VV.Parameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
    else :           from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles     as AngleFuncs
    self['angleFuncs'] = AngleFuncs( cpsi = observables['cpsi'], ctheta = observables['ctheta'], phi = observables['phi']
                                    , DummyCoef = ( anglesEffType and type(angEffMomsFiles) != str ) )

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
        self['amplitudes'] = Amplitudes( ASParameterization = ASParam, AparParameterization = AparParam, **commonArgs )

    elif amplitudeParam == 'phases' :
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet as Amplitudes
        self['amplitudes'] = Amplitudes( ASParameterization = ASParam, **commonArgs )

    elif amplitudeParam == 'bank' :
        from P2VV.Parameterizations.DecayAmplitudes import JpsiVBank_AmplitudeSet as Amplitudes
        self['amplitudes'] = Amplitudes( ASParameterization = ASParam, AparParameterization = AparParam, **commonArgs )

    else :
        raise RuntimeError('P2VV - ERROR: buildBs2JpsiphiSignalPdf(): no valid amplitude parameterization specified')

    # B decay time
    from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
    dGammaVar = dict( Name = 'dGamma' )
    if blind and 'dGamma' in blind :
        if blind['dGamma'] : dGammaVar['Blind'] = blind['dGamma']
        else :               dGammaVar['Blind'] = ( 'UnblindUniform', 'BsDGs2013EPS', 0.02 )
    self['lifetimeParams'] = LifetimeParams( dGamma = dGammaVar, dMConstraint = constrainDeltaM, betaConstraint = constrainBeta )
    if ambiguityPars : self['lifetimeParams']['dGamma'].setVal( -self['lifetimeParams']['dGamma'].getVal() )

    if timeResType.startswith('event') :
        timeResArgs = dict( time = observables['time'], timeResSigma = observables['timeRes'], Cache = False if timeEffType else True )
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
            from P2VV.RooFitWrappers import ConstVar
            timeResArgs['timeResMean']   = ConstVar( Name = namePF + 'timeResMean',   Value = 0. )
            timeResArgs['timeResMeanSF'] = ConstVar( Name = namePF + 'timeResMeanSF', Value = 1. )
            timeResArgs['timeResSFConstraint'] = constrainTResScale
        elif 'constmean' in timeResType.lower() :
            from P2VV.RooFitWrappers import ConstVar
            timeResArgs['timeResMean']   = dict( Value = -0.01, Error = 0.005 )
            timeResArgs['timeResMeanSF'] = ConstVar( Name = namePF + 'timeResMeanSF', Value = 1. )
            timeResArgs['timeResMeanConstraint'] = constrainTResScale
            timeResArgs['timeResSFConstraint'] = constrainTResScale
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
            timeResArgs['timeResMeanConstraint'] = constrainTResScale
            timeResArgs['timeResSFConstraint'] = constrainTResScale

        from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
        resModel = TimeResolution( **timeResArgs )

    else :
        from P2VV.Parameterizations.TimeResolution import Gaussian_TimeResolution as TimeResolution
        resModel = TimeResolution(  time          = observables['time']
                                  , timeResMu     = dict( Value = 0.,    Constant = True )
                                  , timeResSigma  = dict( Value = 0.045, Constant = True )
                                  , PerEventError = False
                                  , Cache         = False if timeEffType else True
                                 )

    print 'P2VV - INFO: buildBs2JpsiphiSignalPdf(): decay time resolution model:'
    resModel['model'].Print()
    for par in resModel.parameters() : par.Print()
    self['timeResModels']     = dict( prototype = resModel )
    self['timeResModelsOrig'] = dict( prototype = resModel )

    # CP violation parameters
    if lambdaCPParam == 'ReIm' : 
        from P2VV.Parameterizations.CPVParams import LambdaCarth_CPParam as CPParam
        ReLambdaCPVar = dict( Name = 'ReLambdaCP' )
        ImLambdaCPVar = dict( Name = 'ImLambdaCP' )
        self['lambdaCP'] = CPParam( ReLambdaCP = ReLambdaCPVar, ImLambdaCP = ImLambdaCPVar )

    elif lambdaCPParam == 'lambPhiRel_CPVDecay' :
        from P2VV.Parameterizations.CPVParams import LambdaAbsArgRel_CPVDecay_CPParam as CPParam
        phiCPVar = dict( Name = 'phiCP_m' )
        self['lambdaCP'] = CPParam( phiCP_m = phiCPVar, AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self['amplitudes'] )
        if ambiguityPars : self['lambdaCP']['phiCP_m'].setVal( pi - self['lambdaCP']['phiCP_m'].getVal() )

    elif lambdaCPParam.startswith('lambPhi_CPVDecay') :
        from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPVDecay_CPParam as CPParam
        if lambdaCPParam.endswith('PSWaves') :
            from ROOT import RooNumber
            RooInf = RooNumber.infinity()
            rhoCP_P = RealVar( namePF + 'rhoCP_P', Title = 'CPV in decay |rho|', Value = 1., Error = 0.04, MinMax = ( 0., 5. ) )
            rhoCP_S = RealVar( namePF + 'rhoCP_S', Title = 'CPV in decay |rho|', Value = 1., Error = 0.04, MinMax = ( 0., 5. ) )
            phiCP_P = RealVar( namePF + 'phiCP_P', Title = 'CPV in decay phi',   Value = 0., Error = 0.1,  MinMax = (-RooInf, +RooInf) )
            phiCP_S = RealVar( namePF + 'phiCP_S', Title = 'CPV in decay phi',   Value = 0., Error = 0.1,  MinMax = (-RooInf, +RooInf) )
            self['lambdaCP'] = CPParam(  AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self['amplitudes']
                                     , rhoCP_A0 = rhoCP_P, rhoCP_Apar = rhoCP_P, rhoCP_Aperp = rhoCP_P, rhoCP_AS = rhoCP_S
                                     , phiCP_A0 = phiCP_P, phiCP_Apar = phiCP_P, phiCP_Aperp = phiCP_P, phiCP_AS = phiCP_S
                                    )

        else :
            self['lambdaCP'] = CPParam( AmplitudeNames = [ 'A0', 'Apar', 'Aperp', 'AS' ], Amplitudes = self['amplitudes'] )

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
        self['lambdaCP'] = CPParam( phiCP = phiCPVar, lambdaCP = lambdaCPVar )
        if ambiguityPars :
            self['lambdaCP']['phiCP'].setVal( pi - self['lambdaCP']['phiCP'].getVal() )

    # coefficients for time functions
    from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
    timeBasisCoefs = TimeBasisCoefs( self['angleFuncs'].functions, self['amplitudes'], self['lambdaCP'], [ 'A0', 'Apar', 'Aperp', 'AS' ] )

    # tagging parameters
    if sigTaggingType == 'histPdf' :
        raise RuntimeError('P2VV - ERROR: buildBs2JpsiphiSignalPdf(): A histogram tagging PDF can only be used when tagging observables are conditional')

    else :
        # get tagging category parameters dictionary/dictionaries
        tagCatsDictOS = tagCatsOS.tagCatsDict()
        if SSTagging : tagCatsDictSS = self['tagCatsSS'].tagCatsDict()

        if sigTaggingType.startswith('tagUntag') :
            # assume products of asymmetries are small and B-Bbar asymmetries are equal for all tagged categories
            if observables['tagCatOS'].numTypes() > 2 or ( SSTagging and observables['tagCatSS'].numTypes() > 2 ) :
                print 'P2VV - INFO: buildBs2JpsiphiSignalPdf(): tagging in signal PDF:\n'\
                    + '    * assuming B-Bbar asymmetries are equal for all tagged categories'

            # provide the same asymmetry for all tagged categories
            from math import sqrt
            asymVal = 0. #asymVal = -self['lambdaCP']['C'].getVal()
            asymErr = 0.1
            from P2VV.RooFitWrappers import RealVar
            avgCEvenSum = RealVar( namePF + 'avgCEvenSum', Title = 'Sum of CP average even coefficients'
                                  , Value = 1., MinMax = (  0., 2. ) )
            avgCOddSum  = RealVar( namePF + 'avgCOddSum', Title = 'Sum of CP average odd coefficients'
                                  , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )
            avgCEvenOS  = RealVar( namePF + 'avgCEvenOSTagged', Title = 'CP average even coefficients OS tagged categories'
                                  , Value = 1., MinMax = (  0., 2. ) )
            avgCOddOS   = RealVar( namePF + 'avgCOddOSTagged' , Title = 'CP average odd coefficients OS tagged categories'
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
                from P2VV.RooFitWrappers import RealVar
                avgCEvenSS = RealVar( namePF + 'avgCEvenSSTagged', Title = 'CP average even coefficients SS tagged categories'
                                     , Value = 1., MinMax = (  0., 2. ) )
                avgCOddSS  = RealVar( namePF + 'avgCOddSSTagged', Title = 'CP average odd coefficients SS tagged categories'
                                     , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ) )
                avgCEven   = RealVar( namePF + 'avgCEvenTagged', Title = 'CP average even coefficients OS+SS tagged categories'
                                     , Value = 1., MinMax = (  0., 2. ) )
                avgCOdd    = RealVar( namePF + 'avgCOddTagged', Title = 'CP average odd coefficients OS+SS tagged categories'
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
            tagCatsDict['ANorm'] = 0. #-self['lambdaCP']['C'].getVal()

        from P2VV.Parameterizations.FlavourTagging import CatDilutionsCoefAsyms_TaggingParams as TaggingParams
        self['taggingParams'] = TaggingParams( **tagCatsDict )

        if condTagging :
            # don't float tagging category coefficients if PDF is conditional on tagging observables
            for coefList in self['taggingParams']['singleTagCatCoefs'] :
                for coef in coefList :
                    if coef and coef.isFundamental() : coef.setConstant(True)

        if not SSTagging :
            args = dict(  tagCat      = observables['tagCatOS']
                        , iTag        = observables['iTagOS']
                        , dilutions   = self['taggingParams']['dilutions']
                        , ADilWTags   = self['taggingParams']['ADilWTags']
                        , avgCEvens   = self['taggingParams']['avgCEvens']
                        , avgCOdds    = self['taggingParams']['avgCOdds']
                        , tagCatCoefs = self['taggingParams']['tagCatCoefs']
                       )
        else :
            args = dict(  tagCat0     = observables['tagCatOS']
                        , tagCat1     = observables['tagCatSS']
                        , iTag0       = observables['iTagOS']
                        , iTag1       = observables['iTagSS']
                        , dilutions0  = self['taggingParams']['dilutions'][0]
                        , dilutions1  = self['taggingParams']['dilutions'][1]
                        , ADilWTags0  = self['taggingParams']['ADilWTags'][0]
                        , ADilWTags1  = self['taggingParams']['ADilWTags'][1]
                        , avgCEvens   = self['taggingParams']['avgCEvens']
                        , avgCOdds    = self['taggingParams']['avgCOdds']
                        , tagCatCoefs = self['taggingParams']['tagCatCoefs']
                       )

    args.update(  time                   = observables['time']
                , tau                    = self['lifetimeParams']['MeanLifetime']
                , dGamma                 = self['lifetimeParams']['dGamma']
                , dm                     = self['lifetimeParams']['dM']
                , coshCoef               = timeBasisCoefs['cosh']
                , sinhCoef               = timeBasisCoefs['sinh']
                , cosCoef                = timeBasisCoefs['cos']
                , sinCoef                = timeBasisCoefs['sin']
                , resolutionModel        = self['timeResModels']['prototype']['model']
                , ConditionalObservables = self['amplitudes'].ConditionalObservables()\
                                           | self['timeResModels']['prototype'].ConditionalObservables()\
                                           | self['taggingParams'].ConditionalObservables()
                , ExternalConstraints    = self['lifetimeParams'].ExternalConstraints()\
                                           | self['timeResModels']['prototype'].ExternalConstraints()\
                                           | self['taggingParams'].ExternalConstraints()
               )

    # build signal PDF
    from P2VV.RooFitWrappers import BTagDecay
    self['sigBTagDecay'] = BTagDecay( namePF + ( 'sig_t_angles' if condTagging else 'sig_t_angles_tagCat_iTag' ), **args )

    # collect python garbage
    import gc
    gc.collect()

    return self['sigBTagDecay']


# function to build B_s^0 -> J/psi phi background decay time PDF
def buildBs2JpsiphiCombBkgTimePdf( self, **kwargs ) :
    observables    = getKWArg( self, kwargs, 'observables' )
    timeResModels  = getKWArg( self, kwargs, 'timeResModels' )
    parNamePrefix  = getKWArg( self, kwargs, 'parNamePrefix', '' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
    setParNamePrefix(parNamePrefix)

    from P2VV.Parameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
    self['backgroundTime'] = BackgroundTime( Name = namePF + 'bkg_t', time = observables['time']
                                            , resolutionModel = timeResModels['prototype']['model'] )
    self['bkgTimePdf'] = self['backgroundTime'].pdf()
    return self['bkgTimePdf']


# function to build B_s^0 -> J/psi phi background decay angle PDF
def buildBs2JpsiphiCombBkgAnglesPdf( self, **kwargs ) :
    observables     = getKWArg( self, kwargs, 'observables' )
    bkgAnglePdfType = getKWArg( self, kwargs, 'bkgAnglePdfType' )
    cbkgData        = getKWArg( self, kwargs, 'cbkgData', None )
    transAngles     = getKWArg( self, kwargs, 'transAngles' )
    numAngleBins    = getKWArg( self, kwargs, 'numAngleBins' )
    fitOptions      = getKWArg( self, kwargs, 'fitOptions', { } )
    parNamePrefix   = getKWArg( self, kwargs, 'parNamePrefix', '' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix, getParNamePrefix
    setParNamePrefix(parNamePrefix)
    namePF = getParNamePrefix(True)

    if cbkgData and bkgAnglePdfType == 'histPdf' :
        # create a histogram from background data and use it for background angular PDF
        print 'P2VV - INFO: buildBs2JpsiphiCombBkgAnglesPdf(): determining angular shape of background from background data'\

        from P2VV.RooFitWrappers import HistPdf
        self['bkg_angles'] = HistPdf(  Name = 'bkg_angles'
                                   , Observables = [ observables['cpsi'], observables['ctheta'], observables['phi'] ]
                                   , Binning = { cpsi : numAngleBins[0], ctheta : numAngleBins[1], phi : numAngleBins[2] }
                                   , Data = cbkgData
                                  )

    else :
        print 'P2VV - INFO: buildBs2JpsiphiCombBkgAnglesPdf(): building a %s angular PDF for background'\
              % ( 'hybrid' if bkgAnglePdfType == 'hybrid' else 'basis functions' if bkgAnglePdfType == 'basis' else 'binned' )

        if bkgAnglePdfType != 'basis' :
            # create a binned PDF for background angular shape
            from array import array
            from math import pi
            baselineBin = -1
            if transAngles :
                nBins = [ 5, 7, 9 ]
                cpsBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 5. * float(i)   for i in range( 1, 5 ) ] + [ 1. ] )
                cthBinBounds = array( 'd', [ -1. ] + [        -1. + 2. / 7. * float(i)   for i in range( 1, 7 ) ] + [ 1. ] )
                phiBinBounds = array( 'd', [ -pi ] + [ pi * ( -1. + 2. / 9. * float(i) ) for i in range( 1, 9 ) ] + [ pi ] )

            elif bkgAnglePdfType == 'hybrid' :
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
            observables['cpsi'].setBinning(   cpsBins, 'bkg_cpsBins' )
            observables['ctheta'].setBinning( cthBins, 'bkg_cthBins' )
            observables['phi'].setBinning(    phiBins, 'bkg_phiBins' )

            # create bin coefficients
            from P2VV.RooFitWrappers import RealVar
            self['bkgAngCoefs'] = [ RealVar(  namePF + 'bkg_angBin_%d_%d_%d' % ( bin0, bin1, bin2 )
                                            , Title    = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                                            , Value    = 1. / cpsNumBins / cthNumBins / phiNumBins
                                            , MinMax   = ( 0., 1. )
                                            , Constant = False if bkgAnglePdfType == 'hybrid' and bin1 != baselineBin else True
                                           )\
                                    if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                                    for bin2 in range( phiNumBins ) for bin1 in range( cthNumBins ) for bin0 in range( cpsNumBins )
                                  ]

            from P2VV.RooFitWrappers import ComplementCoef
            self['bkgAngCoefs'][0] = ComplementCoef( Name = namePF + 'bkg_angBin_0_0_0', Coefficients = self['bkgAngCoefs'][ 1 : ] )

            # create a BinnedPdf
            from P2VV.RooFitWrappers import BinnedPdf
            self['bkg_angBins'] = BinnedPdf(  Name = namePF + ( 'bkg_angBins' if bkgAnglePdfType == 'hybrid' else 'bkg_angles' )
                                            , Observables = [ observables['cpsi'], observables['ctheta'], observables['phi'] ]
                                            , Binnings = [ cpsBins, cthBins, phiBins ]
                                            , Coefficients = self['bkgAngCoefs']
                                            , BinIntegralCoefs = True
                                           )
            self['bkg_angBins'].setForceUnitIntegral(True)

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
                    cpsBin = observables['cpsi'].getBin('bkg_cpsBins')
                    cthBin = observables['ctheta'].getBin('bkg_cthBins')
                    phiBin = observables['phi'].getBin('bkg_phiBins')
                    bin = cpsBin + cthBin * cpsNumBins + phiBin * cpsNumBins * cthNumBins - 1
                    if bin >= 0 : sumBinWeights[ bin ] += cbkgData.weight()
                for name, val in zip( [ 'cpsi', 'ctheta', 'phi' ], angleInitVals ) : observables[name].setVal(val)

                if bkgAnglePdfType == 'hybrid' :
                    baseBinVal = sumBinWeights[ baselineBin - 1 ] / sumWeights / cthBins.binWidth(baselineBin)
                    spikesFrac -= baseBinVal * 2.

                # set bin coefficient values
                #binCoefVals = [ 0.10, 0.07, 0., 0.07, 0.10, 0.37 ]
                #binCoefErrs = [ 0.03, 0.03, 0., 0.03, 0.03, 0.03 ]
                for bin, ( coef, weight ) in enumerate( zip( self['bkgAngCoefs'][ 1 : ], sumBinWeights ) ) :
                    binInt = weight / sumWeights
                    assert binInt >= 0.,\
                           'P2VV - INFO: buildBs2JpsiphiCombBkgAnglesPdf(): background angular PDF coefficient \"%s\" has negative value: %f'\
                           % (coef.GetName(), binInt)
                    if bkgAnglePdfType == 'hybrid' :
                        coef.setVal( ( ( binInt - baseBinVal * cthBins.binWidth(bin + 1) ) / spikesFrac )\
                                     if bin != baselineBin - 1 else 0. )
                        #coef.setVal( binCoefVals[bin] )
                        #coef.setError( binCoefErrs[bin] )
                    else :
                        coef.setVal(binInt)

            if bkgAnglePdfType == 'hybrid' :
                self['bkgAngBinsFrac'] = RealVar(  namePF + 'bkgAngBinsFrac'
                                                 , Title  = 'Binned PDF fraction in angular background'
                                                 , Value  = spikesFrac
                                                 , MinMax = ( 0., 1. )
                                                )

        if bkgAnglePdfType in [ 'basis', 'hybrid' ] :
            if bkgAnglePdfType == 'basis' : nBins = [ 20, 20, 20 ]
            # create an angular PDF with basis functions
            angPDFIndices = [ ( 2, 0, 0 ), ( 0, 2, 0 ) ]
            #angPDFIndices = [ ( 0, 2, 0 ), ( 0, 2, 2 ), ( 2, 0, 0 ), ( 2, 2, 0 ), ( 2, 2, 2 ) ]
            #angPDFIndices = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3)\
            #                  for YIndex1 in range( -YIndex0, YIndex0 + 1 )\
            #                  if PIndex != 0 or YIndex0 != 0 or YIndex1 != 0 ]

            from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms as angularPdfTerms
            cnvrtInd = lambda ind : 'm' + str(abs(ind)) if ind < 0 else str(ind)
            angPdfTerms = angularPdfTerms(  Angles = dict(  cpsi   = observables['cpsi']
                                                          , ctheta = observables['ctheta']
                                                          , phi    = observables['phi']
                                                         )
                                          , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                     , {  'Name'    : namePF + 'bkg_angCoef_%s_%s_%s'\
                                                                      % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                        , 'Value'   : 0.
                                                        , 'Error'   : 0.01
                                                        , 'MinMax'  : ( -0.4, 0.4 )
                                                        , 'Indices' : inds
                                                       }
                                                    ) for inds in angPDFIndices
                                                  )
                                         )
            self['bkg_angFuncs'] = angPdfTerms.buildSumPdf( namePF + 'bkg_angFuncs' )

        # build total angular PDF
        if bkgAnglePdfType == 'hybrid' :
            from P2VV.RooFitWrappers import SumPdf
            self['bkg_angles'] = SumPdf(  Name   = namePF + 'bkg_angles'
                                        , PDFs   = [ self['bkg_angBins'], self['bkg_angFuncs'] ]
                                        , Yields = { self['bkg_angBins'].GetName() : self['bkgAngBinsFrac'] }
                                       )

        elif bkgAnglePdfType == 'basis' :
            self['bkg_angles'] = self['bkg_angFuncs']

        else :
            self['bkg_angles'] = self['bkg_angBins']

        if cbkgData :
            # fit background angular distribution
            print 'P2VV - INFO: buildBs2JpsiphiCombBkgAnglesPdf(): fitting background angular distribution'
            print '  initial parameter values:'
            for par in self['bkg_angles'].getParameters(cbkgData) : par.Print()
            self['bkg_angles'].fitTo( cbkgData, SumW2Error = False, Save = False, **fitOptions )

    return self['bkg_angles']


# function to build a PDF for tagging observables
def buildTaggingPdf( self, **kwargs ) :
    data          = getKWArg( self, kwargs, 'data', None )
    observables   = getKWArg( self, kwargs, 'observables' )
    tagPdfType    = getKWArg( self, kwargs, 'tagPdfType' )
    SSTagging     = getKWArg( self, kwargs, 'SSTagging' )
    tagCatCoefs   = getKWArg( self, kwargs, 'tagCatCoefs', None )
    dilutions     = getKWArg( self, kwargs, 'dilutions', None )
    parNamePrefix = getKWArg( self, kwargs, 'parNamePrefix', '' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix, getParNamePrefix
    setParNamePrefix(parNamePrefix)
    namePF = getParNamePrefix(True)

    if data and pdfType == 'histPdf' :
        # create histogram from signal data and use the (fixed) bin coefficients for the PDF
        print 'P2VV - INFO: buildTaggingPdf(): creating tagging PDF from data'
        from P2VV.RooFitWrappers import HistPdf
        if not SSTagging :
            pdf = HistPdf(  Name = namePF + 'tagCat_iTag'
                          , Observables = [ observables['tagCatOS'], observables['iTagOS'] ]
                          , Data = data
                         )
        else :
            pdf = HistPdf(  Name = namePF + 'tagCat_iTag'
                          , Observables = [observables['tagCatOS'], observables['tagCatSS'], observables['iTagOS'], observables['iTagSS']]
                          , Data = data
                         )

    else :
        # use a PDF with variable bin coefficients
        if pdfType.startswith('tagUntag') or ( observables['tagCatOS'].numTypes() == 2\
                                               and ( not SSTagging or observables['tagCatSS'].numTypes() == 2 ) ) :
            # assume B-Bbar asymmetry is equal for all tagged categories
            if observables['tagCatOS'].numTypes() > 2 or ( SSTagging and observables['tagCatSS'].numTypes() > 2 ) :
                print 'P2VV - INFO: buildTaggingPdf(): tagging PDF:\n'\
                    + '    * assuming B-Bbar asymmetries are equal for all tagged categories'
            else :
                print 'P2VV - INFO: buildTaggingPdf(): building a binned signal tagging PDF'

            from P2VV.Parameterizations.FlavourTagging import TagUntag_BinnedTaggingPdf as TaggingPdf

        else :
            # create independent tagging bin coefficients
            if data :
                print 'P2VV - INFO: buildTaggingPdf(): determining signal tagging coefficients from signal data'
            else :
                print 'P2VV - INFO: buildTaggingPdf(): WARNING:'\
                    + ' no signal data available to determine signal tagging coefficients:\n'\
                    + '    * assuming absence of B-Bbar asymmetries'

            from P2VV.Parameterizations.FlavourTagging import TagCats_BinnedTaggingPdf as TaggingPdf

        # build PDF
        self['taggingPdf'] = TaggingPdf(  namePF + 'tagCat_iTag'
                                        , observables['tagCatOS'], observables['tagCatSS'] if SSTagging else None
                                        , observables['iTagOS'], observables['iTagSS'] if SSTagging else None
                                        , TagCatCoefs   = tagCatCoefs
                                        , TaggedCatName = 'TagCat' if observables['tagCatOS'].numTypes() > 2 else 'Tagged'
                                        , Data          = data
                                        , RelativeCoefs = True if pdfType.endswith('Relative') else False
                                        , Dilutions     = dilutions
                                       )
        pdf = self['taggingPdf']

    return pdf


# function to construct tagging categories
def buildTaggingCategories( self, **kwargs ) :
    observables      = getKWArg( self, kwargs, 'observables' )
    obsDict          = getKWArg( self, kwargs, 'obsDict' )
    tagCatsOS        = getKWArg( self, kwargs, 'tagCatsOS', [ ] )
    tagCatsSS        = getKWArg( self, kwargs, 'tagCatsSS', [ ] )
    constrainTagging = getKWArg( self, kwargs, 'constrainTagging' )
    tagPdfType       = getKWArg( self, kwargs, 'tagPdfType' )
    condTagging      = getKWArg( self, kwargs, 'condTagging' )
    contEstWTag      = getKWArg( self, kwargs, 'contEstWTag' )
    data             = getKWArg( self, kwargs, 'data', None )
    parNamePrefix    = getKWArg( self, kwargs, 'parNamePrefix', '' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
    setParNamePrefix(parNamePrefix)

    print 120 * '='
    print 'P2VV - INFO: buildTaggingCategories(): building tagging categories'

    if data :
        # get category bins
        assert observables['wTagOS'].hasBinning('tagCats'),\
               'P2VV - ERROR: buildTaggingCategories(): binning "tagCats" not found for OS estimated wrong-tag probability'
        assert observables['wTagSS'].hasBinning('tagCats'),\
               'P2VV - ERROR: buildTaggingCategories(): binning "tagCats" not found for SS estimated wrong-tag probability'
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
                assert cat.isValidIndex(it), 'P2VV - ERROR: buildTaggingCategories(): no bin %d found for tagging category "%s" ' % ( it, cat.GetName() )
            binPars = [ ( cat.lookupType(it).GetName(), it, bins.binHigh( bins.numBins()-it-1 ) ) for it in range( bins.numBins() ) ]

            from P2VV.Parameterizations.FlavourTagging import getTagCatParamsFromData as getTagParams
            tagBins.append( getTagParams( data, estWTagName = wTag, tagCats = binPars, numSigmas = 1., SameSide = isSS ) )

        tagCatsOS = tagBins[0]
        tagCatsSS = tagBins[1]

    else :
        assert tagCatsOS and tagCatsSS, 'P2VV - ERROR: buildTaggingCategories(): no tagging category binnings found'

    # build tagging categories
    from P2VV.Parameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
    self['tagCatsOS'] = TaggingCategories(  tagCat = observables['tagCatOS'] if observables['tagCatOS'] else obsDict['tagCatOS'][0]
                                          , estWTag = observables['wTagOS'] if contEstWTag else None
                                          , estWTagName = obsDict['wTagOS'][0], TagCats = tagCatsOS, SameSide = False
                                          , wTagP0Constraint = constrainTagging, wTagP1Constraint = constrainTagging
                                         )
    self['tagCatsSS'] = TaggingCategories(  tagCat = observables['tagCatSS'] if observables['tagCatSS'] else obsDict['tagCatSS'][0]
                                          , estWTag = observables['wTagSS'] if contEstWTag else None
                                          , estWTagName = obsDict['wTagSS'][0], TagCats = tagCatsSS, SameSide = True
                                          , wTagP0Constraint = constrainTagging, wTagP1Constraint = constrainTagging
                                         )

    observables['tagCatOS'] = self['tagCatsOS']['tagCat']
    observables['tagCatSS'] = self['tagCatsSS']['tagCat']

    if tagPdfType or not condTagging :
        self['obsSetP2VV'].append( observables['tagCatOS'] )
        if SSTagging : self['obsSetP2VV'].append( observables['tagCatSS'] )

    if condTagging :
        self['tagCatsOS'].addConditional( observables['tagCatOS'] )
        self['tagCatsOS'].addConditional( observables['iTagOS']   )
        self['tagCatsSS'].addConditional( observables['tagCatSS'] )
        self['tagCatsSS'].addConditional( observables['iTagSS']   )

    return ( self['tagCatsOS'], self['tagCatsSS'] )


# function to multiply a PDF by decay time acceptance function
def multiplyByTimeAcceptance( pdf, self, **kwargs ) :
    histFile          = getKWArg( self, kwargs, 'histFile' )
    histExclBName     = getKWArg( self, kwargs, 'histExclBName' )
    histUBName        = getKWArg( self, kwargs, 'histUBName' )
    data              = getKWArg( self, kwargs, 'data' )
    motherPdf         = getKWArg( self, kwargs, 'motherPdf', pdf )
    observables       = getKWArg( self, kwargs, 'observables' )
    timeEffParameters = getKWArg( self, kwargs, 'timeEffParameters' )
    timeEffType       = getKWArg( self, kwargs, 'timeEffType' )
    parNamePrefix     = getKWArg( self, kwargs, 'parNamePrefix', '' )
    coefNamePF        = getKWArg( self, kwargs, 'coefNamePF', '' )
    timeResModels     = getKWArg( self, kwargs, 'timeResModels' )
    timeResModelsOrig = getKWArg( self, kwargs, 'timeResModelsOrig' )
    resModelKey       = getKWArg( self, kwargs, 'resModelKey', 'prototype' )

    from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
    setParNamePrefix( parNamePrefix + ( '_' if parNamePrefix else '' ) + coefNamePF )

    time      = observables['time']
    hlt1ExclB = observables['hlt1ExclB']
    hlt2B     = observables['hlt2B']
    hlt2UB    = observables['hlt2UB']

    # build new decay-time resolution model that includes the decay-time acceptance function
    print 'P2VV - INFO: multiplyByTimeAcceptance(): multiplying PDF "%s" with decay-time acceptance function' % pdf.GetName()
    if timeEffType == 'fit' :
        hists = {  hlt1ExclB : {  'exclB'    : { 'histogram' : 'hlt1_shape', 'average' : ( 6.285e-01, 1.633e-02 ) }
                                , 'notExclB' : { 'bins'      : time.getRange(), 'heights' : [0.5]                 }
                               }
                 , hlt2B     : { 'B'         : { 'histogram' : 'hlt2_shape', 'average' : ( 6.3290e-01, 1.65e-02 ) } }
                 , hlt2UB    : { 'UB'        : { 'bins'      : time.getRange(), 'heights' : [0.5]                 } }
                }

        from P2VV.Parameterizations.TimeAcceptance import Paper2012_mer_TimeAcceptance as TimeAcceptance
        timeResModels[resModelKey] = TimeAcceptance(  time = time
                                                    , ResolutionModel = timeResModelsOrig[resModelKey]
                                                    , Input = histFile
                                                    , Histograms = hists
                                                    , Data = data
                                                    , Fit = True
                                                    , Original = motherPdf
                                                    , **timeEffParameters
                                                   )

    elif timeEffType == 'paper2012_multi' :
        hists = { hlt1ExclB : {  'exclB'    : { 'histogram' : histExclBName }
                               , 'notExclB' : { 'histogram' : histUBName    }
                              }
                }
        from P2VV.Parameterizations.TimeAcceptance import Paper2012_mer_TimeAcceptance as TimeAcceptance
        timeResModels[resModelKey] = TimeAcceptance(  time = time
                                                    , ResolutionModel = timeResModelsOrig[resModelKey]
                                                    , Input = histFile
                                                    , Histograms = hists
                                                    , Data = data
                                                    , Fit = False
                                                    , Original = motherPdf
                                                    , BinHeightMinMax = ( -RooInf, RooInf )
                                                    , **timeEffParameters
                                                   )

    elif timeEffType == 'paper2012' :
        hists = { hlt1ExclB : {  'exclB'    : { 'histogram' : histExclBName }
                               , 'notExclB' : { 'histogram' : histUBName    }
                              }
                }
        from P2VV.Parameterizations.TimeAcceptance import Paper2012_csg_TimeAcceptance as TimeAcceptance
        timeResModels[resModelKey] = TimeAcceptance(  time = time
                                                    , ResolutionModel = timeResModelsOrig[resModelKey]
                                                    , Input = histFile
                                                    , Histograms = hists
                                                    , **timeEffParameters
                                                   )


    elif timeEffType == 'HLT1Unbiased' or timeEffType == 'HLT1ExclBiased' :
        from P2VV.Parameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
        timeResModels[resModelKey] = TimeAcceptance(  time = time
                                                    , ResolutionModel = timeResModelsOrig[resModelKey]
                                                    , Input = histFile
                                                    , Histogram = histExclBName if timeEffType == 'HLT1ExclBiased'\
                                                                  else histUBName
                                                    , **timeEffParameters
                                                   )
    else:
        raise ValueError( 'P2VV - ERROR: multiplyByTimeAcceptance(): unknown time efficiency type: "%s"' % timeEffType )

    setParNamePrefix(parNamePrefix)

    # multiply PDF with time acceptance
    accPdfs = [ ]
    from ROOT import RooBTagDecay
    for comp in filter( lambda x : type(x) is RooBTagDecay, pdf.getComponents() ) :
        # TODO: don't do this with RooBTagDecay, but make RooAbsAnaConvPdf::changeModel(const RooResolutionModel& newModel) public
        change = comp.changeModel( timeResModels[resModelKey]['model']._target_() )
        if not change :
            accPdfs.append( comp.GetName() )
        else :
            raise AssertionError, 'P2VV - ERROR: multiplyByTimeAcceptance(): failed to multiply "%s" with acceptace' % comp.GetName()
    print 'P2VV - INFO: multiplyByTimeAcceptance(): multiplied the following %s with the acceptance function:'\
          % ( ( '%d components' % len(accPdfs) ) if len(accPdfs) > 1 else 'component' )
    print '             [ ' + ', '.join( name for name in accPdfs ) + ' ]'

    motherPdf['ConditionalObservables'] = motherPdf['ConditionalObservables'] | timeResModels[resModelKey].ConditionalObservables()
    motherPdf['ExternalConstraints']    = motherPdf['ExternalConstraints']    | timeResModels[resModelKey].ExternalConstraints()

    return pdf


# function to multiply a PDF by angular acceptance function
def multiplyByAngularAcceptance( pdf, self, **kwargs ) :
    anglesEffType  = getKWArg( self, kwargs, 'anglesEffType' )
    angEffMomsFile = getKWArg( self, kwargs, 'angEffMomsFile' )
    angleFuncs     = getKWArg( self, kwargs, 'angleFuncs' )
    parNamePrefix  = getKWArg( self, kwargs, 'parNamePrefix', '' )
    coefNamePF     = getKWArg( self, kwargs, 'coefNamePF', '' )

    print 'P2VV - INFO: multiplyByAngularAcceptance(): multiplying PDF "%s" with angular efficiency moments from file "%s"'\
          % ( pdf.GetName(), angEffMomsFile )

    # multiply PDF with angular efficiency
    from P2VV.Utilities.DataMoments import RealMomentsBuilder, angularMomentIndices
    moments = RealMomentsBuilder()
    moments.appendPYList( angleFuncs.angles, angularMomentIndices( anglesEffType, angleFuncs ) )
    moments.read(angEffMomsFile)
    moments.Print()
    moments.multiplyPDFWithEff( pdf, CoefName = parNamePrefix + 'effC' + '_' + coefNamePF )

    return pdf
