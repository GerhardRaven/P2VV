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
        if hasattr( values, '__iter__' ) and len(values) == 7 :
            self._parameters[name] = tuple( [ val for val in values ] )
        else :
            raise KeyError('PdfConfiguration.addParameter(): format of "values" argument should be ( value, error, error low, error high, min, max, floating? )')

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
                                               , par.getAsymErrorLo() if par.hasAsymError() else +1.
                                               , par.getAsymErrorHi() if par.hasAsymError() else -1.
                                               , par.getMin()
                                               , par.getMax()
                                               , False if par.isConstant() else True
                                              )
                             )

    def setParametersInPdf( self, pdf ) :
        for par in pdf.getVariables() :
            if not par.GetName() in self._parameters.keys() : continue
            par.setVal( self._parameters[ par.GetName() ][0] )
            par.setError( self._parameters[ par.GetName() ][1] )
            if self._parameters[ par.GetName() ][2] < 0. :
                par.setError( self._parameters[ par.GetName() ][2], self._parameters[ par.GetName() ][3] )
            par.setMin( self._parameters[ par.GetName() ][4] )
            par.setMax( self._parameters[ par.GetName() ][5] )
            par.setConstant( False if self._parameters[ par.GetName() ][6] else True )

    def readParametersFromFile( self, filePath = 'parameters', **kwargs ) :
        verb = kwargs.pop( 'Verbose', True )

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

        # check for remaining keyword arguments
        assert not kwargs, 'P2VV - ERROR: PdfConfiguration.writeParametersToFile(): superfluous keyword arguments: %s' % kwargs

        # loop over lines and read parameters
        numPars = 0
        fitStatus = ( )
        while True :
            # read next line
            line = parFile.readline()
            if not line : break

            # check for empty or comment lines
            line = line.split()
            if not line :
                continue
            elif line[0][0] == '#' :
                # check if this is the "fit status" line
                if len(line) == 12 and line[1] == 'fit' and line[2] == 'status:' :
                    fitStatus = ( int( line[5][ : -1 ] ), float( line[8][ : -1 ] ), float( line[11] ) )
                continue

            # check moment format
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
            self.addParameter( line[0], ( parVal, parErr, +1., -1., parMin, parMax, parFloat ) )
            numPars += 1

        parFile.close()

        if verb :
            print 'P2VV - INFO: PdfConfiguration.readParametersFromFile(): %d parameter%s read from file \"%s\"'\
                  % ( numPars, '' if numPars == 1 else 's', filePath )

        return fitStatus

    def writeParametersToFile( self, filePath = 'parameters', **kwargs ) :
        # get arguments
        fileFormat = kwargs.pop( 'Format', '' )
        assert fileFormat in [ '', 'common' ]\
               , 'P2VV - ERROR: PdfConfiguration.writeParametersToFile(): unknown file format: "%s"' % fileFormat

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

        # get fit status
        fitStat = kwargs.pop( 'FitStatus', ( -1, 0., 0. ) )

        # check for remaining keyword arguments
        assert not kwargs, 'P2VV - ERROR: PdfConfiguration.writeParametersToFile(): superfluous keyword arguments: %s' % kwargs

        # write parameters to content string
        cont = ''
        if fileFormat != 'common' :
            cont += '# %s: parameters\n' % fileName\
                  + '# name requirement: \'{0}\'\n'.format( names if names else '' )\
                  + '# floating:         \'{0}\'\n'.format( 'True' if floating == True else ( 'False' if floating == False else '' ) )\
                  + '# fit status:       status = {0:d}; NLL = {1:+.16g}; EDM = {2:.3g}\n'.format( fitStat[0], fitStat[1], fitStat[2] )\
                  + '#\n'\
                  + '# ' + '-' * (79 + maxLenName) + '\n'\
                  + ( '# {0:<%s}   {1:<14}   {2:<13}   {3:<14}   {4:<14}   {5:<}\n' % maxLenName )\
                      .format( 'parameter', 'value', 'error', 'min', 'max', 'floating?' )\
                  + '# ' + '-' * (79 + maxLenName) + '\n'

        numPars = 0
        from P2VV.Imports import commonParNames
        for parName in sorted( self._parameters.keys() ) :
            if fileFormat == 'common' and parName not in commonParNames : continue
            if nameExpr and not nameExpr.match(parName) : continue

            parVals = self._parameters[parName]
            if ( floating == True and not parVals[6] ) or ( floating == False and parVals[6] ) : continue

            if fileFormat != 'common' :
                cont += ( '  {0:<%s}   {1:<+14.8g}   {2:<13.8g}   {3:<+14.8g}   {4:<+14.8g}   {5:<}\n' % maxLenName )\
                          .format( parName, parVals[0], parVals[1], parVals[4], parVals[5], 'True' if parVals[6] else 'False' )
            else :
                cont += '{0:s}\t{1:+.8g}\t{2:.8g}\t{3:+.8g}\t{4:+.8g}\n'\
                          .format( commonParNames[parName], parVals[0], parVals[1], parVals[2], parVals[3] )
            numPars += 1

        if fileFormat != 'common' :
            cont += '# ' + '-' * (79 + maxLenName) + '\n'

        # write content to file
        parFile.write(cont)
        parFile.close()

        print 'P2VV - INFO: PdfConfiguration.writeParametersToFile(): %d parameter%s written to file \"%s\"'\
                % ( numPars, '' if numPars == 1 else 's', filePath )


# B_s^0 -> J/psi phi analysis configuration
class Bs2Jpsiphi_PdfConfiguration( PdfConfiguration ) :
    def __init__( self, **kwargs ) :
        from math import pi

        # get keyword arguments
        sFit = kwargs.pop( 'sFit', True )

        # job parameters
        self['sFit']          = sFit            # fit only signal?
        self['blind']         = { }             # { 'phiCP' : ( 'UnblindUniform', 'myString', 0.2 ) }
        self['numEvents']     = 60000           # sum of event yields
        self['sigFrac']       = 0.5             # fraction of signal events
        self['parNamePrefix'] = ''              # prefix for parameter names

        self['obsDict'] = dict(  runPeriod = ( 'runPeriod',            'run period',              { }                                   )
                               , mass      = ( 'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.   )
                               #, mumuMass  = ( 'mdau1',                'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.   )
                               , KKMass    = ( 'mdau2',                'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,   990.,   1050.   )
                               , KKMassCat = ( 'KKMassCat',            'KK mass category',        { }                                   )
                               , time      = ( 'time',                 'Decay time',              'ps',        1.5,    0.3,     14.     )
                               , timeRes   = ( 'sigmat',               '#sigma(t)',               'ps',        0.01,   0.0001,  0.12    )
                               , cpsi      = ( 'helcosthetaK',         'cos(#theta_{K})',         '',          0.,    -1.,     +1.      )
                               , ctheta    = ( 'helcosthetaL',         'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.      )
                               , phi       = ( 'helphi',               '#phi_{h}',                'rad',       0.,    -pi,     +pi      )
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
        self['fitOptions'] = dict( NumCPU = 2, Optimize = 2, Timer = True, Minimizer = 'Minuit2' )

        # PDF parameters
        self['numTimeBins']        = 30               # number of bins for the decay time observable
        self['numTimeResBins']     = 40               # number of bins for the decay-time resolution observable (used for caching)
        self['timeResType']        = ''               # '' / 'event' / 'eventNoMean' / 'eventConstMean' / '3Gauss' / 'event3fb'
        self['constrainTResScale'] = ''               # '' / 'constrain' / 'fixed'
        self['timeEffType']        = 'paper2012'      # 'HLT1Unbiased' / 'HLT1ExclBiased' / 'paper2012' / 'fit'
        self['timeEffConstraintType'] = 'multinomial' # 'multinomial' / 'average'
        self['constrainDeltaM']    = ''               # '' / 'constrain' / 'fixed'

        self['timeEffHistFiles']  = { }
        self['timeEffParameters'] = { }

        self['transAngles']     = False        # use transversity angles?
        self['anglesEffType']   = 'weights'    # '' / 'weights' / 'basis012' / 'basis012Plus' / 'basis012Thetal' / 'basis0123' / 'basis01234' / 'basisSig3' / 'basisSig4'
        self['angEffMomsFiles'] = ''
        self['angularRanges']   = dict( cpsi = [ ], ctheta = [ ], phi = [ ] )

        self['sigTaggingType']   = 'tagUntag'    # 'histPdf' / 'tagUntag' / 'tagCats' / 'tagUntagRelative' / 'tagCatsRelative'
        self['tagPdfType']       = ''
        self['SSTagging']        = True           # use same-side Kaon tagging?
        self['condTagging']      = True           # make tagging categories and B/Bbar tags conditional observables?
        self['contEstWTag']      = True           # use a continuous estimated wrong-tag probability instead of tagging categories?
        self['constrainTagging'] = ''             # '' / 'constrain' / 'fixed'
        self['tagCatsOS']        = [ ]            # [ ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.669, 0.0 ), ( 'Tagged',   1, 0.4999999, 0.392, 0.392, 0.0, 0.331, 0.0 ) ]
        self['tagCatsSS']        = [ ]            # [ ( 'Untagged', 0, 0.5000001, 0.5,   0.5,   0.0, 0.896, 0.0 ), ('Tagged',    1, 0.4999999, 0.359, 0.359, 0.0, 0.104, 0.0 ) ]

        self['amplitudeParam'] = 'phasesSWaveFrac'    # 'phases' / 'phasesSWaveFrac' / 'ReIm' / 'bank'
        self['ASParam']        = 'deltaPerp'          # 'delta0' / 'deltaPerp' / 'ReIm' / 'Mag2ReIm' / 'Mag2ReImPerp'
        self['AparParam']      = 'phase'              # 'phase' / 'ReIm' / 'Mag2ReIm' / 'cos' / 'real'
        self['ambiguityPars']  = False                # set parameters to values of the second minimum?

        self['paramKKMass']     = 'simultaneous'       # '' / 'parameters' / 'simultaneous'
        self['KKMassBinBounds'] = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ]    # KK-mass bin boundaries
        self['CSPValues']       = [ 0.9178, 0.9022, 0.8619, 0.8875, 0.9360, 0.9641 ]  # S-P coupling factors for KK-mass bins

        if not sFit :
            self['bkgAnglePdfType'] = 'hybrid'
            self['numAngleBins']    = ( 10, 24, 5 )

        self['lambdaCPParam'] = 'lambPhi'    # 'ReIm' / 'lambSqPhi' / 'lambPhi' / 'lambPhi_CPVDecay' / 'lambPhiRel_CPVDecay'

        self['splitParams'] = dict( KKMassCat = [ 'f_S', 'ASOddPhase' ] )

        self['externalConstr'] = { }

        # check for remaining keyword arguments
        assert not kwargs, 'P2VV - ERROR: Bs2Jpsiphi_PdfConfiguration: superfluous keyword arguments found: %s' % kwargs

        # initialize PdfConfiguration object
        PdfConfiguration.__init__( self )

    def __setitem__( self, key, val ) :
        if key == 'sFit' :
            if not val : self.addBkgParams()
            else       : self.rmBkgParams()

        return dict.__setitem__( self, key, val )

    def addBkgParams(self) :
        if not 'bkgAnglePdfType' in self : self['bkgAnglePdfType'] = 'hybrid'
        if not 'numAngleBins'    in self : self['numAngleBins']    = ( 10, 24, 5 )

    def rmBkgParams(self) :
        for key in [ 'bkgAnglePdfType', 'numAngleBins' ] : self.pop( key, None )


# B_s^0 -> J/psi phi analysis of 1 fb^-1 2011 data (paper) configuration
class Bs2Jpsiphi_2011Analysis( Bs2Jpsiphi_PdfConfiguration ) :
    def __init__( self, **kwargs ) :
        # initialize PDF configuration
        Bs2Jpsiphi_PdfConfiguration.__init__( self, **kwargs )

        # do 1/fb-analysis-specific configuration
        self['numEvents']  = 54755
        self['sigFrac']    = 0.51
        self['runPeriods'] = [ ]

        self['timeResType']   = 'eventNoMean'
        self['constrainBeta'] = 'noBeta'

        self['obsDict']['wTagOS'] = ( 'tagomega_os', 'OS est. wrong-tag prob.', '', 0.25, 0., 0.50001 )
        self['obsDict']['wTagSS'] = ( 'tagomega_ss', 'SS est. wrong-tag prob.', '', 0.25, 0., 0.50001 )

        self['timeEffHistFiles'] = dict(  file      = 'data/Bs_HltPropertimeAcceptance_Data-20120816.root'
                                        , hlt1UB    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                                        , hlt1ExclB = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                                       )
        self['angEffMomsFiles'] = 'data/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

        self['KKMassBinBounds'] = [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ]
        self['CSPValues']       = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

        self['externalConstr'] = dict(  dM             = (  17.63,  0.11   )
                                      , wTagP0OS       = (  0.392,  0.008  )
                                      , wTagP1OS       = (  1.000,  0.023  )
                                      , wTagDelP0OS    = (  0.0110, 0.0034 )
                                      , wTagDelP1OS    = (  0.,     0.     )
                                      , wTagP0SS       = (  0.350,  0.017  )
                                      , wTagP1SS       = (  1.00,   0.16   )
                                      , wTagDelP0SS    = ( -0.019,  0.005  )
                                      , wTagDelP1SS    = (  0.,     0.     )
                                      , timeResSigmaSF = (  1.45,   0.06   )
                                     )

        from P2VV.Imports import extConstraintValues
        extConstraintValues.setVal( 'DM',      (  17.63,  0.11   ) )
        extConstraintValues.setVal( 'P0OS',    (  0.392,  0.008, 0.392 ) )
        extConstraintValues.setVal( 'P1OS',    (  1.000,  0.023  ) )
        extConstraintValues.setVal( 'DelP0OS', (  0.0110, 0.0034 ) )
        extConstraintValues.setVal( 'DelP1OS', (  0.,     0.     ) )
        extConstraintValues.setVal( 'P0SS',    (  0.350,  0.017, 0.350 ) )
        extConstraintValues.setVal( 'P1SS',    (  1.00,   0.16   ) )
        extConstraintValues.setVal( 'DelP0SS', ( -0.019,  0.005  ) )
        extConstraintValues.setVal( 'DelP1SS', (  0.,     0.     ) )

        # check for remaining keyword arguments
        assert not kwargs, 'P2VV - ERROR: Bs2Jpsiphi_2011Analysis: superfluous keyword arguments found: %s' % kwargs


# B_s^0 -> J/psi phi analysis of 3 fb^-1 run-I data configuration
class Bs2Jpsiphi_RunIAnalysis( Bs2Jpsiphi_PdfConfiguration ) :
    def __init__( self, **kwargs ) :
        # get Run-I-analysis-specific keyword arguments
        runPeriods = kwargs.pop( 'RunPeriods', '3fb' )
        assert runPeriods in [ '2011', '2012', '3fb' ]\
               , 'P2VV - ERROR: Bs2Jpsiphi_RunIAnalysis: "RunPeriod" can be "2011", "2012" or "3fb" (found "%s")' % runPeriod

        # initialize PDF configuration
        Bs2Jpsiphi_PdfConfiguration.__init__( self, **kwargs )

        # do Run-I-analysis-specific configuration
        addStr = '2011' if runPeriods == '2011' else 'NewData' if runPeriods == '2012' else 'Combination'
        self['blind'] = {  'phiCP'       : ( 'UnblindUniform', 'BsPhis%s' % addStr, 0.2  )
                         , 'phiCP_A0'    : ( 'UnblindUniform', 'BsPhiszero%s' % addStr, 0.3  )
                         , 'phiCP_Apar'  : ( 'UnblindUniform', 'BsPhispara%s' % addStr, 0.3  )
                         , 'phiCP_Aperp' : ( 'UnblindUniform', 'BsPhisperp%s' % addStr, 0.3  )
                         , 'phiCP_AS'    : ( 'UnblindUniform', 'BsPhisS%s' % addStr, 0.3  )
                         , 'dGamma'      : ( 'UnblindUniform', 'BsDGs%s'  % addStr, 0.02 )
                        }
        self['numEvents']  = 220000
        self['sigFrac']    = 0.43
        self['runPeriods'] = [ 2011 ] if runPeriods == '2011' else [ 2012 ] if runPeriods == '2012' else [ 2011, 2012 ]

        self['timeResType']   = 'event3fb'
        self['constrainBeta'] = ''

        self['obsDict']['wTagOS'] = ( 'tagomega_os_cb', 'OS est. wrong-tag prob.', '', 0.25, 0., 0.50001 )
        self['obsDict']['wTagSS'] = ( 'tagomega_ss_nn', 'SS est. wrong-tag prob.', '', 0.25, 0., 0.50001 )

        timeEff2011 = dict(  file      = 'data/Bs_HltPropertimeAcceptance_Data_2011_40bins.root'
                           , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
                           , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                          )
        timeEff2012 = dict(  file      = 'data/Bs_HltPropertimeAcceptance_Data_2012_40bins.root'
                           , hlt1UB    = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
                           , hlt1ExclB = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                          )
        if runPeriods == '2011' :
            self['timeEffHistFiles'] = timeEff2011
        elif runPeriods == '2012' :
            self['timeEffHistFiles'] = timeEff2012
        else :
            timeEffHistFiles = SimulCatSettings('timeEffHistFiles')
            timeEffHistFiles.addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ], timeEff2011 )
            timeEffHistFiles.addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ], timeEff2012 )
            self['timeEffHistFiles'] = timeEffHistFiles

        self['angEffMomsFiles'] = 'data/Sim08_20112012_hel_UB_UT_trueTime_BkgCat050_KK30_Phys_moms_norm'

        self['KKMassBinBounds'] = [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ]
        self['CSPValues']       = [ 0.9178, 0.9022, 0.8619, 0.8875, 0.9360, 0.9641 ]

        self['externalConstr'] = dict(  wTagP0OS       = (  0.3853, 0.0044 )
                                      , wTagP1OS       = (  1.000,  0.035  )
                                      , wTagDelP0OS    = (  0.0140, 0.     )
                                      , wTagDelP1OS    = (  0.066,  0.     )
                                      , wTagP0SS       = (  0.450,  0.005  )
                                      , wTagP1SS       = (  1.00,   0.09   )
                                      , wTagDelP0SS    = ( -0.016,  0.     )
                                      , wTagDelP1SS    = (  0.007,  0.     )
                                      , dM             = (  17.768, 0.024  )
                                     )
        from collections import defaultdict
        splitConstr = defaultdict(dict)
        splitConstr['betaTimeEff']['2011']      = ( -0.0083,  0.004 )
        if runPeriods == '3fb':
            splitConstr['betaTimeEff']['2012'] = ( -0.0083, None )
        else:
            splitConstr['betaTimeEff']['2012'] = ( 0., 0. )
        splitConstr['sf_placeholder']['2011']   = (  0.0350,  0. )
        splitConstr['sf_placeholder']['2012']   = (  0.0349,  0. )
        splitConstr['timeResMu']['2011']        = ( -0.00259, 0. )
        splitConstr['timeResMu']['2012']        = ( -0.00333, 0. )
        splitConstr['timeResFrac2']['2011']     = (  0.242,   0. )
        splitConstr['timeResFrac2']['2012']     = (  0.239,   0. )
        splitConstr['sf_mean_slope']['2011']   = (  1.4273,  0. )
        splitConstr['sf_mean_slope']['2012']   = (  1.4887,  0. )
        splitConstr['sf_mean_quad']['2011']    = ( -1.93,    0. )
        splitConstr['sf_mean_quad']['2012']    = ( -3.88,    0. )
        splitConstr['sf_sigma_slope']['2011']  = (  0.3857,  0. )
        splitConstr['sf_sigma_slope']['2012']  = (  0.4143,  0. )
        splitConstr['sf_sigma_quad']['2011']   = ( -0.63,    0. )
        splitConstr['sf_sigma_quad']['2012']   = ( -2.80,    0. )
        if runPeriods in [ '2011', '2012' ] :
            self['externalConstr']['betaTimeEff']     = splitConstr['betaTimeEff']     [runPeriods]
            self['externalConstr']['sf_placeholder']  = splitConstr['sf_placeholder'][runPeriods]
            self['externalConstr']['timeResMu']       = splitConstr['timeResMu']       [runPeriods]
            self['externalConstr']['timeResFrac2']    = splitConstr['timeResFrac2']    [runPeriods]
            self['externalConstr']['sf_mean_slope']  = splitConstr['sf_mean_slope']  [runPeriods]
            self['externalConstr']['sf_mean_quad']   = splitConstr['sf_mean_quad']   [runPeriods]
            self['externalConstr']['sf_sigma_slope'] = splitConstr['sf_sigma_slope'] [runPeriods]
            self['externalConstr']['sf_sigma_quad']  = splitConstr['sf_sigma_quad']  [runPeriods]
        else :
            self['splitParams']['runPeriod'] = [ ]
            for par, vals in splitConstr.iteritems() :
                self['splitParams']['runPeriod'].append(par)
                constr = SimulCatSettings( '%sConstr' % par )
                constr.addSettings( [ 'runPeriod' ], [ [ 'p2011' ] ], vals['2011'] )
                constr.addSettings( [ 'runPeriod' ], [ [ 'p2012' ] ], vals['2012'] )
                self['externalConstr'][par] = constr

        from P2VV.Imports import extConstraintValues
        extConstraintValues.setVal( 'DM',      ( 17.768,  0.024   ) )
        extConstraintValues.setVal( 'P0OS',    (  0.3853, 0.0044, 0.3853 ) )
        extConstraintValues.setVal( 'P1OS',    (  1.00,   0.035  ) )
        extConstraintValues.setVal( 'DelP0OS', (  0.0140, 0.0012 ) )
        extConstraintValues.setVal( 'DelP1OS', (  0.066,  0.012  ) )
        extConstraintValues.setVal( 'P0SS',    (  0.450,  0.005, 0.450 ) )
        extConstraintValues.setVal( 'P1SS',    (  1.00,   0.09   ) )
        extConstraintValues.setVal( 'DelP0SS', ( -0.016,  0.002  ) )
        extConstraintValues.setVal( 'DelP1SS', (  0.007,  0.019  ) )

        # check for remaining keyword arguments
        assert not kwargs, 'P2VV - ERROR: Bs2Jpsiphi_RunIAnalysis: superfluous keyword arguments found: %s' % kwargs


class SimulCatSettings(list) :
    def __init__( self, Name ) :
        self._name = Name
        self._cats = set()
        self._default = None

    def name() :
        return self._name

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

    def getAllSettings( self ):
        return [e[1] for e in self]
            
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
                if ws.var( obsDict[name][0] ) :
                    self['observables'][name] = RealVar( obsDict[name][0] )
                elif ws.cat( obsDict[name][0] ) :
                    self['observables'][name] = Category( obsDict[name][0] )
                else :
                    raise RuntimeError( 'P2VV - ERROR: PdfBuilder._createObservables(): variable "%s" ("%s") not in workspace'\
                                        % ( obsDict[name][0], name ) )
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
        for par in [ 'sFit', 'KKMassBinBounds', 'obsDict', 'CSPValues', 'condTagging', 'contEstWTag', 'SSTagging', 'transAngles'
                    , 'numEvents', 'sigFrac', 'paramKKMass', 'amplitudeParam', 'ASParam', 'signalData', 'fitOptions', 'parNamePrefix'
                    , 'tagPdfType', 'timeEffType', 'timeEffHistFiles', 'timeEffParameters', 'anglesEffType', 'angEffMomsFiles'
                    , 'readFromWS', 'splitParams', 'externalConstr', 'runPeriods', 'timeEffConstraintType' ] :
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

        if not self['sFit'] :
            # build background PDFs
            self['pdfComps'].append( Component( namePF + 'bkg', [ ]
                                    , Yield = ( self['numEvents'] * ( 1. - self['sigFrac'] ), 0., self['numEvents'] ) ) )
            self['pdfComps'][1] += buildBs2JpsiphiCombBkgTimePdf(self)
            self['pdfComps'][1] += buildBs2JpsiphiCombBkgAnglesPdf(self)
            if self['tagPdfType'] or not self['condTagging'] : self['pdfComps'][1] += buildTaggingPdf( self, data = self['signalData'] )

        # build full PDF
        from P2VV.RooFitWrappers import buildPdf
        print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: requesting %s PDF for observables [%s]'\
              % ( 'signal' if self['sFit'] else 'signal + background', ', '.join( str(obs) for obs in self['obsSetP2VV'] ) )
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
        sFit             = getKWArg( self, kwargs, 'sFit' )
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

        if self['runPeriods'] :
            # check run-period category
            assert len( self['obsDict']['runPeriod'][2] ) == observables['runPeriod'].numTypes()\
                   , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: number of run periods specified (%d) does not match number of types of run-period category (%d)'\
                   % ( len( self['obsDict']['runPeriod'][2] ), observables['runPeriod'].numTypes() )
            for name, index in self['obsDict']['runPeriod'][2].iteritems() :
                assert observables['runPeriod'].isValidLabel(name) and observables['runPeriod'].isValidIndex(index)\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: run period "%s" ("%d") not found in run-period category' % ( name, index )

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
        if not sFit :
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
            pdfVars = set(v.GetName() for v in fullPdf.getVariables())
            
            # get splitting categories and parameters
            splitParsDict = { }
            for cat, params in splitParams.iteritems() :
                assert ws.cat(cat), 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: category "%s" not in workspace' % cat
                for par in params :
                    par = namePF + par
                    if par not in pdfVars:
                        continue
                    assert ws.var(par), 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no variable "%s" in workspace' % par
                    if ws[par] not in splitParsDict :
                        splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
                    else :
                        splitParsDict[ ws[par] ].add( ws[cat] )
            self['splitParsDict'] = splitParsDict

            # create list of split parameters and categories
            from P2VV.Utilities.General import createSplitParsList
            splitPars = createSplitParsList(splitParsDict)

            # build simultaneous PDF
            print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: building simultaneous PDF "%s":' % ( fullPdf.GetName() + '_simul' )
            print 13 * ' ' + 'split parameters:'
            for it, pars in enumerate(splitPars) :
                print 13 * ' ' + '%2d: pars: [ %s ]' % ( it, ', '.join( par.GetName() for par in pars[0] ) )
                print 13 * ' ' + '    cats: [ %s ]' % ', '.join( cat.GetName() for cat in pars[1] )
            from P2VV.RooFitWrappers import SimultaneousPdf
            self['simulPdf'] = SimultaneousPdf(  fullPdf.GetName() + '_simul'
                                               , MasterPdf       = fullPdf
                                               , SplitCategories = [ pars[1] for pars in splitPars ]
                                               , SplitParameters = [ pars[0] for pars in splitPars ] )

            # set time resolution models
            from ROOT import RooArgSet, RooArgList
            splitCatPars = RooArgSet()
            self['simulPdf'].treeNodeServerList(splitCatPars)
            splitCat  = self['simulPdf'].indexCat()
            inputCats = RooArgList(splitCat) if splitCat.isFundamental() else splitCat.inputCatList()
            prototype = timeResModelsOrig['prototype']['model']
            from P2VV.RooFitWrappers import AddModel
            ## NOTE, this only works for either a single resolution model or an AddModel which does
            ## not contain further add models
            if isinstance(prototype, AddModel):
                origResParams = {}
                for model in prototype.models():
                    assert(not isinstance(model, AddModel))
                    origResParams[model.GetName()] = model['Parameters']
                origResParams[prototype.GetName()] = prototype.fractions()
            else:
                origResParams = {prototype.GetName() : prototype['Parameters']}

            from P2VV.Utilities.Splitting import replacement_model
            # Figure out which way the time resolution model is to be split
            pfLen = len(namePF)
            resParams = set( var.GetName()[ pfLen : ] for var in timeResModelsOrig['prototype']['model'].getVariables() )
            splitResCats = [c for c, pars in splitParams.iteritems() if set(pars).intersection(resParams)]
            splitResCats = [inputCats.find(c) for c in splitResCats]
            replacements = { ( ) : timeResModelsOrig['prototype'] }
            for splitCatState in splitCat:
                splitCat.setIndex( splitCatState.getVal() )
                k = tuple(c.getLabel() for c in splitResCats)
                if k not in replacements:
                    catPdf = self['simulPdf'].getPdf( splitCatState.GetName() )
                    resModelCount = 0
                    from ROOT import RooBTagDecay
                    for comp in filter( lambda x : isinstance( x, RooBTagDecay ), catPdf.getComponents() ) :
                        # TODO: don't do this with RooBTagDecay, but move RooBTagDecay::resModel() to RooAbsAnaConvPdf
                        assert resModelCount < 1, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: multiple resolution models found for simultaneous category "%s"' % splitCatState.GetName()
                        replacements[k] = replacement_model(prototype, comp.resolutionModel(), inputCats, self['splitParsDict'], splitCatPars, origResParams)
                        resModelCount += 1
                    assert resModelCount > 0, ('P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: no resolution'
                                               + ' model found for simultaneous category "%s"' % splitCatState.GetName())
                timeResModelsOrig[splitCatState.GetName()] = replacements[k]
                    
            if paramKKMass == 'simultaneous' and ASParam != 'Mag2ReIm' :
                # set values for split S-P coupling factors
                if ASParam != 'Mag2ReIm' :
                    print 'P2VV - INFO: Bs2Jpsiphi_PdfBuilder: using S-P-wave coupling factors:',
                    for iter, fac in enumerate( CSPValues ) :
                        print '%d: %.4f%s' % ( iter, fac, '' if iter == len( CSPValues ) - 1 else ',' ),
                    print

                from P2VV.Utilities.General import getSplitPar
                for splitCatState in splitCat:
                    splitCat.setIndex( splitCatState.getVal() )
                    C_SP = getSplitPar( namePF + 'C_SP', observables['KKMassCat'].getLabel(), splitCatPars )
                    C_SP.setVal( CSPValues[ observables['KKMassCat'].getIndex() ] )
                    C_SP.setConstant(True)

        else :
            self['simulPdf'] = None

        return self['simulPdf'] if self['simulPdf'] else self['fullPdf']


    def _multiplyByTimeAcceptance( self, **kwargs ) :
        timeEffType       = getKWArg( self, kwargs, 'timeEffType' )
        timeEffHistFiles  = getKWArg( self, kwargs, 'timeEffHistFiles' )
        pdf               = getKWArg( self, kwargs, 'pdf' )
        simulPdf          = getKWArg( self, kwargs, 'simulPdf' )
        signalData        = getKWArg( self, kwargs, 'signalData' )
        observables       = getKWArg( self, kwargs, 'observables' )
        timeResModels     = getKWArg( self, kwargs, 'timeResModels' )
        timeResModelsOrig = getKWArg( self, kwargs, 'timeResModelsOrig' )

        if not simulPdf :
            # original PDF is not simultaneous: multiply with acceptance from outside PDF
            if type(timeEffHistFiles) == SimulCatSettings :
                timeEffHistFiles = timeEffHistFiles.default()
            multiplyByTimeAcceptance( pdf, self, data = signalData, histFile = timeEffHistFiles['file']
                                     , histUBName = timeEffHistFiles['hlt1UB'], histExclBName = timeEffHistFiles['hlt1ExclB'] )
        else :
            # original PDF is simultaneous: multiply with acceptance separately for all categories
            from ROOT import RooArgList
            splitCat = simulPdf.indexCat()
            inputCats = RooArgList(splitCat) if splitCat.isFundamental() else splitCat.inputCatList()
            splitAccCats = [ ]
            if type(timeEffHistFiles) == SimulCatSettings :
                splitAccCats = [ inputCats.find(cat) for cat in timeEffHistFiles.categories() ]
                assert all( cat for cat in splitAccCats )\
                       , 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: PDF was not split for all categories in "%s"' % timeEffHistFiles.name()
            singleHLT1Eff = timeEffType.startswith('paper2012')\
                            and any( cat.GetName() == observables['hlt1ExclB'].GetName() for cat in inputCats )

            # loop over simultaneous categories
            resModels = { }
            for splitCatState in splitCat:
                # get category state and corresponding acceptance parameters
                splitCat.setIndex( splitCatState.getVal() )
                catPdf = simulPdf.getPdf( splitCatState.GetName() )
                if type(timeEffHistFiles) == SimulCatSettings :
                    effFile = timeEffHistFiles.getSettings( [ ( cat.GetName(), cat.getLabel() ) for cat in splitAccCats ] )
                else :
                    effFile = timeEffHistFiles

                cNamePF = '_'.join( '%s_%s' % ( cat.GetName(), cat.getLabel() ) for cat in splitAccCats )
                if singleHLT1Eff :
                    cNamePF += '_%s_%s' % ( observables['hlt1ExclB'].GetName(), observables['hlt1ExclB'].getLabel() )
                    if observables['hlt1ExclB'].getIndex() == 0 :
                        timeEffType = 'HLT1Unbiased'
                    elif observables['hlt1ExclB'].getIndex() == 1 :
                        timeEffType = 'HLT1ExclBiased'
                    else :
                        raise AssertionError, 'P2VV - ERROR: Bs2Jpsiphi_PdfBuilder: unknown HLT1 exclusively biased state: %s (%d)'\
                                              % ( observables['hlt1ExclB'].getLabel(), observables['hlt1ExclB'].getIndex() )

                if cNamePF in resModels and timeResModelsOrig[ splitCatState.GetName() ] == timeResModelsOrig[ resModels[cNamePF] ] :
                    # use resolution model with acceptance function from a previous category
                    prevResModelKey = resModels[cNamePF]
                else :
                    # create new resolution model with acceptance function
                    resModels[cNamePF] = splitCatState.GetName()
                    prevResModelKey = ''

                # multiply PDF for this category state with acceptance
                multiplyByTimeAcceptance( catPdf, self, data = signalData, timeEffType = timeEffType
                                         , histFile = effFile['file'], histUBName = effFile['hlt1UB']
                                         , histExclBName = effFile['hlt1ExclB'], coefNamePF = cNamePF
                                         , motherPdf = simulPdf, resModelKey = splitCatState.GetName()
                                         , timeResModel = timeResModels[prevResModelKey] if prevResModelKey else None )


    def _multiplyByAngularAcceptance( self, **kwargs ) :
        angEffMomsFiles = getKWArg( self, kwargs, 'angEffMomsFiles' )
        pdf             = getKWArg( self, kwargs, 'pdf' )
        simulPdf        = getKWArg( self, kwargs, 'simulPdf' )

        if type(angEffMomsFiles) == str :
            multiplyByAngularAcceptance( pdf, self, angEffMomsFile = angEffMomsFiles )
        else :
            splitCat      = simulPdf.indexCat()
            inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()
            for splitCatState in splitCat:
                splitCat.setIndex( splitCatState.getVal() )
                effFile = angEffMomsFiles.getSettings( [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats ] )
                catPdf = simulPdf.getPdf( splitCatState.GetName() )
                cNamePF = ( splitCatState.GetName() ).replace( '{', '' ).replace( '}', '' ).replace( ';', '_' )
                multiplyByAngularAcceptance( catPdf, self, angEffMomsFile = effFile, coefNamePF = cNamePF )

    def _createExternalConstraints( self, **kwargs ) :
        externalConstr        = getKWArg( self, kwargs, 'externalConstr' )
        splitParsDict         = getKWArg( self, kwargs, 'splitParsDict' )
        pdf                   = getKWArg( self, kwargs, 'pdf' )
        simulPdf              = getKWArg( self, kwargs, 'simulPdf' )
        timeEffConstraintType = getKWArg( self, kwargs, 'timeEffConstraintType' )
        data                  = getKWArg( self, kwargs, 'signalData' )
        observables           = getKWArg( self, kwargs, 'observables' )

        from P2VV.Parameterizations.GeneralUtils import getParNamePrefix
        namePF = getParNamePrefix(True)

        constraints = set()
        def buildConstraint( par, constrVals ) :
            par.setVal( constrVals[0] )
            try :
                constrErr = float(constrVals[1])
            except :
                constrErr = None

            if constrErr == None :
                par.setConstant(False)
            elif constrErr <= 0. :
                par.setConstant(True)
            else :
                par.setConstant(False)
                from P2VV.RooFitWrappers import Pdf, ConstVar
                from ROOT import RooGaussian as Gaussian
                constrName = par.GetName().replace( '{', '' ).replace( '}', '' ).replace( ';', '_' )
                constraints.add( Pdf(  Name = constrName + '_constraint', Type = Gaussian
                                     , Parameters = [  par
                                                     , ConstVar( Name = constrName + '_mean',  Value = constrVals[0] )
                                                     , ConstVar( Name = constrName + '_sigma', Value = constrErr     )
                                                    ]
                                    )
                               )
        if timeEffConstraintType == 'multinomial':
            constraints |= self.__multinomial_eff_constraints(externalConstr, simulPdf,
                                                              data, observables)
        elif timeEffConstraintType == 'average':
            constraints |= self.__av_eff_constraints(externalConstr)

        ws = pdf.ws()
        pdfVars = pdf.getVariables()

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

                catsStrings = set()
                splitCat = simulPdf.indexCat()
                from ROOT import RooArgList
                inputCats = RooArgList(splitCat) if splitCat.isFundamental() else splitCat.inputCatList()
                for splitCatState in splitCat:
                    splitCat.setIndex( splitCatState.getVal() )
                    catLabels = [ ( cat.GetName(), cat.getLabel() ) for cat in inputCats if cat in splitCats ]
                    catsStr = ';'.join( lab[1] for lab in catLabels )
                    if len(catLabels) > 1 : catsStr = '{' + catsStr + '}'
                    if catsStr in catsStrings : continue

                    catsStrings.add(catsStr)
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

    def __multinomial_eff_constraints(self, externalConstr, simulPdf, data, observables):
        # Time res model constraints include those form the acceptance
        constraints = set()
        acc_settings = externalConstr.pop('acceptance', None)
        if acc_settings and simulPdf:
            for (key, model) in self['timeResModels'].iteritems():
                if key == 'prototype' or not hasattr(model, 'build_multinomial_constraints'):
                    continue
                constraints |= model.build_multinomial_constraints(data, observables)
        elif acc_settings:
            acceptance = timeResModels['prototype']
            if hasattr(acceptance, 'build_multinomial_constraints'):
                constraints |= acceptance.build_multinomial_constraints(data, observables)
        return constraints
    
    def __av_eff_constraints(self, externalConstr, simulPdf):
        # Time res model constraints include those form the acceptance
        from P2VV.Parameterizations.TimePDFs import Single_Exponent_Time
        from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution
        truth_res = Truth_TimeResolution(time = self['obsSetP2VV'][0])
        simple_time = Single_Exponent_Time(Name = 'constraint_time', time = self['obsSetP2VV'][0],
                                           resolutionModel = truth_res.model(),
                                           tau = self['lifetimeParams']['MeanLifetime'])
        simple_time_pdf = simple_time.pdf()
        acc_settings = externalConstr.pop('acceptance', None)
        if simulPdf:
            splitCat = simulPdf.indexCat()
            from ROOT import RooArgList
            inputCats = RooArgList(splitCat) if splitCat.isFundamental() else splitCat.inputCatList()

        constraints = set()
            
        if acc_settings and simulPdf:
            splitAccCats = [inputCats.find(c if type(c) == str else c.GetName()) for c in acc_settings.categories()]

            ## Create a simultaneous pdf for the constraints.
            from P2VV.Utilities.General import createSplitParsList

            splitSet = set(k.GetName() for k in splitParsDict.iterkeys())
            accSplitDict = dict((var, set(splitAccCats)) for var in simple_time_pdf.getVariables() if var.GetName() in splitSet)
            accSplitPars = createSplitParsList(accSplitDict)

            # build simultaneous PDF
            from P2VV.RooFitWrappers import SimultaneousPdf
            simul_time_pdf = SimultaneousPdf(simple_time_pdf.GetName() + '_simul',
                                             MasterPdf       = simple_time_pdf,
                                             SplitCategories = [pars[1] for pars in accSplitPars],
                                             SplitParameters = [pars[0] for pars in accSplitPars])

            from ROOT import RooRealVar
            simulVars = simul_time_pdf.getVariables()
            for v in simulPdf.getVariables():
                if v.getAttribute('Observable'):
                    continue
                if not isinstance(v, RooRealVar):
                    continue
                sv = simulVars.find(v.GetName())
                if sv:
                    sv.setVal(v.getVal())
                    sv.setError(v.getError())
                            
            simple_split_cat = simul_time_pdf.indexCat()
            for (key, model) in self['timeResModels'].iteritems():
                if key == 'prototype':
                    continue
                if not hasattr(model, 'build_av_constraints'):
                    continue
                splitCat.setLabel(key)
                state = simple_split_cat.getLabel()
                split_pdf = simul_time_pdf.getPdf(state)

                sk = tuple((cat.GetName(), cat.getLabel()) for cat in splitAccCats)
                values = acc_settings.getSettings(sk)
                # I use a simple PDF with only the average lifetime here. I
                # think this is a valid approximation.
                constraints |= model.build_av_constraints(split_pdf, values)
        elif acc_settings:
            acceptance = timeResModels['prototype']
            if hasattr(acceptance, 'build_av_constraints'):
                constraints |= acceptance.build_av_constraints(simple_time_pdf, acc_settings)

        return constraints

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
        else :               dGammaVar['Blind'] = ( 'UnblindUniform', 'BsDGs', 0.02 )
    self['lifetimeParams'] = LifetimeParams( dGamma = dGammaVar, dMConstraint = constrainDeltaM, betaConstraint = constrainBeta )
    if ambiguityPars : self['lifetimeParams']['dGamma'].setVal( -self['lifetimeParams']['dGamma'].getVal() )

    if timeResType.startswith('event') :
        timeResArgs = dict( time = observables['time'], Cache = False if timeEffType else True )
        if 'nomean' in timeResType.lower() :
            from P2VV.RooFitWrappers import ConstVar
            from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
            timeResArgs['timeResMean']   = ConstVar( Name = namePF + 'timeResMean',   Value = 0. )
            timeResArgs['timeResMeanSF'] = ConstVar( Name = namePF + 'timeResMeanSF', Value = 1. )
            timeResArgs['timeResSFConstraint'] = constrainTResScale
            timeResArgs['timeResSigma'] = observables['timeRes']
        elif 'constmean' in timeResType.lower() :
            from P2VV.RooFitWrappers import ConstVar
            from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
            timeResArgs['timeResMean']   = dict( Value = -0.01, Error = 0.005 )
            timeResArgs['timeResMeanSF'] = ConstVar( Name = namePF + 'timeResMeanSF', Value = 1. )
            timeResArgs['timeResMeanConstraint'] = constrainTResScale
            timeResArgs['timeResSFConstraint'] = constrainTResScale
            timeResArgs['timeResSigma'] = observables['timeRes']
        elif '3fb' in timeResType.lower() :
            from P2VV.Parameterizations.TimeResolution import Multi_Gauss_TimeResolution as TimeResolution
            timeResArgs = dict( time = observables['time'], sigmat = observables['timeRes'], Cache = True, PerEventError = True
                               , Parameterise = 'RMS', TimeResSFParam = 'linear'
                               , ScaleFactors     = [ ( 2, 1.817 ), ( 1, 1.131 ) ]
                               , timeResMu        = dict( Value = 0., Constant = True ) #dict( Value = -0.00298, Constant = True )
                               , Fractions        = [ ( 2, 0.168 ) ]
                               , timeResFrac2     = dict( Value =  0.4,  Constant = True ) #dict( Value =  0.24295,   Constant = True )
                               , sf_mean_slope    = dict( Value =  1.45, Constant = True ) #dict( Value =  1.42479,   Constant = True )
                               , sf_mean_quad     = dict( Value =  0.,   Constant = True ) #dict( Value = -0.0526273, Constant = True )
                               , sf_sigma_slope   = dict( Value =  0.,   Constant = True ) #dict( Value =  0.381861,  Constant = True )
                               , sf_sigma_quad    = dict( Value =  0.,   Constant = True ) #dict( Value = -0.0147151, Constant = True )
                               , sf_placeholder   = dict( Value =  0.,   Constant = True ) #dict( Value =  0.032,     Constant = True )
                              )
        else :
            from P2VV.Parameterizations.TimeResolution import Paper2012_TimeResolution as TimeResolution
            timeResArgs['timeResMean'] = dict( Value = -4.0735e-03, Error = 1.33e-04 )
            timeResArgs['timeResMeanConstraint'] = constrainTResScale
            timeResArgs['timeResSFConstraint'] = constrainTResScale
            timeResArgs['timeResSigma'] = observables['timeRes']

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
            from P2VV.RooFitWrappers import RealVar
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
            ampNames = [ 'A0', 'Apar', 'Aperp', 'AS' ]
            phiCPPars = dict( [ ( 'phiCP_%s' % amp, dict( Name = 'phiCP_%s' % amp ) ) for amp in ampNames ] )
            for name, parDict in phiCPPars.iteritems() :
                if blind and name in blind :
                    if blind[name] : parDict['Blind'] = blind[name]
                    else           : parDict['Blind'] = ( 'UnblindUniform', 'BsPhis',  0.2 )
            self['lambdaCP'] = CPParam( AmplitudeNames = ampNames, Amplitudes = self['amplitudes'], **phiCPPars )

    else :
        if lambdaCPParam == 'lambPhi' :
            from P2VV.Parameterizations.CPVParams import LambdaAbsArg_CPParam as CPParam
        else :
            from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam as CPParam

        phiCPVar = dict( Name = 'phiCP' )
        if blind and 'phiCP' in blind :
            if blind['phiCP'] : phiCPVar['Blind'] = blind['phiCP']
            else              : phiCPVar['Blind'] = ( 'UnblindUniform', 'BsPhis',  0.2 )
        self['lambdaCP'] = CPParam( phiCP = phiCPVar )
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
                                  , Value = 1., MinMax = (  0., 2. ), Constant = True )
            avgCOddSum  = RealVar( namePF + 'avgCOddSum', Title = 'Sum of CP average odd coefficients'
                                  , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ), Constant = True )
            avgCEvenOS  = RealVar( namePF + 'avgCEvenOSTagged', Title = 'CP average even coefficients OS tagged categories'
                                  , Value = 1., MinMax = (  0., 2. ), Constant = True )
            avgCOddOS   = RealVar( namePF + 'avgCOddOSTagged' , Title = 'CP average odd coefficients OS tagged categories'
                                  , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ), Constant = True )
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
                                     , Value = 1., MinMax = (  0., 2. ), Constant = True )
                avgCOddSS  = RealVar( namePF + 'avgCOddSSTagged', Title = 'CP average odd coefficients SS tagged categories'
                                     , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ), Constant = True )
                avgCEven   = RealVar( namePF + 'avgCEvenTagged', Title = 'CP average even coefficients OS+SS tagged categories'
                                     , Value = 1., MinMax = (  0., 2. ), Constant = True )
                avgCOdd    = RealVar( namePF + 'avgCOddTagged', Title = 'CP average odd coefficients OS+SS tagged categories'
                                     , Value = asymVal, Error = asymErr, MinMax = ( -2., 2. ), Constant = True )

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
    timeResModel      = getKWArg( self, kwargs, 'timeResModel', None )
    timeResModels     = getKWArg( self, kwargs, 'timeResModels' )
    timeResModelsOrig = getKWArg( self, kwargs, 'timeResModelsOrig' )
    resModelKey       = getKWArg( self, kwargs, 'resModelKey', 'prototype' )

    print 'P2VV - INFO: multiplyByTimeAcceptance(): multiplying PDF "%s" with decay-time acceptance function' % pdf.GetName()
    if timeResModel :
        # use resolution model with acceptance function from arguments
        timeResModels[resModelKey] = timeResModel

    else :
        # build new resolution model with acceptance function
        from P2VV.Parameterizations.GeneralUtils import setParNamePrefix
        setParNamePrefix( parNamePrefix + ( '_' if parNamePrefix else '' ) + coefNamePF )

        time      = observables['time']
        hlt1ExclB = observables['hlt1ExclB']
        hlt2B     = observables['hlt2B']
        hlt2UB    = observables['hlt2UB']

        # get index categories if mother PDF is a simultaneous PDF
        if hasattr( motherPdf, 'indexCat' ) :
            indexCat = motherPdf.indexCat()
            indexCatNames = [ indexCat.GetName() ] if indexCat.isFundamental() else [ cat.GetName() for cat in indexCat.inputCatList() ]
        else :
            indexCatNames = [ ]

        # build new decay-time resolution model that includes the decay-time acceptance function
        if timeEffType == 'fit' :
            assert all( cat.GetName() not in indexCatNames for cat in [ hlt1ExclB, hlt2B, hlt2UB ] )\
                   , 'P2VV - ERROR: multiplyByTimeAcceptance(): acceptance function depends on the index category of the simultaneous mother PDF'

            hists = {  hlt1ExclB : {  'exclB'    : { 'histogram' : histExclBName }
                                    , 'notExclB' : { 'bins'      : time.getRange(), 'heights' : [0.7] }
                                   }
                     , hlt2B     : { 'B'  : { 'histogram' : histUBName } }
                     , hlt2UB    : { 'UB' : { 'histogram' : histUBName } }
                    }

            from P2VV.Parameterizations.TimeAcceptance import Paper2012_csg_TimeAcceptance as TimeAcceptance
            timeResModels[resModelKey] = TimeAcceptance(time = time
                                                        , ResolutionModel = timeResModelsOrig[resModelKey]
                                                        , Input = histFile
                                                        , Histograms = hists
                                                        , Fit = True
                                                        , **timeEffParameters
                                                        )

        elif timeEffType == 'paper2012_multi' :
            assert hlt1ExclB.GetName() not in indexCatNames\
                   , 'P2VV - ERROR: multiplyByTimeAcceptance(): acceptance function depends on the index category of the simultaneous mother PDF'

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
            assert hlt1ExclB.GetName() not in indexCatNames\
                   , 'P2VV - ERROR: multiplyByTimeAcceptance(): acceptance function depends on the index category of the simultaneous mother PDF'

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
        new_model = timeResModels[resModelKey]['model']._target_()
        change = comp.changeModel( new_model )
        if not change :
            # remove ORIGNAME:<"old model name"> attribute; this attribute gives problems with customization of the PDF later on
            for att in new_model.attributes() :
                if att.startswith( 'ORIGNAME:' + timeResModelsOrig[resModelKey]['model'].GetName().split('{')[0] ) :
                    new_model.setAttribute( att, False )
            accPdfs.append( comp.GetName() )
        else :
            raise AssertionError, 'P2VV - ERROR: multiplyByTimeAcceptance(): failed to multiply "%s" with acceptace' % comp.GetName()
    print 'P2VV - INFO: multiplyByTimeAcceptance(): multiplied the following %s with acceptance function "%s":'\
          % ( ( '%d components' % len(accPdfs) ) if len(accPdfs) > 1 else 'component', timeResModels[resModelKey].acceptance() )
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
    moments.multiplyPDFWithEff( pdf, CoefName = parNamePrefix + ( '_' if parNamePrefix else '' ) + 'effC'\
                                                + ( '_' if coefNamePF else '' ) + coefNamePF )

    return pdf



