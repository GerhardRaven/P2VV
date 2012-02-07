###########################################################################################################################################
## P2VVParameterizations.FullPDFs: Parameterizations of complete PDFs that are used in an analysis                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

Bs2Jpsiphi_Winter2012 = {}

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

        angEffMomentsFile = pdfConfig.pop('angEffMomentsFile')

        nTupleName = pdfConfig.pop('nTupleName')
        nTupleFile = pdfConfig.pop('nTupleFile')

        timeEffHistFile = pdfConfig.pop('timeEffHistFile')
        timeEffHistName = pdfConfig.pop('timeEffHistName')

        fitOpts = pdfConfig.pop('fitOptions')

        # PDF options
        components        = pdfConfig.pop('components')
        transAngles       = pdfConfig.pop('transversityAngles')
        bkgAnglePdf       = pdfConfig.pop('bkgAnglePdf')
        perEventTimeRes   = pdfConfig.pop('perEventTimeRes')
        multiplyByTimeEff = pdfConfig.pop('multiplyByTimeEff')

        tagConds = pdfConfig.pop('taggingConditionals')

        sigFrac = pdfConfig.pop('signalFraction')
        massRangeBackground = pdfConfig.pop('massRangeBackground')

        # transversity amplitudes
        pdfConfig['amplitudeParam'] = pdfConfig.pop('amplitudeParam')

        A0Mag2    = pdfConfig.pop('A0Mag2')
        AperpMag2 = pdfConfig.pop('AperpMag2')
        AparMag2  = pdfConfig.pop('AparMag2')

        A0Ph      = pdfConfig.pop('A0Phase')
        AperpPh   = pdfConfig.pop('AperpPhase')
        AparPh    = pdfConfig.pop('AparPhase')

        ASMag2 = pdfConfig.pop('ASMag2')
        ASPh   = pdfConfig.pop('ASPhase')
        if 'f_S' in pdfConfig                : fS = pdfConfig.pop('f_S')
        elif isinstance( ASMag2, RooObject ) : fS = ASMag2.GetVal() / ( 1. + ASMag2.GetVal() )
        else                                 : fS = ASMag2 / ( 1. + ASMag2 )

        # CP violation parameters
        carthLambdaCP = pdfConfig.pop('carthLambdaCP')
        phiCP         = pdfConfig.pop('phiCP')
        lambdaCPSq    = pdfConfig.pop('lambdaCPSq')

        # B lifetime parameters
        Gamma  = pdfConfig.pop('Gamma')
        dGamma = pdfConfig.pop('deltaGamma')
        dM     = pdfConfig.pop('deltaM')

        # asymmetries
        AProd = pdfConfig.pop('AProd')

        # plots
        angleNames = pdfConfig.pop('angleNames')
        numBkgAngleBins = pdfConfig['numBkgAngleBins']


        ###################################################################################################################################
        ## create variables (except for tagging category) and read real data ##
        #######################################################################

        # RooObject wrappers
        from RooFitWrappers import RooObject, RealVar, Category

        # angular functions
        if nominalPdf :
            from P2VVParameterizations.AngularFunctions import JpsiphiTransversityAngles as AngleFuncs
            self._angleFuncs = AngleFuncs( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
        else :
            from P2VVParameterizations.AngularFunctions import JpsiphiHelicityAngles as AngleFuncs
            self._angleFuncs = AngleFuncs( cpsi = 'helcosthetaK', ctheta = 'helcosthetaL', phi = 'helphi' )

        # variables in PDF (except for tagging category)
        time = RealVar(  'time', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. )
                             , Ranges = dict( Bulk = ( None, 5. ) )
                            )
        cpsi   = self._angleFuncs.angles['cpsi']
        ctheta = self._angleFuncs.angles['ctheta']
        phi    = self._angleFuncs.angles['phi']
        iTag = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'B' : +1, 'Bbar' : -1 } )

        BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True
                        , Value = 5368., MinMax = ( 5200., 5550. ), nBins = 48
                        ,  Ranges = dict(  LeftSideBand  = ( None,  5330. )
                                         , Signal        = ( 5330., 5410. )
                                         , RightSideBand = ( 5410., None  )
                                        )
                       )

        angles = [ cpsi, ctheta, phi ]
        obsSetP2VV = [ time ] + angles + [ iTag ]
        if not SFit : obsSetP2VV += [ BMass ]

        # ntuple variables
        mpsi = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. ), nBins =  32 )
        mphi = RealVar( 'mdau2', Title = 'M(KK)',     Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. ), nBins =  16 )

        timeRes = RealVar( 'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = ( 0., 0.12 ), nBins = 50 )
        timeRes.setBins( 20, 'cache' )

        tagDecision = Category( 'tagdecision_os', Title = 'Tag decision', Observable = True
                               , States = { 'B' : +1, 'Bbar' : -1 , 'Untagged' : 0 }
                              )
        estWTag = RealVar( 'tagomega_os', Title = 'Estimated wrong tag probability', Observable = True
                           , Value = 0.25, MinMax = ( 0., 0.50001 ) )
        estWTag.setBins( 20, 'cache' )
        tagCat = Category( 'tagcat_os',   Title = 'Tagging Category', Observable = True
                          , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                         )

        sel   = Category( 'sel',             Title = 'Selection',        Observable = True, States = { 'selected' : +1 } )
        trig  = Category( 'triggerDecision', Title = 'Trigger Decision', Observable = True, States = { 'selected' : +1 } )

        observables = {  time.GetName()        : time
                       , angles[0].GetName()   : angles[0]
                       , angles[1].GetName()   : angles[1]
                       , angles[2].GetName()   : angles[2]
                       , iTag.GetName()        : iTag
                       , tagDecision.GetName() : tagDecision
                       , estWTag.GetName()     : estWTag
                       , tagCat.GetName()      : tagCat
                       , BMass.GetName()       : BMass
                       , mpsi.GetName()        : mpsi
                       , mphi.GetName()        : mphi
                       , timeRes.GetName()     : timeRes
                       , sel.GetName()         : sel
                       , trig.GetName()        : trig
                      }

        obsSetNTuple = [ time ] + angles +  [ BMass, mpsi, mphi, timeRes ] + [ tagDecision, estWTag, tagCat ] + [ sel, trig ]

        # read ntuple
        from P2VVGeneralUtils import readData
        self._data = readData(  filePath = nTupleFile, dataSetName = nTupleName, NTuple = True, observables = obsSetNTuple
                              , Rename = 'JpsiphiData' )

        # get data in signal and side band ranges
        self._sigRangeData = self._data.reduce( CutRange = 'Signal'       )
        self._bkgRangeData = self._data.reduce( CutRange = 'LeftSideBand' )
        self._bkgRangeData.append( self._data.reduce( CutRange = 'RightSideBand' ) )


        ###################################################################################################################################
        ## build tagging categories ##
        ##############################

        # build tagging categories
        from P2VVParameterizations.FlavourTagging import Linear_TaggingCategories as TaggingCategories
        if tagConds == 'estWTag' :
            self._tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = self._data, estWTag = estWTag
                                              , wTagP0Constraint = True, wTagP1Constraint = True )
        else :
            self._tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = self._data, estWTagName = estWTag.GetName()
                                              , TagCats = [ ], NumSigmaTagBins = 1., wTagP0Constraint = True, wTagP1Constraint = True )

        tagCatP2VV = self._tagCats['tagCat']
        observables[tagCatP2VV.GetName()] = tagCatP2VV
        obsSetP2VV.append( tagCatP2VV )

        # add tagging category to data set
        from P2VVGeneralUtils import addTaggingObservables
        addTaggingObservables(  self._data, iTag.GetName(), tagCatP2VV.GetName(), tagDecision.GetName(), estWTag.GetName()
                              , self._tagCats['tagCats'] )


        ###################################################################################################################################
        ## initialize PDF component objects ##
        ######################################

        nEvents     = self._data.numEntries()
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
            massPdf = buildPdf( [ sigMassComps, bkgMassComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )

        else :
            signalComps     += signalBMass.pdf()
            backgroundComps += backgroundBMass.pdf()
            massPdf = buildPdf( [ signalComps, backgroundComps ], Observables = [ BMass ], Name = 'JpsiphiMass' )


        ###################################################################################################################################
        ## compute S-weights ##
        #######################

        print 120 * '='
        print 'Bs2Jpsiphi_PdfBuilder: computing S-weights'

        from P2VVGeneralUtils import SData, splot
        massPdf.fitTo( self._data, **fitOpts )
        #for par in massPdf.Parameters() : par.setConstant( not par.getAttribute('Yield') )
        self._SWeightData = SData( Pdf = massPdf, Data = self._data, Name = 'massSData' )
        self._sigSWeightData = self._SWeightData.data( 'sigMass' if SFit else 'signal' )
        self._bkgSWeightData = self._SWeightData.data( 'bkgMass' if SFit else 'bkg'    )

        print 120 * '=' + '\n'


        ###################################################################################################################################
        ## time acceptance function ##
        ##############################

        if multiplyByTimeEff in [ 'all', 'signal', 'background' ] :
            from P2VVParameterizations.TimeAcceptance import Moriond2012_TimeAcceptance as TimeAcceptance
            timeAcceptance = TimeAcceptance(  time = time, Input = timeEffHistFile, Histogram = timeEffHistName )


        ###################################################################################################################################
        ## build the B_s -> J/psi phi signal time, angular and tagging PDF ##
        #####################################################################

        # transversity amplitudes
        if nominalPdf :
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(  A0Mag2    = A0Mag2
                                    , A0Phase   = A0Ph
                                    , AperpMag2 = AperpMag2
                                    , AparPhase = AparPh
                                    , f_S       = fS
                                    , ASPhase   = ASPh
                                   )

        elif amplitudeParam == 'phasesSWaveFrac' :
            from math import cos, sin, sqrt
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(  A0Mag2    = A0Mag2
                                    , A0Phase   = A0Ph
                                    , AperpMag2 = AperpMag2
                                    , AparPhase = AparPh
                                    , sqrtfS_Re = sqrt(fS)* cos(ASPh)
                                    , sqrtfS_Im = sqrt(fS)* sin(ASPh)
                                   )

        else :
            from math import sqrt, cos, sin
            from P2VVParameterizations.DecayAmplitudes import JpsiVCarthesian_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(  ReApar  = sqrt(AparMag2  / A0Mag2) * cos(AparPh)
                                    , ImApar  = sqrt(AparMag2  / A0Mag2) * sin(AparPh)
                                    , ReAperp = sqrt(AperpMag2 / A0Mag2) * cos(AperpPh)
                                    , ImAperp = sqrt(AperpMag2 / A0Mag2) * sin(AperpPh)
                                    , ReAS    = sqrt(ASMag2    / A0Mag2) * cos(ASPh)
                                    , ImAS    = sqrt(ASMag2    / A0Mag2) * sin(ASPh)
                                   )

        # B lifetime
        from P2VVParameterizations.LifetimeParams import Gamma_LifetimeParams as LifetimeParams
        dGammaVar = dict( Name = 'dGamma', Value = dGamma )
        if blind : dGammaVar['Blind'] = ( 'UnblindUniform', 'BsRooBarbMoriond2012', 0.02 )
        lifetimeParams = LifetimeParams( Gamma = dict(Value = Gamma), deltaGamma = dGammaVar, deltaMConstraint = True )

        from P2VVParameterizations.TimeResolution import LP2011_TimeResolution as TimeResolution
        timeResModel = TimeResolution( time = time, timeResSFConstraint = True )

        # CP violation parameters
        from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
        phiCPVar = dict( Value = phiCP )
        if blind: phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
        self._lambdaCP = CPParam( lambdaCPSq = lambdaCPSq, phiCP = phiCPVar )

        # tagging parameters
        from P2VVParameterizations.FlavourTagging import CatDilutionsCoefAsyms_TaggingParams as TaggingParams
        self._taggingParams = TaggingParams( AProd = AProd, ANorm = -self._lambdaCP['C'].getVal(), **self._tagCats.tagCatsDict() )

        # coefficients for time functions
        from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

        # build signal PDF
        args = dict(  time                   = time
                    , iTag                   = iTag
                    , tagCat                 = tagCatP2VV
                    , tau                    = lifetimeParams['MeanLifetime']
                    , dGamma                 = lifetimeParams['deltaGamma']
                    , dm                     = lifetimeParams['deltaM']
                    , dilutions              = self._taggingParams['dilutions']
                    , ADilWTags              = self._taggingParams['ADilWTags']
                    , avgCEvens              = self._taggingParams['avgCEvens']
                    , avgCOdds               = self._taggingParams['avgCOdds']
                    , tagCatCoefs            = self._taggingParams['tagCatCoefs']
                    , coshCoef               = timeBasisCoefs['cosh']
                    , sinhCoef               = timeBasisCoefs['sinh']
                    , cosCoef                = timeBasisCoefs['cos']
                    , sinCoef                = timeBasisCoefs['sin']
                    , resolutionModel        = timeResModel['model']
                    , ConditionalObservables = timeResModel.conditionalObservables() + self._tagCats.conditionalObservables()
                    , ExternalConstraints    = lifetimeParams.externalConstraints()\
                                               + timeResModel.externalConstraints()\
                                               + self._tagCats.externalConstraints()
                   )

        from RooFitWrappers import BTagDecay
        sig_t_angles_tagCat_iTag = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

        if angEffMomentsFile :
            # multiply signal PDF with angular efficiency
            print 120 * '='
            print 'Bs2Jpsiphi_PdfBuilder: multiplying signal PDF with angular efficiency moments from file "%s"' % angEffMomentsFile

            from P2VVGeneralUtils import RealMomentsBuilder
            moments = RealMomentsBuilder()
            moments.appendPYList( self._angleFuncs.angles, [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ) ] if not nominalPdf \
                                                else [ ( 0, 0, 0 ), ( 2, 0, 0 ), ( 0, 2, 0 ), ( 0, 2, 2 ) ]
                                )
            moments.read(angEffMomentsFile)
            moments.Print()
            sig_t_angles_tagCat_iTag = moments * sig_t_angles_tagCat_iTag

            print 120 * '=' + '\n'

        if multiplyByTimeEff in [ 'all', 'signal' ] :
            # multiply signal PDF with time acceptance
            sig_t_angles_tagCat_iTag = timeAcceptance * sig_t_angles_tagCat_iTag

        signalComps += sig_t_angles_tagCat_iTag

        if tagConds == 'estWTag' :
            # tagging category is used as conditional observable: register dummy PDF for tagging category
            from RooFitWrappers import UniformPdf
            signalComps += UniformPdf( 'tagCatDummyPdf', Arguments = [ tagCatP2VV ] )


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        if not SFit :
            from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
            backgroundTime = BackgroundTime( Name = 'bkg_t', time = time, resolutionModel = timeResModel['model'] )
            bkg_t = backgroundTime.pdf()

            if multiplyByTimeEff in [ 'all', 'background' ] :
                # multiply background time PDF with time acceptance
                bkg_t = timeAcceptance * bkg_t

            backgroundComps += bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit :
            print 120 * '='
            print 'Bs2Jpsiphi_PdfBuilder: determining angular shape of background from %s data'\
                  % ( 'B mass side band' if massRangeBackground else 'B mass S-weight' )

            bkgAngleData = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
            if bkgAnglePdf == 'histPdf' :
                # create a HistPdf from background data
                from RooFitWrappers import HistPdf
                self._bkg_angles = HistPdf(  Name = 'bkg_angles'
                                           , Observables = angles
                                           , Binning = { cpsi : numBkgAngleBins[0], ctheta : numBkgAngleBins[1], phi : numBkgAngleBins[2] }
                                           , Data = bkgAngleData
                                          )

            else :
                # define angular bins
                from math import pi
                from array import array
                if nominalPdf :
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
                                              , Title = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                                              , Value = 1. / cpsNumBins / cthNumBins / phiNumBins
                                              , MinMax = ( 0., 1. )
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

                # determine bin coefficient values
                sumWeights = 0.
                sumBinWeights = ( cpsNumBins * cthNumBins * phiNumBins - 1 ) * [ 0. ]
                for obsSet in bkgAngleData :
                    for angle in angles : angle.setVal( obsSet.getRealValue( angle.GetName() ) )
                    sumWeights += bkgAngleData.weight()
                    cpsBin = cpsi.getBin('bkg_cpsBins')
                    cthBin = ctheta.getBin('bkg_cthBins')
                    phiBin = phi.getBin('bkg_phiBins')
                    bin = cpsBin + cthBin * cpsNumBins + phiBin * cpsNumBins * cthNumBins - 1
                    if bin >= 0 : sumBinWeights[ bin ] += bkgAngleData.weight()

                # set bin coefficient values
                for coef, weight in zip( self._bkgAngCoefs, sumBinWeights ) :
                    value = weight / sumWeights
                    assert value >= 0.,\
                        'Bs2Jpsiphi_PdfBuilder: background angular PDF coefficient \"%s\" has negative value: %f' % (coef.GetName(), value)
                    coef.setVal(value)
                    coef.setConstant()

            backgroundComps += self._bkg_angles

            if makePlots :
                # import plotting tools
                from P2VVLoad import ROOTStyle
                from P2VVGeneralUtils import plot
                from ROOT import TCanvas, kBlue, kRed, kGreen, kDashed

                # plot background angles
                self._bkgAnglesCanv = TCanvas( 'bkgAnglesCanv', 'Background Decay Angles' )
                for ( pad, obs, data, nBins, plotTitle, xTitle )\
                      in zip(  self._bkgAnglesCanv.pads( 3, 2 )
                             , 2 * angles
                             , 3 * ( self._bkgRangeData, ) + 3 * ( self._bkgSWeightData, )
                             , 2 * numBkgAngleBins
                             , [ angle.GetTitle() + ' - mass side bands' for angle in angles ]
                             + [ angle.GetTitle() + ' - mass S-weights'  for angle in angles ]
                             , 2 * ( angleNames[0][1], angleNames[1][1], angleNames[2][1] )
                            ) :
                    plot(  pad, obs, data, self._bkg_angles, xTitle = xTitle
                         , frameOpts  = dict( Bins = nBins, Title = plotTitle   )
                         , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4 )
                         , pdfOpts    = dict( LineColor = kBlue, LineWidth = 2  )
                        )

            print 120 * '=' + '\n'


        ###################################################################################################################################
        ## build background tagging PDF ##
        ##################################

        if not SFit :
            from RooFitWrappers import BinnedPdf
            backgroundComps += BinnedPdf(  'bkg_tagCat_iTag'
                                         , Categories   = ( tagCatP2VV, iTag )
                                         , Coefficients = [  self._taggingParams['tagCatCoefs']
                                                           , [  RealVar( 'bkg_BbarFrac'
                                                                        , Title    = 'Anti-B fraction in background'
                                                                        , Value    = 0.5
                                                                        , MinMax = ( 0., 1. )
                                                                        , Constant = True
                                                                       )
                                                             ]
                                                          ]
                                        )


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

        PdfBuilder.__init__( self, pdf, observables, { } )

