###########################################################################################################################################
## P2VVParameterizations.FullPDFs: Parameterizations of complete PDFs that are used in an analysis                                       ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################

Bs2Jpsiphi_Winter2012 = {}

class PdfBuilder ( object ) :
    def __init__( self, pdf ) : self._pdf = pdf
    def __getitem__( self, kw ) : return getattr( self, '_' + kw )
    def pdf(self) : return self._pdf


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
        if 'f_S' in pdfConfig                   : fS = pdfConfig.pop('f_S')
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
        self._time = RealVar(  'time', Title = 'Decay time', Unit = 'ps', Observable = True, Value = 0.5, MinMax = ( 0.3, 14. )
                             , Ranges = dict( Bulk = ( None, 5. ) )
                            )
        self._iTag = Category( 'iTag', Title = 'Initial state flavour tag', Observable = True, States = { 'B' : +1, 'Bbar' : -1 } )
        angles = [ self._angleFuncs.angles['cpsi'], self._angleFuncs.angles['ctheta'], self._angleFuncs.angles['phi'] ]

        BMass = RealVar( 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV', Observable = True
                        , Value = 5368., MinMax = ( 5200., 5550. ), nBins = 48
                        ,  Ranges = dict(  LeftSideBand  = ( None,  5330. )
                                         , Signal        = ( 5330., 5410. )
                                         , RightSideBand = ( 5410., None  )
                                        )
                       )

        self._obsSetP2VV = [ self._time ] + angles + [ self._iTag ]
        if not SFit : self._obsSetP2VV += [ BMass ]

        # ntuple variables
        mpsi = RealVar( 'mdau1', Title = 'M(#mu#mu)', Unit = 'MeV', Observable = True, MinMax = ( 3090. - 60., 3090. + 60. ), nBins =  32 )
        mphi = RealVar( 'mdau2', Title = 'M(KK)',     Unit = 'MeV', Observable = True, MinMax = ( 1020. - 12., 1020. + 12. ), nBins =  16 )

        timeRes = RealVar( 'sigmat', Title = '#sigma(t)', Unit = 'ps', Observable = True, MinMax = ( 0., 0.12 ), nBins =  50 )

        tagDecision = Category( 'tagdecision_os', Title = 'Tag decision', Observable = True
                               , States = { 'B' : +1, 'Bbar' : -1 , 'Untagged' : 0 }
                              )
        tagOmega = RealVar( 'tagomega_os', Title = 'Estimated wrong tag', Observable = True, Value = 0.25, MinMax = ( 0., 0.50001 ) )
        tagCat = Category( 'tagcat_os',   Title = 'Tagging Category', Observable = True
                          , States = [ 'Untagged' ] + [ 'TagCat%d' % cat for cat in range( 1, 6 ) ]
                         )

        sel   = Category( 'sel',             Title = 'Selection',        Observable = True, States = { 'selected' : +1 } )
        trig  = Category( 'triggerDecision', Title = 'Trigger Decision', Observable = True, States = { 'selected' : +1 } )

        obsSetNTuple = [ self._time ] + angles +  [ BMass, mpsi, mphi, timeRes ] + [ tagDecision, tagOmega, tagCat ] + [ sel, trig ]

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
        tagCats = TaggingCategories(  tagCat = 'tagCatP2VV', DataSet = self._data, estWTagName = tagOmega.GetName(), TagCats = [ ]
                                    , NumSigmaTagBins = 1., wTagP0Constraint = True, wTagP1Constraint = True )
        tagCatP2VV = tagCats['tagCat']
        self._obsSetP2VV.append( tagCatP2VV )

        # add tagging category to data set
        from P2VVGeneralUtils import addTaggingObservables
        addTaggingObservables(  self._data, self._iTag.GetName(), tagCatP2VV.GetName(), tagDecision.GetName(), tagOmega.GetName()
                              , tagCats['tagCats'] )

        # tagging parameters
        numTagCats    = tagCats['numTagCats']
        tagCat5Min    = tagCats.traditionalCatRange(5)[0]
        taggedCatsStr = ','.join( [ 'TagCat%d' % cat for cat in range( 1,          numTagCats ) ] )
        tagCat5Str    = ','.join( [ 'TagCat%d' % cat for cat in range( tagCat5Min, numTagCats ) ] )

        # tagging category ranges
        tagCatP2VV.setRange( 'UntaggedRange', 'Untagged'    )
        tagCatP2VV.setRange( 'TaggedRange',   taggedCatsStr )
        tagCatP2VV.setRange( 'TagCat5Range',  tagCat5Str    )


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
            timeAcceptance = TimeAcceptance(  time = self._time, Input = timeEffHistFile, Histogram = timeEffHistName )


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
            from math import cos, sin
            from P2VVParameterizations.DecayAmplitudes import JpsiVPolarSWaveFrac_AmplitudeSet as Amplitudes
            amplitudes = Amplitudes(  A0Mag2    = A0Mag2
                                    , A0Phase   = A0Ph
                                    , AperpMag2 = AperpMag2
                                    , AparPhase = AparPh
                                    , f_S_Re    = fS * cos(ASPh)
                                    , f_S_Im    = fS * sin(ASPh)
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
        timeResModel = TimeResolution( time = self._time, timeResSFConstraint = True )

        # CP violation parameters
        from P2VVParameterizations.CPVParams import LambdaSqArg_CPParam as CPParam
        phiCPVar = dict( Value = phiCP )
        if blind: phiCPVar['Blind'] = ( 'UnblindUniform', 'BsCustardMoriond2012', 0.3 )
        self._lambdaCP = CPParam( lambdaCPSq = lambdaCPSq, phiCP = phiCPVar )

        # tagging parameters
        from P2VVParameterizations.FlavourTagging import WTagCatsCoefAsyms_TaggingParams as TaggingParams
        self._taggingParams = TaggingParams( AProd = AProd, ANorm = -self._lambdaCP['C'].getVal(), **tagCats.tagCatsDict() )

        # coefficients for time functions
        from P2VVParameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients as TimeBasisCoefs
        timeBasisCoefs = TimeBasisCoefs( self._angleFuncs.functions, amplitudes, self._lambdaCP, [ 'A0', 'Apar', 'Aperp', 'AS' ] ) 

        # build signal PDF
        args = dict(  time                = self._time
                    , iTag                = self._iTag
                    , tagCat              = tagCatP2VV
                    , tau                 = lifetimeParams['MeanLifetime']
                    , dGamma              = lifetimeParams['deltaGamma']
                    , dm                  = lifetimeParams['deltaM']
                    , dilutions           = self._taggingParams['dilutions']
                    , ADilWTags           = self._taggingParams['ADilWTags']
                    , avgCEvens           = self._taggingParams['avgCEvens']
                    , avgCOdds            = self._taggingParams['avgCOdds']
                    , tagCatCoefs         = self._taggingParams['tagCatCoefs']
                    , coshCoef            = timeBasisCoefs['cosh']
                    , sinhCoef            = timeBasisCoefs['sinh']
                    , cosCoef             = timeBasisCoefs['cos']
                    , sinCoef             = timeBasisCoefs['sin']
                    , resolutionModel     = timeResModel['model']
                    , ExternalConstraints = lifetimeParams.externalConstraints()\
                                            + timeResModel.externalConstraints()\
                                            + tagCats.externalConstraints()
                   )

        from RooFitWrappers import BTagDecay
        sig_t_angles_tagCat_iTag = BTagDecay( 'sig_t_angles_tagCat_iTag', **args )

        if angEffMomentsFile :
            # multiply signal PDF with angular efficiency
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


        ###################################################################################################################################
        ## build background time PDF ##
        ###############################

        if not SFit :
            from P2VVParameterizations.TimePDFs import LP2011_Background_Time as BackgroundTime
            backgroundTime = BackgroundTime( Name = 'bkg_t', time = self._time, resolutionModel = timeResModel['model'] )
            bkg_t = backgroundTime.pdf()

            if multiplyByTimeEff in [ 'all', 'background' ] :
                # multiply background time PDF with time acceptance
                bkg_t = timeAcceptance * bkg_t

            backgroundComps += bkg_t


        ###################################################################################################################################
        ## build background angular PDF ##
        ##################################

        if not SFit :
            if bkgAnglePdf == 'histPdf' :
                # create a HistPdf from background data
                from RooFitWrappers import HistPdf
                self._bkg_angles = HistPdf(  Name = 'bkg_angles'
                                           , Observables = angles
                                           , Binning =  {  self._angleFuncs.angles['cpsi']   : numBkgAngleBins[0]
                                                         , self._angleFuncs.angles['ctheta'] : numBkgAngleBins[1]
                                                         , self._angleFuncs.angles['phi' ]   : numBkgAngleBins[2]
                                                        }
                                           , Data = self._bkgRangeData if massRangeBackground else self._bkgSWeightData
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

                from ROOT import RooBinning
                cpsBins = RooBinning( len(cpsBinBounds) - 1, cpsBinBounds, 'bkg_cpsBins' )
                cthBins = RooBinning( len(cthBinBounds) - 1, cthBinBounds, 'bkg_cthBins' )
                phiBins = RooBinning( len(phiBinBounds) - 1, phiBinBounds, 'bkg_phiBins' )
                self._angleFuncs.angles['cpsi'].setBinning(   cpsBins, 'bkg_cpsBins' )
                self._angleFuncs.angles['ctheta'].setBinning( cthBins, 'bkg_cthBins' )
                self._angleFuncs.angles['phi'].setBinning(    phiBins, 'bkg_phiBins' )

                # create bin coefficients
                self._bkgAngCoefs = [ RealVar(  'bkg_angBin_%d_%d_%d' % ( bin0, bin1, bin2 )
                                              , Title = 'Background angles bin %d-%d-%d' % ( bin0, bin1, bin2 )
                                              , Value = 1. / ( len(cpsBinBounds) - 1 ) / ( len(cthBinBounds) - 1 ) / ( len(phiBinBounds) - 1 )
                                              , MinMax = ( 0., 1. )
                                             )\
                                      if bin0 != 0 or bin1 != 0 or bin2 != 0 else None\
                                      for bin2 in range( len(phiBinBounds) - 1 )\
                                      for bin1 in range( len(cthBinBounds) - 1 )\
                                      for bin0 in range( len(cpsBinBounds) - 1 )
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
                self._bkg_angles.fitTo( self._bkgRangeData if massRangeBackground else self._bkgSWeightData, SumW2Error = False )

                ## set bin coefficient values
                #for bin, coef in enumerate(self._bkgAngCoefs) :
                #    cpsBin, cthBin, phiBin = coef.GetName()[ 11 : ].split('_')
                #    selection = dict(  cps = angles[0].GetName(), cth = angles[1].GetName(), phi = angles[2].GetName()
                #                     , cpsMin = cpsBinBounds[ int(cpsBin) ], cpsMax = cpsBinBounds[ int(cpsBin) + 1 ]
                #                     , cthMin = cthBinBounds[ int(cthBin) ], cthMax = cthBinBounds[ int(cthBin) + 1 ]
                #                     , phiMin = phiBinBounds[ int(phiBin) ], phiMax = phiBinBounds[ int(phiBin) + 1 ]
                #                    )

                #    value = float( sideBandDataTree.GetEntries( (     '%(cps)s >= %(cpsMin)f && %(cps)s < %(cpsMax)f'\
                #                                                 + '&& %(cth)s >= %(cthMin)f && %(cth)s < %(cthMax)f'\
                #                                                 + '&& %(phi)s >= %(phiMin)f && %(phi)s < %(phiMax)f'
                #                                                ) % selection
                #                                              )
                #                 ) / sideBandDataTree.GetEntries()
                #    coef.setVal(value)
                #    coef.setError( 0.1 * value )
                #    coef.setConstant()

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
                             , angles
                             , 3 * ( self._bkgRangeData, ) + 3 * ( self._bkgSWeightData, )
                             , 2 * numBkgAngleBins
                             , [ angle.GetTitle() + ' - mass side bands' for angle in angles ]
                             + [ angle.GetTitle() + ' - mass S-weights'  for angle in angles ]
                             , 2 * angleNames
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
            from RooFitWrappers import BinnedPdf
            backgroundComps += BinnedPdf(  'bkg_tagCat_iTag'
                                         , Categories   = ( tagCatP2VV, self._iTag )
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
            pdf = buildPdf( [ signalComps ], Observables = self._obsSetP2VV, Name = 'Jpsiphi' )

        else :
            if components == 'signal' :
                pdf = buildPdf( [ signalComps                  ], Observables = self._obsSetP2VV, Name = 'JpsiphiSig' )
            elif components == 'background' :
                pdf = buildPdf( [ backgroundComps              ], Observables = self._obsSetP2VV, Name = 'JpsiphiBkg' )
            else :
                pdf = buildPdf( [ signalComps, backgroundComps ], Observables = self._obsSetP2VV, Name = 'Jpsiphi'    )

        PdfBuilder.__init__( self, pdf )

