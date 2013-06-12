###########################################################################################################################################
## script settings ##
#####################

nTupleFilePath   = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
nTupleName       = 'DecayTree'
dataSetsFilePath = 'P2VVDataSets.root'
plotsFilePath    = 'P2VVMassPlots.ps'

selection        = 'paper2012'
dataSample       = ''
addTaggingObs    = True
KKMassBinBounds  = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ]

ntupleCuts = 'sel == 1'\
             ' && muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
selections = dict(  HLT1Unbiased   = 'hlt1_unbiased_dec==1 && hlt2_biased==1'
                  , HLT1ExclBiased = 'hlt1_excl_biased_dec==1 && hlt2_biased==1'
                  , paper2012      = '(hlt1_biased==1 || hlt1_unbiased_dec==1) && hlt2_biased==1'
                  , timeEffFit     = '(hlt1_biased==1 || hlt1_unbiased_dec==1) && (hlt2_biased==1 || hlt2_unbiased==1)'
                 )

sigFrac      = 0.504
sigMassModel = ''
bkgMassModel = ''
SWeightsType = 'simultaneousFreeBkg'
numMassBins  = [ 70, 40, 20, 20, 20 ]

from math import pi
from ROOT import RooNumber
RooInf  = RooNumber.infinity()
KKMMin  = KKMassBinBounds[0]
KKMMax  = KKMassBinBounds[-1]
obsKeys = [  'mass', 'mumuMass', 'KKMass', 'time', 'timeRes', 'ctk', 'ctl', 'phih', 'cpsi', 'cttr', 'phitr'
           , 'wTag', 'wTagOS', 'wTagSS', 'tagDec', 'tagDecOS', 'tagDecSS', 'tagCatOS'
           #, 'sel', 'selA', 'selB'
           , 'hlt1ExclB', 'hlt1B', 'hlt1UB', 'hlt2B', 'hlt2UB'
          ]

obsDict = dict(  mass      = ( True,  'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.       )
               , mumuMass  = ( True,  'mdau1',                'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.       )
               , KKMass    = ( True,  'mdau2',                'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,  KKMMin,  KKMMax      )
               , time      = ( True,  'time',                 'Decay time',              'ps',        1.5,    0.3,     14.         )
               , timeRes   = ( True,  'sigmat',               '#sigma(t)',               'ps',        0.01,   0.0001,  0.12        )
               , ctk       = ( True,  'helcthetaK',           'cos(#theta_{K})',         '',          0.,    -1.,     +1.          )
               , ctl       = ( True,  'helcthetaL',           'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.          )
               , phih      = ( True,  'helphi',               '#phi_{h}',                'rad',       0.,    -pi,     +pi          )
               , cpsi      = ( True,  'trcpsi',               'cos(#psi_{tr})',          '',          0.,    -1.,     +1.          )
               , cttr      = ( True,  'trctheta',             'cos(#theta_{tr})',        '',          0.,    -1.,     +1.          )
               , phitr     = ( True,  'trphi',                '#phi_{tr}',               'rad',       0.,    -pi,     +pi          )
               , wTag      = ( True,  'tagomega',             'est. wrong-tag prob.',    '',          0.25,   0.,      0.50001     )
               , wTagOS    = ( True,  'tagomega_os',          'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , wTagSS    = ( True,  'tagomega_ss',          'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , tagDec    = ( True,  'tagdecision',          'tag decision',    { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecOS  = ( True,  'tagdecision_os',       'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecSS  = ( True,  'tagdecision_ss',       'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagCatOS  = ( True,  'tagcat_os',            'OS tag category', [ 'unt' ] + [ 'cat%d' % c for c in range(1, 6) ]  )
               , sel       = ( True,  'sel',                  'selection',       { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selA      = ( True,  'selA',                 'selection A',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selB      = ( True,  'selB',                 'selection B',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , hlt1ExclB = ( True,  'hlt1_excl_biased_dec', 'HLT1 excl. B.',    { 'exclB' : 1, 'notExclB' : 0 }                  )
               , hlt1B     = ( True,  'hlt1_biased',          'HLT1 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt1UB    = ( True,  'hlt1_unbiased_dec',    'HLT1 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , hlt2B     = ( True,  'hlt2_biased',          'HLT2 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt2UB    = ( True,  'hlt2_unbiased',        'HLT2 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , sel_onecand                = ( True, 'sel_onecand',                { 'sel' : 1, 'notSel' : 0 }              )
               , sel_one_gl                 = ( True, 'sel_one_gl',                 { 'sel' : 1, 'notSel' : 0 }              )
               , muPTrChi2                  = ( True, 'muplus_track_chi2ndof',      'mu+ chi^2/#dof', 1.,    0.,      4.     )
               , muMTrChi2                  = ( True, 'muminus_track_chi2ndof',     'mu- chi^2/#dof', 1.,    0.,      4.     )
               , KPTrChi2                   = ( True, 'Kplus_track_chi2ndof',       'K+ chi^2/#dof',  1.,    0.,      4.     )
               , KMTrChi2                   = ( True, 'Kminus_track_chi2ndof',      'K- chi^2/#dof',  1.,    0.,      4.     )
               , GLsb                       = ( True, 'GLsb',                       'GLsb',           0.,    0.,      1.     )
               , muplus_PIDmu               = ( True, 'muplus_PIDmu',               'muplus_PIDmu',   0.,   -RooInf, +RooInf )
               , muminus_PIDmu              = ( True, 'muminus_PIDmu',              'muminus_PIDmu',  0.,   -RooInf, +RooInf )
               , Kplus_pidK                 = ( True, 'Kplus_pidK',                 'Kplus_pidK',     0.,    0.,     +RooInf )
               , Kminus_pidK                = ( True, 'Kminus_pidK',                'Kminus_pidK',    0.,    0.,     +RooInf )
               , muplus_PX                  = ( True, 'muplus_PX',                  'muplus_PX',      0.,   -RooInf, +RooInf )
               , muplus_PY                  = ( True, 'muplus_PY',                  'muplus_PY',      0.,   -RooInf, +RooInf )
               , muplus_PZ                  = ( True, 'muplus_PZ',                  'muplus_PZ',      0.,   -RooInf, +RooInf )
               , muminus_PX                 = ( True, 'muminus_PX',                 'muminus_PX',     0.,   -RooInf, +RooInf )
               , muminus_PY                 = ( True, 'muminus_PY',                 'muminus_PY',     0.,   -RooInf, +RooInf )
               , muminus_PZ                 = ( True, 'muminus_PZ',                 'muminus_PZ',     0.,   -RooInf, +RooInf )
               , Kplus_PX                   = ( True, 'Kplus_PX',                   'Kplus_PX',       0.,   -RooInf, +RooInf )
               , Kplus_PY                   = ( True, 'Kplus_PY',                   'Kplus_PY',       0.,   -RooInf, +RooInf )
               , Kplus_PZ                   = ( True, 'Kplus_PZ',                   'Kplus_PZ',       0.,   -RooInf, +RooInf )
               , Kminus_PX                  = ( True, 'Kminus_PX',                  'Kminus_PX',      0.,   -RooInf, +RooInf )
               , Kminus_PY                  = ( True, 'Kminus_PY',                  'Kminus_PY',      0.,   -RooInf, +RooInf )
               , Kminus_PZ                  = ( True, 'Kminus_PZ',                  'Kminus_PZ',      0.,   -RooInf, +RooInf )
               , B_Pt                       = ( True, 'B_Pt',                       'B_Pt',           0.,    0.,      RooInf )
               , phi_1020_pt                = ( True, 'phi_1020_pt',                'phi_1020_pt',    500.,  500.,    RooInf )
               , B_s0_LOKI_CosPolAngle_Dau1 = ( True, 'B_s0_LOKI_CosPolAngle_Dau1', 'mumu cos(th)',   0.,   -1.,     +1.     )
               , B_s0_IP_OWNPV              = ( True, 'B_s0_IP_OWNPV',              'B_s0_IP_OWNPV',  0.,   -RooInf, +RooInf )
               , B_s0_IPCHI2_OWNPV          = ( True, 'B_s0_IPCHI2_OWNPV',          'IP chi2 PV',     0.,   -RooInf, +RooInf )
               , B_s0_MINIPCHI2NEXTBEST     = ( True, 'B_s0_MINIPCHI2NEXTBEST',     'IP chi2 next',   0.,   -RooInf, +RooInf )
               , B_s0_LOKI_DTF_VCHI2NDOF    = ( True, 'B_s0_LOKI_DTF_VCHI2NDOF',    'DTF chi2/#dof',  0.,   -RooInf, +RooInf )
               , B_s0_ENDVERTEX_CHI2        = ( True, 'B_s0_ENDVERTEX_CHI2',        'Bs0 vert chi2',  0.,   0.,       50.    )
               , phi_1020_ENDVERTEX_CHI2    = ( True, 'phi_1020_ENDVERTEX_CHI2',    'mumu vert chi2', 0.,   0.,       16.    )
               , J_psi_1S_ENDVERTEX_CHI2    = ( True, 'J_psi_1S_ENDVERTEX_CHI2',    'KK vert chi2',   0.,   0.,       16.    )
              )

massRanges = dict(  LeftSideBand  = ( 5200., 5320. )
                  , Signal        = ( 5320., 5420. )
                  , RightSideBand = ( 5420., 5550. )
                  , PeakBkg       = ( 5390., 5440. )
                 )

fitOpts = dict(  NumCPU    = 6
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Offset    = True
              )


###########################################################################################################################################
## read data ##
###############

from P2VV.Load import RooFitOutput, LHCbStyle

# create workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# create observables
observables  = { }
obsSetNTuple = [ ]
from P2VV.RooFitWrappers import RealVar, Category
for obs in obsKeys :
    if type( obsDict[obs][3] ) == dict or type( obsDict[obs][3] ) == list :
        observables[obs] = Category( obsDict[obs][1], Title = obsDict[obs][2], Observable = True, States = obsDict[obs][3] )
    else :
        observables[obs] = RealVar( obsDict[obs][1], Title = obsDict[obs][2], Unit = obsDict[obs][3], Observable = True
                                   , Value = obsDict[obs][4], MinMax = ( obsDict[obs][5], obsDict[obs][6] ) )

    if obsDict[obs][0] : obsSetNTuple.append( observables[obs] )

# add mass ranges
observables['mass'].setRanges(massRanges)

# build cuts string
if dataSample == 'Summer2011' :
    ntupleCuts = 'runNumber > 87219 && runNumber < 94386' + ( ' && ' if ntupleCuts else '' ) + ntupleCuts

ntupleCuts += ( ' && ' if ntupleCuts else '' ) + selections[selection]

from P2VV.GeneralUtils import readData
dataSets = dict( data = readData( filePath = nTupleFilePath, dataSetName = nTupleName, NTuple = True, observables = obsSetNTuple
                                 , Rename = 'JpsiKK', ntupleCuts = ntupleCuts ) )

if len(KKMassBinBounds) > 2 :
    # create KK mass binning
    from array import array
    KKMassBinsArray = array( 'd', KKMassBinBounds )

    from ROOT import RooBinning
    KKMassBinning = RooBinning( len(KKMassBinBounds) - 1, KKMassBinsArray, 'KKMassBinning' )
    observables['KKMass'].setBinning( KKMassBinning, 'KKMassBinning' )

    # add KK mass split category to data
    from P2VV.RooFitWrappers import BinningCategory
    observables['KKMassCat'] = BinningCategory( 'KKMassCat'
                                               , Observable = observables['KKMass']
                                               , Binning = KKMassBinning
                                               , Fundamental = True
                                               , Data = [ dataSets['data'] ]
                                               , CatTypeName = 'bin'
                                              )


###########################################################################################################################################
## build J/psiKK mass PDFs ##
#############################

# initialize PDF components
nEvents     = dataSets['data'].sumEntries()
nSignal     = nEvents * sigFrac
nBackground = nEvents * ( 1. - sigFrac )

from P2VV.RooFitWrappers import Component
sigMassComps = Component( 'sigMass', [ ], Yield = ( nSignal,     0., nEvents ) )
bkgMassComps = Component( 'bkgMass', [ ], Yield = ( nBackground, 0., nEvents ) )

# build the signal mass PDF
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

signalBMass = SignalBMass( **sigMassArgs )

# build the background mass PDF
bkgMassArgs = dict( Name = 'bkg_m', mass = observables['mass'] )
if bkgMassModel.startswith('linear') :
    from P2VV.Parameterizations.MassPDFs import Linear_Background_Mass as BackgroundBMass
    if bkgMassModel.endswith('Constant') : bkgMassArgs['Constant'] = True
else :
    from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass

backgroundBMass = BackgroundBMass( **bkgMassArgs )

# build mass PDF
from P2VV.RooFitWrappers import buildPdf
sigMassComps += signalBMass.pdf()
bkgMassComps += backgroundBMass.pdf()
massPdf = buildPdf( [ sigMassComps, bkgMassComps ], Observables = [ observables['mass'] ], Name = 'JpsiKKMass' )


###########################################################################################################################################
## fit J/psiKK mass distributions ##
####################################

# determine mass parameters with a fit
print 120 * '='
print 'P2VV - INFO: createB2CCNTuple: fitting with mass PDF'
massFitResult = massPdf.fitTo( dataSets['data'], Save = True, **fitOpts )

from P2VV.Imports import parNames
massFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
massFitResult.covarianceMatrix().Print()
massFitResult.correlationMatrix().Print()

splitCats = [ [ ] ]
if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit'] or len(KKMassBinBounds) > 2 ) :
    # categories for splitting the PDF
    splitCats[0] += [ observables['hlt1ExclB'] ] if selection == 'paper2012'\
                     else [ observables['hlt1ExclB'], observables['hlt2B'] ] if selection == 'timeEffFit' else [ ]
    splitCats[0] += [ observables['KKMassCat'] ] if len(KKMassBinBounds) > 2 else [ ]

    # get mass parameters that are split
    splitParams = [ [ par for par in massPdf.Parameters() if par.getAttribute('Yield') ] ]
    if 'FreeBkg' in SWeightsType and len(KKMassBinBounds) > 2 :
        splitCats.append( [ observables['KKMassCat'] ] )
        splitParams.append( [ par for par in backgroundBMass.parameters() if not par.isConstant() ] )

    # build simultaneous mass PDF
    from P2VV.RooFitWrappers import SimultaneousPdf
    sWeightMassPdf = SimultaneousPdf(  massPdf.GetName() + '_simul'
                                     , MasterPdf       = massPdf
                                     , SplitCategories = splitCats
                                     , SplitParameters = splitParams
                                    )

    # set yields for categories
    splitCat      = sWeightMassPdf.indexCat()
    splitCatIter  = splitCat.typeIterator()
    splitCatState = splitCatIter.Next()
    massPdfPars   = sWeightMassPdf.getVariables()
    from P2VV.GeneralUtils import getSplitPar
    from math import sqrt
    while splitCatState :
        sigYield = getSplitPar( 'N_sigMass', splitCatState.GetName(), massPdfPars )
        bkgYield = getSplitPar( 'N_bkgMass', splitCatState.GetName(), massPdfPars )

        if splitCat.isFundamental() :
            selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
        else :
            splitCat.setLabel( splitCatState.GetName() )
            selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
        nEv    = dataSets['data'].sumEntries()
        nEvBin = dataSets['data'].sumEntries(selStr)

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
        fixedMassPars = [ par for par in sWeightMassPdf.Parameters()\
                          if not ( par.getAttribute('Yield') or par.isConstant() or 'Category' in par.ClassName() ) ]
        #from P2VV.Imports import parValues
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
        from P2VV.Imports import parNames, parValues
        fixedMassPars = [ par for par in sWeightMassPdf.Parameters()\
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

        massNLL = sWeightMassPdf.createNLL( dataSets['data'] )
        simMassFitResult = sWeightMassPdf.fitTo( dataSets['data'], Save = True, **fitOpts )
        simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
        massNLLValNom = massNLL.getVal()
        for par in fixedMassPars :
            par.setConstant(False)
            #par.setVal( fixedMassParVals[par.GetName()][0] )
            #par.setError( fixedMassParVals[par.GetName()][1] )

    # determine mass parameters in each subsample with a fit
    print 120 * '='
    print 'P2VV - INFO: createB2CCNTuple: fitting with simultaneous mass PDF'
    simMassFitResult = sWeightMassPdf.fitTo( dataSets['data'], Save = True, **fitOpts )

    from P2VV.Imports import parValues
    simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    simMassFitResult.covarianceMatrix().Print()
    simMassFitResult.correlationMatrix().Print()

    # hack to do fit with fixed shape parameters first
    if fixedShapeFit :
        massNLLValThis = massNLL.getVal()
        print '  mass NLL values:'
        print '  nominal: %f;  this: %f; 2*DeltaNLL = %f'\
              % ( massNLLValNom, massNLLValThis, 2. * ( massNLLValNom - massNLLValThis ) )

    if SWeightsType.endswith( 'Fixed' ) :
        # free parameters that were fixed for mass fit
        print 'P2VV - INFO: createB2CCNTuple: constant parameters in mass fit:'
        for par in fixedMassPars :
            par.Print()
            par.setConstant(False)

else :
    massPdfPars    = massPdf.getVariables()
    sWeightMassPdf = massPdf


###########################################################################################################################################
## compute S-weights and create signal and background data sets ##
##################################################################

print 120 * '='
print 'P2VV - INFO: createB2CCNTuple: computing S-weights'

# open ROOT file for data sets
from ROOT import TFile

# create sWeigthed data sets
from P2VV.GeneralUtils import SData
SData = SData( Pdf = sWeightMassPdf, Data = dataSets['data'], Name = 'JpsiKK' )
dataSets['SWeightData']    = SData.data()
dataSets['sigSWeightData'] = SData.data('sigMass')
dataSets['bkgSWeightData'] = SData.data('bkgMass')

# print signal/background info to screen
allCats = [  dataSets['data'].get().find( observables['hlt1ExclB'].GetName() )
           , dataSets['data'].get().find( observables['hlt2B'].GetName() )
          ]
if len(KKMassBinBounds) > 2 : allCats.append( dataSets['data'].get().find( observables['KKMassCat'].GetName() ) )
allCats = [ cat for cat in allCats if cat ]

from P2VV.GeneralUtils import printEventYields, printEventYieldsData
for dataSet in dataSets.itervalues() : dataSet.Print()
printEventYields(  ParameterSet        = massPdfPars
                 , YieldNames          = [ 'N_sigMass', 'N_bkgMass' ]
                 , SplittingCategories = [ cat for catList in splitCats for cat in catList ]
                )
printEventYieldsData(  FullDataSet         = dataSets['SWeightData']
                     , WeightedDataSets    = [ dataSets[name] for name in [ 'sigSWeightData', 'bkgSWeightData' ] ]
                     , DataSetNames        = [ 'Signal', 'Background' ]
                     , SplittingCategories = allCats
                    )

# create signal and background data sets with side band ranges
dataSets['sigRangeData'] = dataSets['data'].reduce( Name = 'JpsiKKSigRange', Title = 'JpsiKKSigRange', CutRange = 'Signal'       )
dataSets['bkgRangeData'] = dataSets['data'].reduce( Name = 'JpsiKKBkgRange', Title = 'JpsiKKBkgRange', CutRange = 'LeftSideBand' )
dataSets['bkgRangeData'].append( dataSets['data'].reduce( CutRange = 'RightSideBand' ) )

# create n-tuple containing signal and background weights
dataSets['dataTree'] = dataSets['SWeightData'].buildTree( Name = 'DecayTree', Title = 'DecayTree', RooFitFormat = False )

if plotsFilePath :
    # import plotting tools
    from P2VV.GeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kFullCircle, TPaveText
    from P2VV.GeneralUtils import _P2VVPlotStash

    LHCbLabel = TPaveText( 0.24, 0.81, 0.37, 0.89, 'BRNDC' )
    LHCbLabel.AddText('LHCb')
    LHCbLabel.SetFillColor(0)
    LHCbLabel.SetTextAlign(12)
    LHCbLabel.SetTextSize(0.072)
    LHCbLabel.SetBorderSize(0)
    _P2VVPlotStash.append(LHCbLabel)

    if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit'] or len(KKMassBinBounds) > 2 ) :
        # create projection data set
        indexCat = sWeightMassPdf.indexCat()
        if indexCat.isFundamental() :
            projWDataSet = [ indexCat ]
        else :
            projWDataSet = [ cat for cat in indexCat.getObservables( dataSets['data'] ) ]

        projWData = dict( ProjWData = ( dataSets['data'].reduce( ArgSet = projWDataSet ), False ) )
        print 'P2VV - INFO: createB2CCNTuple: projection data set for mumuKK mass plots:'
        projWData['ProjWData'][0].Print()

    else :
        # don't use projection data set
        projWData = dict()

    # plot mumuKK mass distributions
    massCanvs = [  TCanvas( 'massCanvLog',     'B mass logarithmic scale'  )
                 , TCanvas( 'massCanvSig',     'B mass signal range'       )
                 , TCanvas( 'massCanvLeft',    'B mass left side band'     )
                 , TCanvas( 'massCanvRight',   'B mass right side band'    )
                 , TCanvas( 'massCanvPeakBkg', 'B mass peaking background' )
                ]
    for index, ( pad, frameRange, nBins, plotTitle, plotName, logy, scale, yTitleOffset, markSize, markLineWidth )\
          in enumerate ( zip(  massCanvs
                             , [ '', 'Signal', 'LeftSideBand', 'RightSideBand', 'PeakBkg' ]
                             , numMassBins
                             , [  observables['mass'].GetTitle()
                                , observables['mass'].GetTitle() + ' mass fit - signal'
                                , observables['mass'].GetTitle() + ' mass fit - left side band'
                                , observables['mass'].GetTitle() + ' mass fit - right side band'
                                , observables['mass'].GetTitle() + ' mass fit - peaking background'
                               ]
                             , [  observables['mass'].GetName()
                                , observables['mass'].GetName() + ' fit - signal'
                                , observables['mass'].GetName() + ' fit - left side band'
                                , observables['mass'].GetName() + ' fit - right side band'
                                , observables['mass'].GetName() + ' fit - peaking background'
                               ]
                             , [ True, False, False, False, False ]
                             , [ ( 1.9e2, 1.2e4 ), ( None, None ), ( None, None ), ( None, None ), ( None, None ) ]
                             #, [ ( 1.e3, 2.5e4 ), ( None, None ), ( None, None ), ( None, None ), ( None, None ) ]
                             , [ 1.00, 1.20, 1.15, 1.15, 1.15 ]
                             , [ 0.6,  0.7,  0.8,  0.8,  0.8  ]
                             , [ 2,    3,    3,    3,    3    ]
                       ) ) :
        pad.SetLeftMargin(0.18)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.18)
        pad.SetTopMargin(0.05)

        binWidth = ( observables['mass'].getMax(frameRange) - observables['mass'].getMin(frameRange) ) / float(nBins)
        plot(  pad, observables['mass'], dataSets['data'], sWeightMassPdf, logy = logy, yScale = scale
             , xTitle = 'm(J/#psi K^{+}K^{-}) [MeV/c^{2}]', yTitle = 'Candidates / (%.1f MeV/c^{2})' % binWidth
             , xTitleOffset = 1.10, yTitleOffset = yTitleOffset
             #, plotResidHist = True, normalize = True, symmetrize = True
             , frameOpts  = dict( Range = frameRange, Bins = nBins, Title = plotTitle, Name = plotName )
             , dataOpts   = dict( MarkerStyle = kFullCircle, MarkerSize = markSize, LineWidth = markLineWidth )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 3 )
             , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = 7, LineWidth = 3 )
                             , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = 9, LineWidth = 3 )
                            }
            )
        if index < 2 : LHCbLabel.Draw()

    if SWeightsType.startswith('simultaneous') and ( selection in ['paper2012', 'timeEffFit'] or len(KKMassBinBounds) > 2 ) :
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
                for cat in indexCat.getObservables( dataSets['data'] ) :
                    bins[-1].append( ( cat.GetName(), cat.getIndex(), cat.getLabel() ) )

            pdfs.append( sWeightMassPdf.getPdf( indexCatState.GetName() ) )
            indexCatState = indexCatIter.Next()

        # plot mumuKK mass distributions in KK mass bins
        if   len(bins) <= 4 : nPads = ( 2, 2 )
        elif len(bins) <= 6 : nPads = ( 3, 2 )
        elif len(bins) <= 9 : nPads = ( 3, 3 )
        else :                nPads = ( 4, 3 )
        massCanvs.append( TCanvas( 'massCanvSigBins', 'B mass signal range bins' ) )
        for ( pad, pdf, plotTitle, dataCuts, norm )\
                in zip(  massCanvs[-1].pads( nPads[0], nPads[1] )
                       , pdfs
                       , [ observables['mass'].GetTitle() + ' bin %d - signal' % bin[0][1] for bin in bins ]
                       , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                       , [ dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                           / dataSets['data'].sumEntries() for bin in bins ]
                      ) :
            plot(  pad, observables['mass'], dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                 , frameOpts  = dict( Range = 'Signal', Bins = numMassBins[0], Title = plotTitle )
                 , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                 , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
                                }
                )

        massCanvs.append( TCanvas( 'massCanvLeftBins', 'B mass left side band bins' ) )
        for ( pad, pdf, plotTitle, dataCuts, norm )\
                in zip(  massCanvs[-1].pads( nPads[0], nPads[1] )
                       , pdfs
                       , [ observables['mass'].GetTitle() + ' bin %d - left side band' % bin[0][1] for bin in bins ]
                       , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                       , [ dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                           / dataSets['data'].sumEntries() for bin in bins ]
                      ) :
            plot(  pad, observables['mass'], dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                 , frameOpts  = dict( Range = 'LeftSideBand', Bins = numMassBins[1], Title = plotTitle )
                 , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                 , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
                                }
                )

        massCanvs.append( TCanvas( 'massCanvRightBins', 'B mass right side band bins' ) )
        for ( pad, pdf, plotTitle, dataCuts, norm )\
                in zip(  massCanvs[-1].pads( nPads[0], nPads[1] )
                       , pdfs
                       , [ observables['mass'].GetTitle() + ' bin %d - right side band' % bin[0][1] for bin in bins ]
                       , [ dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) ) for bin in bins ]
                       , [ dataSets['data'].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                           / dataSets['data'].sumEntries() for bin in bins ]
                      ) :
            plot(  pad, observables['mass'], dataSets['data'], pdf#, logy = True, yScale = ( 1., None )
                 , frameOpts  = dict( Range = 'RightSideBand', Bins = numMassBins[2], Title = plotTitle )
                 , dataOpts   = dict( MarkerStyle = 8, MarkerSize = 0.4, **dataCuts               )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                 , components = {  'sig*' : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'bkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
                                }
                )

    for it, canv in enumerate(massCanvs) :
        canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(massCanvs) - 1 else '' ) )


# store data sets in file
dataSetsFile = TFile.Open( dataSetsFilePath, 'RECREATE' )
for data in [ 'SWeightData', 'sigSWeightData', 'bkgSWeightData', 'sigRangeData', 'bkgRangeData', 'dataTree' ] :
    dataSetsFile.Append( dataSets[data] )
dataSetsFile.Write()

dataSetsFile.Close()