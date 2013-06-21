###########################################################################################################################################
## script settings ##
#####################

nTupleFilePath   = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
#nTupleFilePath   = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_2012_20130425_tupleB.root'
nTupleName       = 'DecayTree'
#dataSetsFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets_6KKMassBins_noTagCats.root'
#dataSetsFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets_4KKMassBins_noTagCats.root'
#dataSetsFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/P2VVDataSets_4KKMassBins_freeTagCats.root'
dataSetsFilePath = 'P2VVDataSets_temp.root'
plotsFilePath    = 'plots/P2VVMassPlots2011.ps'

selection        = 'paper2012' # 'HLT1Unbiased' # 'paper2012'
dataSample       = ''
addTaggingObs    = ( 2, 2 ) # ( 0, 0 )
KKMassBinBounds  = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 1008., 1020., 1032. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ]

ntupleCuts = 'sel == 1 && sel_cleantail == 1'\
             ' && muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4. && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
selections = dict(  HLT1Unbiased   = 'hlt1_unbiased_dec==1 && hlt2_biased==1'
                  , HLT1ExclBiased = 'hlt1_excl_biased_dec==1 && hlt2_biased==1'
                  , paper2012      = '(hlt1_biased==1 || hlt1_unbiased_dec==1) && hlt2_biased==1'
                  , timeEffFit     = '(hlt1_biased==1 || hlt1_unbiased_dec==1) && (hlt2_biased==1 || hlt2_unbiased==1)'
                 )

sigFrac          = 0.504
sigMassModel     = ''
cbkgMassModel    = ''
SWeightsType     = 'simultaneousFreeCBkg'
numMassBins      = [ 70, 40, 20, 20, 20 ]
massLogPlotRange = ( 1.9e2, 1.2e4 ) # ( 1.9e2, 1.2e4 ) # ( 8.e2, 2.5e4 )

fitOpts = dict(  NumCPU    = 6
               , Optimize  = 2
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Offset    = True
              )

from math import pi
from ROOT import RooNumber
RooInf  = RooNumber.infinity()
KKMMin  = KKMassBinBounds[0]
KKMMax  = KKMassBinBounds[-1]
obsKeys = [  'mass', 'KKMass', 'mumuMass'
           , 'time', 'timeRes'
           , 'ctk', 'ctl', 'phih'
           #, 'cpsi', 'cttr', 'phitr'
           #, 'wTag', 'tagDec'
           , 'wTagOS', 'tagDecOS'#, 'tagCatOS'
           , 'wTagSS', 'tagDecSS'
           #, 'sel', 'selA', 'selB'
           , 'hlt1ExclB', 'hlt2B', 'hlt2UB'#, 'hlt1B', 'hlt1UB'
          ]

obsDict = dict(  mass      = ( 'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.       )
               , mumuMass  = ( 'mdau1',                'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.       )
               , KKMass    = ( 'mdau2',                'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,  KKMMin,  KKMMax      )
               , time      = ( 'time',                 'Decay time',              'ps',        1.5,    0.3,     14.         )
               , timeRes   = ( 'sigmat',               '#sigma(t)',               'ps',        0.01,   0.0001,  0.12        )
               , ctk       = ( 'helcosthetaK',         'cos(#theta_{K})',         '',          0.,    -1.,     +1.          )
               , ctl       = ( 'helcosthetaL',         'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.          )
               , phih      = ( 'helphi',               '#phi_{h}',                'rad',       0.,    -pi,     +pi          )
               , cpsi      = ( 'trcospsi',             'cos(#psi_{tr})',          '',          0.,    -1.,     +1.          )
               , cttr      = ( 'trcostheta',           'cos(#theta_{tr})',        '',          0.,    -1.,     +1.          )
               , phitr     = ( 'trphi',                '#phi_{tr}',               'rad',       0.,    -pi,     +pi          )
               , wTag      = ( 'tagomega',             'est. wrong-tag prob.',    '',          0.25,   0.,      0.50001     )
               , wTagOS    = ( 'tagomega_os',          'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , wTagSS    = ( 'tagomega_ss',          'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , tagDec    = ( 'tagdecision',          'tag decision',    { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecOS  = ( 'tagdecision_os',       'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecSS  = ( 'tagdecision_ss',       'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagCatOS  = ( 'tagcat_os',            'OS tag category', [ 'unt' ] + [ 'cat%d' % c for c in range(1, 6) ]  )
               , sel       = ( 'sel',                  'selection',       { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selA      = ( 'selA',                 'selection A',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selB      = ( 'selB',                 'selection B',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , hlt1ExclB = ( 'hlt1_excl_biased_dec', 'HLT1 excl. B.',    { 'exclB' : 1, 'notExclB' : 0 }                  )
               , hlt1B     = ( 'hlt1_biased',          'HLT1 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt1UB    = ( 'hlt1_unbiased_dec',    'HLT1 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , hlt2B     = ( 'hlt2_biased',          'HLT2 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt2UB    = ( 'hlt2_unbiased',        'HLT2 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , sel_onecand                = ( 'sel_onecand',                { 'sel' : 1, 'notSel' : 0 }              )
               , sel_one_gl                 = ( 'sel_one_gl',                 { 'sel' : 1, 'notSel' : 0 }              )
               , muPTrChi2                  = ( 'muplus_track_chi2ndof',      'mu+ chi^2/#dof', 1.,    0.,      4.     )
               , muMTrChi2                  = ( 'muminus_track_chi2ndof',     'mu- chi^2/#dof', 1.,    0.,      4.     )
               , KPTrChi2                   = ( 'Kplus_track_chi2ndof',       'K+ chi^2/#dof',  1.,    0.,      4.     )
               , KMTrChi2                   = ( 'Kminus_track_chi2ndof',      'K- chi^2/#dof',  1.,    0.,      4.     )
               , GLsb                       = ( 'GLsb',                       'GLsb',           0.,    0.,      1.     )
               , muplus_PIDmu               = ( 'muplus_PIDmu',               'muplus_PIDmu',   0.,   -RooInf, +RooInf )
               , muminus_PIDmu              = ( 'muminus_PIDmu',              'muminus_PIDmu',  0.,   -RooInf, +RooInf )
               , Kplus_pidK                 = ( 'Kplus_pidK',                 'Kplus_pidK',     0.,    0.,     +RooInf )
               , Kminus_pidK                = ( 'Kminus_pidK',                'Kminus_pidK',    0.,    0.,     +RooInf )
               , muplus_PX                  = ( 'muplus_PX',                  'muplus_PX',      0.,   -RooInf, +RooInf )
               , muplus_PY                  = ( 'muplus_PY',                  'muplus_PY',      0.,   -RooInf, +RooInf )
               , muplus_PZ                  = ( 'muplus_PZ',                  'muplus_PZ',      0.,   -RooInf, +RooInf )
               , muminus_PX                 = ( 'muminus_PX',                 'muminus_PX',     0.,   -RooInf, +RooInf )
               , muminus_PY                 = ( 'muminus_PY',                 'muminus_PY',     0.,   -RooInf, +RooInf )
               , muminus_PZ                 = ( 'muminus_PZ',                 'muminus_PZ',     0.,   -RooInf, +RooInf )
               , Kplus_PX                   = ( 'Kplus_PX',                   'Kplus_PX',       0.,   -RooInf, +RooInf )
               , Kplus_PY                   = ( 'Kplus_PY',                   'Kplus_PY',       0.,   -RooInf, +RooInf )
               , Kplus_PZ                   = ( 'Kplus_PZ',                   'Kplus_PZ',       0.,   -RooInf, +RooInf )
               , Kminus_PX                  = ( 'Kminus_PX',                  'Kminus_PX',      0.,   -RooInf, +RooInf )
               , Kminus_PY                  = ( 'Kminus_PY',                  'Kminus_PY',      0.,   -RooInf, +RooInf )
               , Kminus_PZ                  = ( 'Kminus_PZ',                  'Kminus_PZ',      0.,   -RooInf, +RooInf )
               , B_Pt                       = ( 'B_Pt',                       'B_Pt',           0.,    0.,      RooInf )
               , phi_1020_pt                = ( 'phi_1020_pt',                'phi_1020_pt',    500.,  500.,    RooInf )
               , B_s0_LOKI_CosPolAngle_Dau1 = ( 'B_s0_LOKI_CosPolAngle_Dau1', 'mumu cos(th)',   0.,   -1.,     +1.     )
               , B_s0_IP_OWNPV              = ( 'B_s0_IP_OWNPV',              'B_s0_IP_OWNPV',  0.,   -RooInf, +RooInf )
               , B_s0_IPCHI2_OWNPV          = ( 'B_s0_IPCHI2_OWNPV',          'IP chi2 PV',     0.,   -RooInf, +RooInf )
               , B_s0_MINIPCHI2NEXTBEST     = ( 'B_s0_MINIPCHI2NEXTBEST',     'IP chi2 next',   0.,   -RooInf, +RooInf )
               , B_s0_LOKI_DTF_VCHI2NDOF    = ( 'B_s0_LOKI_DTF_VCHI2NDOF',    'DTF chi2/#dof',  0.,   -RooInf, +RooInf )
               , B_s0_ENDVERTEX_CHI2        = ( 'B_s0_ENDVERTEX_CHI2',        'Bs0 vert chi2',  0.,   0.,       50.    )
               , phi_1020_ENDVERTEX_CHI2    = ( 'phi_1020_ENDVERTEX_CHI2',    'mumu vert chi2', 0.,   0.,       16.    )
               , J_psi_1S_ENDVERTEX_CHI2    = ( 'J_psi_1S_ENDVERTEX_CHI2',    'KK vert chi2',   0.,   0.,       16.    )
              )

massRanges = dict(  LeftSideBand  = ( 5200., 5320. )
                  , Signal        = ( 5320., 5420. )
                  , RightSideBand = ( 5420., 5550. )
                  , PeakBkg       = ( 5390., 5440. )
                 )

if addTaggingObs[0] == 2 :
    tagCatsOS = [  ( 'Untagged', 0, 0.5000001 )
                 , ( 'Tagged',   1, 0.4999999 )
                ]
elif addTaggingObs[0] > 2 :
    tagCatsOS = [  ( 'Untagged', 0, 0.5000001 )
                 , ( 'TagCat1',  1, 0.4999999 )
                 , ( 'TagCat2',  2, 0.38      )
                 , ( 'TagCat3',  3, 0.31      )
                 , ( 'TagCat4',  4, 0.24      )
                 , ( 'TagCat5',  5, 0.17      )
                ]
else :
    tagCatsOS = [ ]

if addTaggingObs[1] == 2 :
    tagCatsSS = [  ( 'Untagged', 0, 0.5000001 )
                 , ( 'Tagged',   1, 0.4999999 )
                ]
elif addTaggingObs[1] > 2 :
    tagCatsSS = [  ( 'Untagged', 0, 0.5000001 )
                 , ( 'TagCat1',  1, 0.4999999 )
                 , ( 'TagCat2',  2, 0.32      )
                 , ( 'TagCat3',  3, 0.25      )
                ]
else :
    tagCatsSS = [ ]


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
    if type( obsDict[obs][2] ) == dict or type( obsDict[obs][2] ) == list :
        observables[obs] = Category( obsDict[obs][0], Title = obsDict[obs][1], Observable = True, States = obsDict[obs][2] )
    else :
        observables[obs] = RealVar( obsDict[obs][0], Title = obsDict[obs][1], Unit = obsDict[obs][2], Observable = True
                                   , Value = obsDict[obs][3], MinMax = ( obsDict[obs][4], obsDict[obs][5] ) )

    obsSetNTuple.append( observables[obs] )

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
sigMassComps  = Component( 'sigMass',  [ ], Yield = ( nSignal,     0., nEvents ) )  # signal
cbkgMassComps = Component( 'cbkgMass', [ ], Yield = ( nBackground, 0., nEvents ) )  # combinatorial background

massComps  = [ sigMassComps, cbkgMassComps ]
yieldNames = [ comp.getYield().GetName() for comp in massComps ]

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

# build the combinatorial background mass PDF
cbkgMassArgs = dict( Name = 'cbkg_m', mass = observables['mass'] )
if cbkgMassModel.startswith('linear') :
    from P2VV.Parameterizations.MassPDFs import Linear_Background_Mass as BackgroundBMass
    if cbkgMassModel.endswith('Constant') : cbkgMassArgs['Constant'] = True
else :
    from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass

backgroundBMass = BackgroundBMass( **cbkgMassArgs )

# build mass PDF
from P2VV.RooFitWrappers import buildPdf
sigMassComps  += signalBMass.pdf()
cbkgMassComps += backgroundBMass.pdf()
massPdf = buildPdf( massComps, Observables = [ observables['mass'] ], Name = 'JpsiKKMass' )


###########################################################################################################################################
## fit J/psiKK mass distributions ##
####################################

# determine mass parameters with a fit
print 120 * '='
print 'P2VV - INFO: createB2CCDataSet: fitting with mass PDF'
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
    if 'FreeCBkg' in SWeightsType and len(KKMassBinBounds) > 2 :
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
        if splitCat.isFundamental() :
            selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
        else :
            splitCat.setLabel( splitCatState.GetName() )
            selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
        nEv    = dataSets['data'].sumEntries()
        nEvBin = dataSets['data'].sumEntries(selStr)

        for yieldVar in [ getSplitPar( name,  splitCatState.GetName(), massPdfPars ) for name in yieldNames ] :
            yieldVar.setVal( yieldVar.getVal() * nEvBin / nEv )
            yieldVar.setError( sqrt( yieldVar.getVal() ) )
            yieldVar.setMin(0.)
            yieldVar.setMax(nEvBin)

        splitCatState = splitCatIter.Next()

    if SWeightsType.endswith( 'Fixed' ) :
        # fix mass shape parameters for fit
        fixedMassPars = [ par for par in sWeightMassPdf.Parameters()\
                          if not ( par.getAttribute('Yield') or par.isConstant() or 'Category' in par.ClassName() ) ]
        #from P2VV.Imports import parValues
        #parVals = {  'm_cbkg_exp_bin0'  : -2.1814e-03
        #           , 'm_cbkg_exp_bin1'  : -4.6151e-03
        #           , 'm_cbkg_exp_bin2'  : -1.4166e-03
        #           , 'm_cbkg_exp_bin3'  : -2.5203e-03
        #           , 'm_cbkg_exp_bin4'  : -1.3963e-03
        #           , 'm_cbkg_exp_bin5'  : -2.0610e-03
        #           , 'm_sig_frac'      :  7.4479e-01
        #           , 'm_sig_mean'      :  5.36809e+03
        #           , 'm_sig_sigma_1'   :  6.1690e+00
        #           , 'm_sig_sigma_sf'  :  2.0769e+00
        #          }
        for par in fixedMassPars :
            #par.setVal( parValues[ par.GetName() ][0]\
            #            + ( 3. * parValues[ par.GetName() ][1] if par.GetName().startswith('m_cbkg_exp') else 0. ) )
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
        #             , 'm_cbkg_exp'      : ( -1.6249e-03,  9.67e-05, -1.       )
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
    print 'P2VV - INFO: createB2CCDataSet: fitting with simultaneous mass PDF'
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
        print 'P2VV - INFO: createB2CCDataSet: constant parameters in mass fit:'
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
print 'P2VV - INFO: createB2CCDataSet: computing S-weights'

# open ROOT file for data sets
from ROOT import TFile

# create sWeigthed data sets
from P2VV.GeneralUtils import SData
SData = SData( Pdf = sWeightMassPdf, Data = dataSets['data'], Name = 'JpsiKK' )
dataSets['SWeightData']     = SData.data()
dataSets['sigSWeightData']  = SData.data( sigMassComps.GetName()  )
dataSets['cbkgSWeightData'] = SData.data( cbkgMassComps.GetName() )

# print signal/background info to screen
allCats = [  dataSets['data'].get().find( obsDict['hlt1ExclB'][0] )
           , dataSets['data'].get().find( obsDict['hlt2B'][0] )
          ]
if len(KKMassBinBounds) > 2 : allCats.append( dataSets['data'].get().find( observables['KKMassCat'].GetName() ) )
allCats = [ cat for cat in allCats if cat ]

from P2VV.GeneralUtils import printEventYields, printEventYieldsData
for dataSet in dataSets.itervalues() : dataSet.Print()
printEventYields(  ParameterSet        = massPdfPars
                 , YieldNames          = yieldNames
                 , SplittingCategories = [ cat for catList in splitCats for cat in catList ]
                )
printEventYieldsData(  FullDataSet         = dataSets['SWeightData']
                     , WeightedDataSets    = [ dataSets[name] for name in [ 'sigSWeightData', 'cbkgSWeightData' ] ]
                     , DataSetNames        = [ 'Signal', 'Combinatorial background' ]
                     , SplittingCategories = allCats
                    )


###########################################################################################################################################
## make J/psiKK mass plots ##
#############################

if plotsFilePath :
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: plotting J/psiKK invariant mass distribution'

    # import plotting tools
    from P2VV.GeneralUtils import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kFullDotLarge, TPaveText
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
        print 'P2VV - INFO: createB2CCDataSet: projection data set for mumuKK mass plots:'
        projWData['ProjWData'][0].Print()

    else :
        # don't use projection data set
        projWData = dict()

    # plot J/psiKK mass distributions
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
                             , [  obsDict['mass'][1]
                                , obsDict['mass'][1] + ' mass fit - signal'
                                , obsDict['mass'][1] + ' mass fit - left side band'
                                , obsDict['mass'][1] + ' mass fit - right side band'
                                , obsDict['mass'][1] + ' mass fit - peaking background'
                               ]
                             , [  obsDict['mass'][0]
                                , obsDict['mass'][0] + ' fit - signal'
                                , obsDict['mass'][0] + ' fit - left side band'
                                , obsDict['mass'][0] + ' fit - right side band'
                                , obsDict['mass'][0] + ' fit - peaking background'
                               ]
                             , [ True, False, False, False, False ]
                             , [ massLogPlotRange, ( None, None ), ( None, None ), ( None, None ), ( None, None ) ]
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
             , plotResidHist = 'E3', normalize = True, symmetrize = True
             , frameOpts  = dict( Range = frameRange, Bins = nBins, Title = plotTitle, Name = plotName )
             , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = markSize, LineWidth = markLineWidth )
             , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 3, Precision = 1.e-4 )
             , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7, LineWidth = 3 )
                             , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9, LineWidth = 3 )
                            }
            )
        if index < 2 :
            pad.cd()
            LHCbLabel.Draw()

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
                 , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.4, **dataCuts  )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm     )
                 , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
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
                 , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.4, **dataCuts        )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                 , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
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
                 , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.4, **dataCuts         )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm      )
                 , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7 )
                                 , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9 )
                                }
                )

    # plot sWeights
    massCanvs += [ TCanvas('sWeightsAll'), TCanvas('sWeightsSignal') ]
    for canv in massCanvs[ -2 : ] :
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)

    allMassVals = [ ]
    sigMassVals = [ ]
    allSWeights = [ ]
    sigSWeights = [ ]
    for obsSet in dataSets['SWeightData'] :
        allMassVals.append( obsSet.getRealValue( obsDict['mass'][0] ) )
        allSWeights.append( obsSet.getRealValue('N_sigMass_sw') )
        if obsSet.find( obsDict['mass'][0] ).inRange('Signal') :
            sigMassVals.append( obsSet.getRealValue( obsDict['mass'][0] ) )
            sigSWeights.append( obsSet.getRealValue('N_sigMass_sw') )

    from array import array
    allMassArr = array( 'd', allMassVals )
    sigMassArr = array( 'd', sigMassVals )
    allSWArr   = array( 'd', allSWeights )
    sigSWArr   = array( 'd', sigSWeights )

    from ROOT import TGraph
    allSWeights = TGraph( len(allMassArr), allMassArr, allSWArr )
    sigSWeights = TGraph( len(sigMassArr), sigMassArr, sigSWArr )
    for graph, range in zip( [ allSWeights, sigSWeights ], [ '', 'Signal' ] ) :
        graph.SetMarkerStyle(kFullDotLarge)
        graph.SetMarkerSize(0.2)
        graph.SetMarkerColor(kBlue)

        #graph.SetMinimum(-0.55)
        #graph.SetMaximum(1.15)
        graph.GetXaxis().SetLimits( observables['mass'].getMin(range), observables['mass'].getMax(range) )
        graph.GetXaxis().SetTitle('m(J/#psi K^{+}K^{-}) [MeV/c^{2}]')
        graph.GetYaxis().SetTitle('sWeight')
        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset(1.0)

    from ROOT import TLine
    allZeroLine = TLine( allSWeights.GetXaxis().GetXmin(), 0., allSWeights.GetXaxis().GetXmax(), 0. )
    sigZeroLine = TLine( sigSWeights.GetXaxis().GetXmin(), 0., sigSWeights.GetXaxis().GetXmax(), 0. )

    massCanvs[-2].cd()
    allSWeights.Draw('AP')
    allZeroLine.Draw()
    massCanvs[-1].cd()
    sigSWeights.Draw('AP')
    sigZeroLine.Draw()

    for it, canv in enumerate(massCanvs) :
        canv.Print( plotsFilePath + ( '(' if it == 0 else ')' if it == len(massCanvs) - 1 else '' ) )


###################################################################################################################################
## add tagging observables to data sets ##
##########################################

if addTaggingObs :
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: building tagging categories'

    # tagging observable names
    wTagOSName   = obsDict['wTagOS'][0]
    wTagSSName   = obsDict['wTagSS'][0]
    tagDecOSName = obsDict['tagDecOS'][0]
    tagDecSSName = obsDict['tagDecSS'][0]

    # get tagging category bins
    from P2VV.Parameterizations.FlavourTagging import getTagCatParamsFromData as getTagParams
    tagBinsOS = getTagParams( dataSets['sigSWeightData'], estWTagName = wTagOSName, tagCats = tagCatsOS, numSigmas = 1., SameSide = False )
    tagBinsSS = getTagParams( dataSets['sigSWeightData'], estWTagName = wTagSSName, tagCats = tagCatsSS, numSigmas = 1., SameSide = True  )

    # add tagging categories to data sets
    from P2VV.GeneralUtils import addTaggingObservables
    for dataKey, data in dataSets.iteritems() :
        if data and not data.get().find('iTagOS') :
            addTaggingObservables( data, 'iTagOS', 'tagCatP2VVOS', tagDecOSName, wTagOSName, tagBinsOS )
            addTaggingObservables( data, 'iTagSS', 'tagCatP2VVSS', tagDecSSName, wTagSSName, tagBinsSS )

    observables['iTagOS']       = Category( ws.put( dataSets['SWeightData'].get().find('iTagOS')       ).GetName() )
    observables['iTagSS']       = Category( ws.put( dataSets['SWeightData'].get().find('iTagSS')       ).GetName() )
    observables['tagCatP2VVOS'] = Category( ws.put( dataSets['SWeightData'].get().find('tagCatP2VVOS') ).GetName() )
    observables['tagCatP2VVSS'] = Category( ws.put( dataSets['SWeightData'].get().find('tagCatP2VVSS') ).GetName() )

    # print tagging categories distributions for signal and background
    from P2VV.RooFitWrappers import ArgSet
    print 'P2VV - INFO: createB2CCDataSet: distribution in opposite side tagging category for signal:'
    dataSets['sigSWeightData'].table(  ArgSet( 'sigOSTagSet',  [ observables['tagCatP2VVOS'], observables['iTagOS'] ] ) ).Print('v')
    print 'P2VV - INFO: createB2CCDataSet: distribution in opposite side tagging category for combinatorial background:'
    dataSets['cbkgSWeightData'].table( ArgSet( 'cbkgOSTagSet', [ observables['tagCatP2VVOS'], observables['iTagOS'] ] ) ).Print('v')
    print 'P2VV - INFO: createB2CCDataSet: distribution in same side tagging category for signal:'
    dataSets['sigSWeightData'].table(  ArgSet( 'sigSSTagSet',  [ observables['tagCatP2VVSS'], observables['iTagSS'] ] ) ).Print('v')
    print 'P2VV - INFO: createB2CCDataSet: distribution in same side tagging category for combinatorial background:'
    dataSets['cbkgSWeightData'].table( ArgSet( 'cbkgSSTagSet', [ observables['tagCatP2VVSS'], observables['iTagSS'] ] ) ).Print('v')


###########################################################################################################################################
## store data sets in ROOT file ##
##################################

# create signal and background data sets with side band ranges
dataSets['sigRangeData']  = dataSets['data'].reduce( Name = 'JpsiKKSigRange',  Title = 'JpsiKKSigRange',  CutRange = 'Signal'       )
dataSets['cbkgRangeData'] = dataSets['data'].reduce( Name = 'JpsiKKCBkgRange', Title = 'JpsiKKCBkgRange', CutRange = 'LeftSideBand' )
dataSets['cbkgRangeData'].append( dataSets['data'].reduce( CutRange = 'RightSideBand' ) )

# create n-tuple containing signal and background weights
dataSets['dataTree'] = dataSets['SWeightData'].buildTree( Name = 'DecayTree', Title = 'DecayTree', RooFitFormat = False )

# save data sets to file
print 120 * '='
print 'P2VV - INFO: createB2CCDataSet: saving data sets to ROOT file %s:' % dataSetsFilePath
dataSetsFile = TFile.Open( dataSetsFilePath, 'RECREATE' )

for data in [ 'SWeightData', 'sigSWeightData', 'cbkgSWeightData', 'sigRangeData', 'cbkgRangeData', 'dataTree' ] :
    dataSetsFile.Append( dataSets[data] )
    print
    if dataSets[data].ClassName() == 'TTree' :
        print 'TTree::%s[%s] = %d entries' % (  dataSets[data].GetName()
                                              , ','.join( br.GetName() for br in dataSets[data].GetListOfBranches() )
                                              , dataSets[data].GetEntries()
                                             )
    else :
        dataSets[data].Print()

dataSetsFile.Write()

dataSetsFile.Close()
