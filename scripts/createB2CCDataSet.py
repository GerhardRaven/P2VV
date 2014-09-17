###########################################################################################################################################
## script settings ##
#####################

#nTupleFilePath  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2011_2012_s20_dv33r6p1_20140213_tupleB_selTrig_TOS.root'
nTupleFilePath  = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/nTupleC_merged_add_20140822.root'
#nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/Bs2JpsiPhi_MC2011_Sim08a_ntupleB_20130909_add.root'
#nTupleFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/MC_Reco14/Bs2JpsiPhi_MC2012_ntupleB_20130904_add.root'
nTupleName       = 'DecayTree'
dataSetsFilePath = 'temp.root' #'/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_6KKMassBins_2TagCats.root'
appendToFile     = False
savedObjects     = [ 'sigSWeight' ] # [ 'main', 'sigSWeight', 'cbkgSWeight' ]
plotsFilePath    = 'temp_plots' #'/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_6KKMassBins_2TagCats_plots'
parFileIn        = 'eventYields6KKBins_HLT2.par' #'eventYields6KKBins.par' #'eventYields6KKBinsNoMC.par' #'eventYields6KKBins.par'
parFileOut       = ''
plotLabelText    = ''

simulation       = False
weightName       = 'pbkgWeight' # 'wMC' # 'pbkgWeight'
runPeriods       = [ 2011, 2012 ]
triggerSel       = 'timeEffFit_tos' # 'paper2012_tos'
dataCuts         = 'nominal2011'
addCuts          = '' # 'hlt2_biased == 1' # 'runPeriod == 2011 && hlt1_excl_biased_dec == 1' # 'wMC > 0.'
dataSample       = '(bkgcat==0 || bkgcat==50)' if simulation else ''
sWeightName      = 'sWeights_ipatia' # 'sigWeight_I2'
addSWeights      = False
addKKMassCat     = True
addTrackMomenta  = False
addTaggingObs    = ( 2, 2 ) # ( 0, 0 )
createRangeData  = False
createNTuple     = False
splitDataSet     = [ ] #[ 'tagCatP2VVOS', 'tagCatP2VVSS' ]
KKMassBinBounds  = [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 1008., 1020., 1032. ]

sigMassModel     = 'Ipatia2'
cbkgMassModel    = ''
SWeightsType     = 'simultaneousFreeCBkg'

eventFracs       = [  dict( N_sigMass = 0.44, N_cbkgMass = None )
                    #, dict( N_sigMass = 0.11,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.70,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.64,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.18,  N_cbkgMass = None )
                   ]

constSplitPars = [ ]
constSplitCats = [ ]
constSplitVals = { }
if sigMassModel == 'Ipatia2' :
    #constSplitCats = [ 'hlt1ExclB', 'KKMassCat' ]
    constSplitCats = [ 'KKMassCat' ]
    constSplitPars = [ 'm_sig_alpha_1', 'm_sig_alpha_2', 'm_sig_n_1', 'm_sig_n_2', 'm_sig_lambda', 'm_sig_zeta' ]
    constSplitVals = {  ( 'bin0', ) : [ 3.1265598805552086, 2.4478555551857504, 0.01,                0.4626575371693939,  -1.940905888712173, 0.0010225614571429675  ]
                      , ( 'bin1', ) : [ 2.235667390313107,  2.296350006990024,  1.1725685774531813,  1.696945683056572,   -2.590902278024503, 0.2967019395179116     ]
                      , ( 'bin2', ) : [ 2.8232366360188683, 2.570677068540935,  2.046626071513164,   2.6043347350995187,  -2.76448170618875,  0.00019473384424262274 ]
                      , ( 'bin3', ) : [ 3.1064191306717484, 2.096424131956603,  1.69939291075719,    2.7310345435871,     -2.903782724519466, 0.0007918143518637066  ]
                      , ( 'bin4', ) : [ 2.5929196863701067, 2.1195375110354044, 1.2511053016584157,  1.432994401017839,   -2.626462929014679, 0.19046609538221904    ]
                      , ( 'bin5', ) : [ 3.0236236294665466, 3.033593188175826,  0.39183360830444464, 0.43646654437298116, -2.216215361029956, 0.00019353774759423015 ]
                     }

numMassBins      = [ 140, 56, 40, 44, 24 ]
massLogPlotRange = ( 5.e2, 2.e4 ) # ( 1.e2, 6.e3 ) # ( 4.e2, 1.3e4 )

fitOpts = dict(  NumCPU    = 8
               , Optimize  = 1
               , Timer     = True
#               , Verbose   = True
#               , Minos     = True
#               , Hesse     = False
               , Minimizer = 'Minuit2'
               , Strategy  = 1
               , Offset    = True
              )

from math import pi
from ROOT import RooNumber
RooInf  = RooNumber.infinity()
KKMMin  = KKMassBinBounds[0]
KKMMax  = KKMassBinBounds[-1]

obsKeys = [ # 'sWeights_ipatia',
            # 'sWeight_orig'
             'wMC'
           #, 'hlt2_prescale', 'polarity', 'nPVCat', 'BpTCat'
           , 'runPeriod'
           , 'mass', 'KKMass', 'mumuMass'
           , 'time', 'timeRes'
           #, 'truetime'
           , 'ctk', 'ctl', 'phih'
           #, 'cpsi', 'cttr', 'phitr'
           #, 'wTag', 'tagDec'
           , 'wTagOS'#, 'tagDecOS'
           , 'wTagSS'#, 'tagDecSS'
           #, 'tagCatOS', 'tagCatSS'
           #, 'sel', 'selClTail'#, 'selA', 'selB'
           , 'hlt1ExclB', 'hlt2B', 'hlt2UB'#, 'hlt1B', 'hlt1UB'
           #, 'trigDecUnb', 'trigDecExclB'
           #, 'B_s0_ENDVERTEX_X', 'B_s0_ENDVERTEX_Y', 'B_s0_ENDVERTEX_Z'
           #, 'B_P', 'B_Pt', 'B_eta', 'B_phi'
           #, 'Kplus_PX', 'Kplus_PY', 'Kplus_PZ', 'Kplus_LOKI_ETA'
           #, 'Kminus_PX', 'Kminus_PY', 'Kminus_PZ', 'Kminus_LOKI_ETA'
           #, 'muplus_PX', 'muplus_PY', 'muplus_PZ', 'muplus_LOKI_ETA'
           #, 'muminus_PX', 'muminus_PY', 'muminus_PZ', 'muminus_LOKI_ETA'
          ]
if simulation :
    obsKeys += [ 'truetime', 'bkgcat' ]
if addTrackMomenta :
    obsKeys += [ '%s_P%s' % ( part, comp ) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'X', 'Y', 'Z' ) ]

obsDict = dict(  runPeriod = ( 'runPeriod',            'run period', dict( [ ( 'p%d' % period, period ) for period in runPeriods ] ) )
               , wMC       = ( 'wMC',                  'pbkgWeight',              '',          0.,     -RooInf, +RooInf     )
               , firstData = ( 'firstData',            'first data',      { 'first' : +1, 'later' : 0 }                     )
               , hlt2_prescale = ( 'hlt2_prescale',    'HLT2 prescale',   { 'presc' : +1, 'noPresc' : 0, 'notDef1' : -1, 'notDef2' : -2 } )
               , polarity  = ( 'polarity',             'magnet polarity', { 'up' : +1, 'down' : -1, 'noPol' : 0 }           )
               , nPVCat    = ( 'nPVCat',               'num. PVs',        dict( [ ( 'c%d' % cat, cat ) for cat in range(11) ] )  )
               , BpTCat    = ( 'BpTCat',               'B p_T cat',       dict( [ ( 'c%d' % cat, cat ) for cat in range(5) ] )   )
               , mass      = ( 'mass',                 'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.       )
               , mumuMass  = ( 'mdau1',                'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.       )
               , KKMass    = ( 'mdau2',                'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,  KKMMin,  KKMMax      )
               , time      = ( 'time',                 'Decay time',              'ps',        1.5,    0.3,     14.         )
               , truetime  = ( 'truetime',             'true time',               'ps',        0.,     0.,      30.         )
               , timeRes   = ( 'sigmat',               '#sigma(t)',               'ps',        0.01,   0.0001,  0.12        )
               , ctk       = ( 'helcosthetaK',         'cos(#theta_{K})',         '',          0.,    -1.,     +1.          )
               , ctl       = ( 'helcosthetaL',         'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.          )
               , phih      = ( 'helphi',               '#phi_{h}',                'rad',       0.,    -pi,     +pi          )
               , cpsi      = ( 'trcospsi',             'cos(#psi_{tr})',          '',          0.,    -1.,     +1.          )
               , cttr      = ( 'trcostheta',           'cos(#theta_{tr})',        '',          0.,    -1.,     +1.          )
               , phitr     = ( 'trphi',                '#phi_{tr}',               'rad',       0.,    -pi,     +pi          )
               , wTag      = ( 'tagomega',             'est. wrong-tag prob.',    '',          0.25,   0.,      0.50001     )
               #, wTagOS    = ( 'tagomega_os',          'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               #, wTagSS    = ( 'tagomega_ss',          'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , wTagOS    = ( 'tagomega_os_cb',       'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               #, wTagSS    = ( 'tagomega_ss_cb',       'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               #, wTagOS    = ( 'tagomega_os_nn',       'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , wTagSS    = ( 'tagomega_ss_nn',       'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , tagDec    = ( 'tagdecision',          'tag decision',    { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               #, tagDecOS  = ( 'tagdecision_os',       'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               #, tagDecSS  = ( 'tagdecision_ss',       'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecOS  = ( 'tagdecision_os_cb',    'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               #, tagDecOS  = ( 'iTagOS',               'OS tag decision', { 'B' : +1, 'Bbar' : -1 }                          )
               #, tagDecSS  = ( 'tagdecision_ss_cb',    'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               #, tagDecOS  = ( 'tagdecision_os_nn',    'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecSS  = ( 'tagdecision_ss_nn',    'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               #, tagDecSS  = ( 'iTagSS',               'SS tag decision', { 'B' : +1, 'Bbar' : -1 }                          )
               #, tagCatOS  = ( 'tagcat_os',            'OS tag category', [ 'unt' ] + [ 'cat%d' % c for c in range(1, 6) ]  )
               , tagCatOS  = ( 'tagCatP2VVOS',         'OS tag category', { 'Untagged' : 0, 'Tagged' : 1 }                  )
               , tagCatSS  = ( 'tagCatP2VVSS',         'SS tag category', { 'Untagged' : 0, 'Tagged' : 1 }                  )
               , bkgcat    = ( 'bkgcat',               'background cat',  { 'signal' : 0, 'lowMass' : 50 }                  )
               , trueid    = ( 'trueid',               'true B ID',       { 'B' : 531, 'Bbar' : 531,  'noMatch' : 0 }       )
               , sel       = ( 'sel',                  'selection',       { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selClTail = ( 'sel_cleantail',        'clean tail sel.', { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selA      = ( 'selA',                 'selection A',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selB      = ( 'selB',                 'selection B',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , hlt1ExclB = ( 'hlt1_excl_biased',     'HLT1 excl. B.',   { 'exclB' : 1, 'notExclB' : 0 }                   )
               #, hlt1ExclB = ( 'hlt1_excl_biased_dec', 'HLT1 excl. B.',    { 'exclB' : 1, 'notExclB' : 0 }                  )
               , hlt1B     = ( 'hlt1_biased',          'HLT1 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt1UB    = ( 'hlt1_unbiased',        'HLT1 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               #, hlt1UB    = ( 'hlt1_unbiased_dec',    'HLT1 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , hlt2B     = ( 'hlt2_biased',          'HLT2 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt2UB    = ( 'hlt2_unbiased',        'HLT2 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
               , trigDecUnb                 = ( 'triggerDecisionUnbiased',   'trigDecUnb',   { 'UB'    : 1, 'notUB'    : 0 } )
               , trigDecExclB               = ( 'triggerDecisionBiasedExcl', 'trigDecExclB', { 'exclB' : 1, 'notExclB' : 0 } )
               , sel_onecand                = ( 'sel_onecand',                { 'sel' : 1, 'notSel' : 0 }              )
               , sel_one_gl                 = ( 'sel_one_gl',                 { 'sel' : 1, 'notSel' : 0 }              )
               , muPTrChi2                  = ( 'muplus_track_chi2ndof',      'mu+ chi^2/#dof', 1.,    0.,      4.     )
               , muMTrChi2                  = ( 'muminus_track_chi2ndof',     'mu- chi^2/#dof', 1.,    0.,      4.     )
               , KPTrChi2                   = ( 'Kplus_track_chi2ndof',       'K+ chi^2/#dof',  1.,    0.,      4.     )
               , KMTrChi2                   = ( 'Kminus_track_chi2ndof',      'K- chi^2/#dof',  1.,    0.,      4.     )
               , GLsb                       = ( 'GLsb',                       'GLsb',           0.,    0.,      1.     )
               , B_s0_ENDVERTEX_X           = ( 'B_s0_ENDVERTEX_X',           'B_s0_ENDVERTEX_X', '', 0., -RooInf, +RooInf )
               , B_s0_ENDVERTEX_Y           = ( 'B_s0_ENDVERTEX_Y',           'B_s0_ENDVERTEX_Y', '', 0., -RooInf, +RooInf )
               , B_s0_ENDVERTEX_Z           = ( 'B_s0_ENDVERTEX_Z',           'B_s0_ENDVERTEX_Z', '', 0., -RooInf, +RooInf )
               , muplus_PIDmu               = ( 'muplus_PIDmu',               'muplus_PIDmu',   0.,   -RooInf, +RooInf )
               , muminus_PIDmu              = ( 'muminus_PIDmu',              'muminus_PIDmu',  0.,   -RooInf, +RooInf )
               , Kplus_pidK                 = ( 'Kplus_pidK',                 'Kplus_pidK',     0.,    0.,     +RooInf )
               , Kminus_pidK                = ( 'Kminus_pidK',                'Kminus_pidK',    0.,    0.,     +RooInf )
               , muplus_PX                  = ( 'muplus_PX',                  'muplus_PX',        'MeV/c', 0.,   -RooInf, +RooInf )
               , muplus_PY                  = ( 'muplus_PY',                  'muplus_PY',        'MeV/c', 0.,   -RooInf, +RooInf )
               , muplus_PZ                  = ( 'muplus_PZ',                  'muplus_PZ',        'MeV/c', 0.,   -RooInf, +RooInf )
               , muplus_LOKI_ETA            = ( 'muplus_LOKI_ETA',            'muplus_LOKI_ETA',  '',      0.,   -RooInf, +RooInf )
               , muminus_PX                 = ( 'muminus_PX',                 'muminus_PX',       'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_PY                 = ( 'muminus_PY',                 'muminus_PY',       'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_PZ                 = ( 'muminus_PZ',                 'muminus_PZ',       'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_LOKI_ETA           = ( 'muminus_LOKI_ETA',           'muminus_LOKI_ETA', '',      0.,   -RooInf, +RooInf )
               , Kplus_PX                   = ( 'Kplus_PX',                   'Kplus_PX',         'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_PY                   = ( 'Kplus_PY',                   'Kplus_PY',         'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_PZ                   = ( 'Kplus_PZ',                   'Kplus_PZ',         'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_LOKI_ETA             = ( 'Kplus_LOKI_ETA',             'Kplus_LOKI_ETA',   '',      0.,   -RooInf, +RooInf )
               , Kminus_PX                  = ( 'Kminus_PX',                  'Kminus_PX',        'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_PY                  = ( 'Kminus_PY',                  'Kminus_PY',        'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_PZ                  = ( 'Kminus_PZ',                  'Kminus_PZ',        'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_LOKI_ETA            = ( 'Kminus_LOKI_ETA',            'Kminus_LOKI_ETA',  '',      0.,   -RooInf, +RooInf )
               , B_Pt                       = ( 'B_Pt',                       'B_Pt',             'MeV/c', 0.,    0.,      RooInf )
               , B_P                        = ( 'B_P',                        'B_P',              'MeV/c', 0.,    0.,      RooInf )
               , B_eta                      = ( 'B_eta',                      'B_eta',            '',      0.,   -RooInf, +RooInf )
               , B_phi                      = ( 'B_phi',                      'B_phi',            'rad',   0.,    -pi,     +pi    )
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

massRanges = dict(  LeftSideBand  = ( 5200., 5300. )
                  , Signal        = ( 5300., 5440. )
                  , RightSideBand = ( 5440., 5550. )
                  , PeakBkg       = ( 5380., 5440. )
                 )

if addTaggingObs and addTaggingObs[0] == 1 :
    tagCatsOS = [ ( 'Untagged', 0, 0.500000001 ) ]
elif addTaggingObs and addTaggingObs[0] == 2 :
    tagCatsOS = [  ( 'Untagged', 0, 0.500000001 )
                 , ( 'Tagged',   1, 0.499999999 )
                ]
elif addTaggingObs and addTaggingObs[0] == 4 :
    tagCatsOS = [  ( 'Untagged', 0, 0.500000001 )
                 , ( 'TagCat1',  1, 0.499999999 )
                 , ( 'TagCat2',  2, 0.34        )
                 , ( 'TagCat3',  3, 0.22        )
                ]
elif addTaggingObs and addTaggingObs[0] != 0 :
    raise AssertionError( 'P2VV - ERROR: createB2CCDataSet: unknown number of opposite-side tagging categories: %s' % addTaggingObs[0] )
else :
    tagCatsOS = [ ]

if addTaggingObs and addTaggingObs[1] == 1 :
    tagCatsSS = [ ( 'Untagged', 0, 0.500000001 ) ]
elif addTaggingObs and addTaggingObs[1] == 2 :
    tagCatsSS = [  ( 'Untagged', 0, 0.500000001 )
                 , ( 'Tagged',   1, 0.499999999 )
                ]
elif addTaggingObs and addTaggingObs[1] == 4 :
    tagCatsSS = [  ( 'Untagged', 0, 0.500000001 )
                 , ( 'TagCat1',  1, 0.499999999 )
                 , ( 'TagCat2',  2, 0.40        )
                 , ( 'TagCat3',  3, 0.30        )
                ]
elif addTaggingObs and addTaggingObs[1] != 0 :
    raise AssertionError( 'P2VV - ERROR: createB2CCDataSet: unknown number of same-side tagging categories: %s' % addTaggingObs[1] )
else :
    tagCatsSS = [ ]

from P2VV.Imports import extConstraintValues
extConstraintValues.setVal( 'DM',      ( 17.768,  0.024   ) )
extConstraintValues.setVal( 'P0OS',    (  0.3791, 0.0044, 0.3791 ) )
extConstraintValues.setVal( 'P1OS',    (  1.00,   0.035  ) )
extConstraintValues.setVal( 'DelP0OS', (  0.0140, 0.0012 ) )
extConstraintValues.setVal( 'DelP1OS', (  0.066,  0.012  ) )
extConstraintValues.setVal( 'P0SS',    (  0.445,  0.005, 0.445 ) )
extConstraintValues.setVal( 'P1SS',    (  1.00,   0.09   ) )
extConstraintValues.setVal( 'DelP0SS', ( -0.016,  0.002  ) )
extConstraintValues.setVal( 'DelP1SS', (  0.007,  0.019  ) )

plotsFilePath = plotsFilePath.split('.')[0]

###########################################################################################################################################
## read data ##
###############

from math import sqrt
from P2VV.Load import RooFitOutput, LHCbStyle

# create list of required observables
reqObsList = [ 'index', 'mass', 'KKMass', 'tagDecOS', 'tagDecSS', 'wTagOS', 'wTagSS', 'hlt1ExclB', 'hlt2B' ]
reqObsList += [ weightName ] if weightName else [ ]
reqObsList += ['runPeriod'] if runPeriods else [ ]
reqObsList += [ '%s_P%s' % ( part, comp ) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'X', 'Y', 'Z' ) ]\
              if addTrackMomenta else [ ]

# create workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# create observables
from P2VV.RooFitWrappers import RealVar, Category
observables  = { }
observables['index'] = RealVar( 'index', Title = 'event index', Observable = True, Value = 0., MinMax = ( -RooInf, +RooInf ) )
obsSetPreDS  = [ observables['index'] ]
if weightName :
    observables[weightName] = RealVar( weightName, Title = 'event weight', Observable = True, Value = 0., MinMax = ( -RooInf, +RooInf ) )
    obsSetPreDS += [ observables[weightName] ]
obsSetNTuple = [ ]
for obs in obsKeys + reqObsList :
    if obs in observables : continue

    if obs.startswith('sigWeight') or obs.startswith('sWeight') :
        observables[obs] = RealVar( obs, Title = 'signal sWeight', Observable = True, Value = 1. )

    elif type( obsDict[obs][2] ) == dict or type( obsDict[obs][2] ) == list :
        observables[obs] = Category( obsDict[obs][0], Title = obsDict[obs][1], Observable = True, States = obsDict[obs][2] )

    else :
        observables[obs] = RealVar( obsDict[obs][0], Title = obsDict[obs][1], Unit = obsDict[obs][2], Observable = True
                                   , Value = obsDict[obs][3], MinMax = ( obsDict[obs][4], obsDict[obs][5] ) )

    if obs in obsKeys    : obsSetNTuple.append( observables[obs] )
    if obs in reqObsList : obsSetPreDS.append(  observables[obs] )

# add mass ranges
observables['mass'].setRanges(massRanges)

# build cuts string
from P2VV.Imports import cutSelStrings
ntupleCuts = cutSelStrings[dataCuts]
if dataSample == 'Summer2011' :
    ntupleCuts = 'runNumber >= 87219 && runNumber <= 94386' + ( ' && ' if ntupleCuts else '' ) + ntupleCuts
elif dataSample and type(dataSample) == str : ntupleCuts += ( ' && ' if ntupleCuts else '' ) + dataSample

from P2VV.Imports import triggerSelStrings
trigSelStr = triggerSelStrings[triggerSel]
if trigSelStr : ntupleCuts += ( ' && ' if ntupleCuts else '' ) + trigSelStr

if addCuts : ntupleCuts += ( ' && ' if ntupleCuts else '' ) + addCuts

# read n-tuple file(s)
print 'P2VV - INFO: createB2CCDataSet: reading n-tuple "%s" from file "%s"' % ( nTupleName, nTupleFilePath )
print 'P2VV - INFO: createB2CCDataSet: selection cuts:\n    %s' % ntupleCuts
from ROOT import TFile, TTree
dataTreeFile = TFile.Open(nTupleFilePath)
assert dataTreeFile, 'P2VV - ERROR: createB2CCDataSet: could not open file "%s"' % nTupleFilePath
dataTree = dataTreeFile.Get(nTupleName)
assert dataTree, 'P2VV - ERROR: createB2CCDataSet: could not locate tree "%s" in file "%s"' % ( nTupleName, nTupleFilePath )
from ROOT import gROOT
gROOT.cd('PyROOT:/')

# create data set from n-tuple
dataSets = dict( pre = ( dataTree.buildDataSet( Observables = obsSetPreDS, Name = 'JpsiKK', Title = 'JpsiKK', Cuts = ntupleCuts
                        , IndexName = 'index', WeightName = weightName if weightName else '' ), [ ] ) )
print 'P2VV - INFO: createB2CCDataSet: data set from n-tuple:\n' + ' ' * 13,
dataSets['pre'][0].Print()
print

#selDataFilePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/Bs2JpsiPhi_2011_2012_s20_dv33r6p1_20140213_tupleB_selTrig_TOS.root'
#selDataFile = TFile.Open( selDataFilePath, 'RECREATE' )
#selDataTree = dataTree.CopyTree(ntupleCuts)
#print 'P2VV - INFO: selected number of entries in tree = %d' % selDataTree.GetEntries()
#from ROOT import TObject
#selDataFile.Write( selDataFilePath, TObject.kOverwrite )
#selDataFile.Close()
#assert False

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
                                           , Data = [ dataSets['pre'][0] ]
                                           , CatTypeName = 'bin'
                                          )
if addKKMassCat : obsSetNTuple.append( observables['KKMassCat'] )


###########################################################################################################################################
## build J/psiKK mass PDFs ##
#############################

if not simulation :
    # initialize PDF components
    sigFrac     = eventFracs[0]['N_sigMass'] if eventFracs[0]['N_sigMass'] else 1. - sum(val for val in eventFracs[0].itervalues() if val)
    nEvents     = dataSets['pre'][0].sumEntries()
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

    elif sigMassModel.startswith('Ipatia2') :
        from P2VV.Parameterizations.MassPDFs import Ipatia2_Signal_Mass as SignalBMass
        sigMassArgs['m_sig_lambda']  = -2.5
        sigMassArgs['m_sig_zeta']    =  0.01
        sigMassArgs['m_sig_alpha_1'] =  3.0
        sigMassArgs['m_sig_alpha_2'] =  2.5
        sigMassArgs['m_sig_n_1']     =  1.5
        sigMassArgs['m_sig_n_1']     =  2.0

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

splitPars = [ ]
if not simulation :
    # determine mass parameters with a fit
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: fitting with mass PDF'
    massFitResult = massPdf.fitTo( dataSets['pre'][0], Save = True, SumW2Error = False, **fitOpts )

    from P2VV.Imports import parNames
    massFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
    massFitResult.covarianceMatrix().Print()
    massFitResult.correlationMatrix().Print()

    splitCats = { }
    if SWeightsType.startswith('simultaneous') :
        # get splitting categories and parameters
        yieldNames  = [ par.GetName() for par in massPdf.Parameters() if par.getAttribute('Yield') ]
        allParNames = [ par.GetName() for par in massPdf.Parameters() if not par.getAttribute('Yield') if not par.isConstant() ]
        sigParNames = [ par.GetName() for par in signalBMass.parameters() if not par.isConstant() ]
        bkgParNames = [ par.GetName() for par in backgroundBMass.parameters() if not par.isConstant() ]
        if runPeriods :
            # split parameters for run period
            splitCats[ observables['runPeriod'].GetName() ] = set( yieldNames + bkgParNames )
            #if sigMassModel.startswith('Ipatia2') :
            #    for parName in sigParNames :
            #        if 'mean' in parName or 'sigma' in parName : splitCats[ observables['runPeriod'].GetName() ].add(parName)

        if len(KKMassBinBounds) > 1 :
            # split parameters for KK-mass category
            splitCats[ observables['KKMassCat'].GetName() ] = set(yieldNames)
            if sigMassModel.startswith('Ipatia2') :
                for parName in sigParNames : splitCats[ observables['KKMassCat'].GetName() ].add(parName)
            #if 'FreeCBkg' in SWeightsType :
            #    for parName in bkgParNames : splitCats[ observables['KKMassCat'].GetName() ].add(parName)

        # split parameters for HLT1 category
        #splitCats[ observables['hlt1ExclB'].GetName() ] = set(yieldNames)
        #if sigMassModel.startswith('Ipatia2') :
        #    for parName in sigParNames : splitCats[ observables['hlt1ExclB'].GetName() ].add(parName)
        #if 'FreeCBkg' in SWeightsType :
        #    for parName in bkgParNames : splitCats[ observables['hlt1ExclB'].GetName() ].add(parName)

        # split parameters for HLT2 category
        splitCats[ observables['hlt2B'].GetName() ] = set(yieldNames)
        if sigMassModel.startswith('Ipatia2') :
            for parName in sigParNames :
                if 'sigma' in parName : splitCats[ observables['hlt2B'].GetName() ].add(parName)
        if 'FreeCBkg' in SWeightsType :
            for parName in bkgParNames : splitCats[ observables['hlt2B'].GetName() ].add(parName)

        for cat in constSplitCats :
            # split specified constant shape parameters
            catName = observables[cat].GetName()
            if catName not in splitCats : splitCats[catName] = set()
            for parName in constSplitPars : splitCats[catName].add(parName)

        splitParsDict = { }
        pdfVars = massPdf.getVariables()
        for cat, params in splitCats.iteritems() :
            assert ws.cat(cat), 'P2VV - ERROR: createB2CCDataSet: category "%s" not in workspace' % cat
            for par in params :
                assert ws.var(par), 'P2VV - ERROR: createB2CCDataSet: no variable "%s" in workspace' % par
                assert pdfVars.find(par), 'P2VV - ERROR: createB2CCDataSet: variable "%s" not in PDF' % par
                if ws[par] not in splitParsDict :
                    splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
                else :
                    splitParsDict[ ws[par] ].add( ws[cat] )

        # create list of split parameters and categories
        from P2VV.Utilities.General import createSplitParsList
        splitPars = createSplitParsList(splitParsDict)

        # build simultaneous mass PDF
        print 'P2VV - INFO: createB2CCDataSet: building simultaneous PDF "%s":' % ( massPdf.GetName() + '_simul' )
        print 13 * ' ' + 'split parameters:'
        for it, pars in enumerate(splitPars) :
            print 13 * ' ' + '%2d: pars: [ %s ]' % ( it, ', '.join( par.GetName() for par in pars[0] ) )
            print 13 * ' ' + '    cats: [ %s ]' % ', '.join( cat.GetName() for cat in pars[1] )
        from P2VV.RooFitWrappers import SimultaneousPdf
        sWeightMassPdf = SimultaneousPdf(  massPdf.GetName() + '_simul'
                                         , MasterPdf       = massPdf
                                         , SplitCategories = [ pars[1] for pars in splitPars ]
                                         , SplitParameters = [ pars[0] for pars in splitPars ]
                                        )

        # set yields and constant splitted parameters
        splitCat      = sWeightMassPdf.indexCat()
        inputCats     = [ splitCat ] if splitCat.isFundamental() else splitCat.inputCatList()
        splitCatIter  = splitCat.typeIterator()
        splitCatState = splitCatIter.Next()
        massPdfPars   = sWeightMassPdf.getVariables()
        from P2VV.Utilities.General import getSplitPar
        while splitCatState :
            splitCat.setIndex( splitCatState.getVal() )
            KKMassState = -1
            if splitCat.isFundamental() :
                selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
                if splitCat.GetName() == observables['KKMassCat'].GetName() : KKMassState = splitCatState.getVal()
            else :
                splitCat.setLabel( splitCatState.GetName() )
                selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
                if observables['KKMassCat'] in [ par for pars in splitPars for par in pars[1] ] :
                    for cat in splitCat.inputCatList() :
                        if cat.GetName() == observables['KKMassCat'].GetName() : KKMassState = cat.getIndex()

            nEv    = dataSets['pre'][0].sumEntries()
            nEvBin = dataSets['pre'][0].sumEntries(selStr)
            for yName in yieldNames :
                stateName = splitCatState.GetName()
                yieldVar = getSplitPar( yName, stateName, massPdfPars )

                if len(eventFracs) > 1 and KKMassState >= 0 :
                    sigFrac = eventFracs[ KKMassState + 1 ][yName] if eventFracs[ KKMassState + 1 ][yName]\
                              else 1. - sum( val for val in eventFracs[ KKMassState + 1 ].itervalues() if val )
                    yieldVar.setVal( sigFrac * nEvBin )
                else :
                    yieldVar.setVal( yieldVar.getVal() * nEvBin / nEv )

                yieldVar.setError( sqrt( yieldVar.getVal() ) )
                yieldVar.setMin(0.)
                yieldVar.setMax(nEvBin)

            catStates = tuple( [ observables[catName].getLabel() for catName in constSplitCats ] )
            stateName = catStates[0] if len(catStates) == 1 else '{' + ';'.join( state for state in catStates ) + '}'
            for it, parName in enumerate(constSplitPars) :
                par = getSplitPar( parName, stateName, massPdfPars )
                par.setVal( constSplitVals[catStates][it] )

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

            massNLL = sWeightMassPdf.createNLL( dataSets['pre'][0] )
            simMassFitResult = sWeightMassPdf.fitTo( dataSets['pre'][0], Save = True, SumW2Error = False, **fitOpts )
            simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
            massNLLValNom = massNLL.getVal()
            for par in fixedMassPars :
                par.setConstant(False)
                #par.setVal( fixedMassParVals[par.GetName()][0] )
                #par.setError( fixedMassParVals[par.GetName()][1] )

        if parFileIn or parFileOut :
            # create PDF configuration object
            from P2VV.Parameterizations.FullPDFs import PdfConfiguration as PdfConfig
            pdfConfig = PdfConfig()

        print 120 * '='
        if parFileIn :
            # read parameters from file
            pdfConfig.readParametersFromFile( filePath = parFileIn )
            pdfConfig.setParametersInPdf(sWeightMassPdf)
        print 'P2VV - INFO: createB2CCDataSet: parameters in simultaneous mass PDF:'
        sWeightMassPdf.getParameters( dataSets['pre'][0] ).Print('v')

        # determine mass parameters in each subsample with a fit
        print 'P2VV - INFO: createB2CCDataSet: fitting with simultaneous mass PDF'
        #from P2VV.Utilities.DataHandling import correctWeights
        #fitData = correctWeights( dataSets['pre'][0], [ 'runPeriod', 'hlt2_biased', 'KKMassCat' ] )
        fitData = dataSets['pre'][0]
        simMassFitResult = sWeightMassPdf.fitTo( fitData, Save = True, SumW2Error = False, **fitOpts )
        #del fitData

        if parFileOut :
            # write parameters to file
            pdfConfig.getParametersFromPdf( sWeightMassPdf, dataSets['pre'][0] )
            pdfConfig.writeParametersToFile( filePath = parFileOut )

        simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames )
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
## make J/psiKK mass plots ##
#############################

if not simulation and plotsFilePath :
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: plotting J/psiKK invariant mass distribution'

    # import plotting tools
    from P2VV.Utilities.Plotting import plot
    from ROOT import TCanvas, kBlue, kRed, kGreen, kFullDotLarge, TPaveText, TLatex
    from P2VV.Utilities.Plotting import _P2VVPlotStash

    if plotLabelText :
        plotLabel = TLatex()
        plotLabel.SetTextAlign(32)
        plotLabel.SetTextSize(0.072)
        _P2VVPlotStash.append(plotLabel)

    binLabel = TLatex()
    binLabel.SetTextAlign(12)
    binLabel.SetTextSize(0.05)
    _P2VVPlotStash.append(binLabel)

    if SWeightsType.startswith('simultaneous') :
        # create projection data set
        indexCat = sWeightMassPdf.indexCat()
        if indexCat.isFundamental() :
            projWDataSet = [ indexCat ]
        else :
            projWDataSet = [ cat for cat in indexCat.inputCatList() ]

        projWData = dict( ProjWData = ( dataSets['pre'][0].reduce( ArgSet = projWDataSet ), False ) )
        print 'P2VV - INFO: createB2CCDataSet: projection data set for J/psiKK mass plots:'
        projWData['ProjWData'][0].Print()

    else :
        # don't use projection data set
        projWData = dict()

    # plot J/psiKK mass distributions
    massPlots = [ ]
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
                                , obsDict['mass'][0] + '_signal'
                                , obsDict['mass'][0] + '_left'
                                , obsDict['mass'][0] + '_right'
                                , obsDict['mass'][0] + '_peakBkg'
                               ]
                             , [ True, False, False, False, False ]
                             , [ massLogPlotRange, ( 0., None ), ( 0., None ), ( 0., None ), ( 0., None ) ]
                             #, [ False, False, False, False, False ]
                             #, [ ( None, None ), ( None, None ), ( None, None ), ( None, None ), ( None, None ) ]
                             , [ 1.00, 1.30, 1.15, 1.15, 1.15 ]
                             , [ 0.6,  0.7,  0.8,  0.8,  0.8  ]
                             , [ 2,    3,    3,    3,    3    ]
                       ) ) :
        pad.SetLeftMargin(0.18)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.18)
        pad.SetTopMargin(0.05)

        binWidth = ( observables['mass'].getMax(frameRange) - observables['mass'].getMin(frameRange) ) / float(nBins)
        massPlots.append( plot(  pad, observables['mass'], dataSets['pre'][0], sWeightMassPdf, logy = logy, yScale = scale
                               , xTitle = 'm(J/#psi K^{+}K^{-}) [MeV/c^{2}]', yTitle = 'Candidates / (%.1f MeV/c^{2})' % binWidth
                               , xTitleOffset = 1.10, yTitleOffset = yTitleOffset
                               #, plotResidHist = 'E3', normalize = True, symmetrize = True
                               , plotResidHist = 'BX', normalize = True, symmetrize = True
                               , frameOpts  = dict( Range = frameRange, Bins = nBins, Title = plotTitle, Name = plotName )
                               , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = markSize, LineWidth = markLineWidth )
                               , pdfOpts    = dict( list( projWData.items() ), LineColor = kBlue, LineWidth = 3, Precision = 1.e-4 )
                               , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7, LineWidth = 3 )
                                               , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9, LineWidth = 3 )
                                              }
                              )
        )
        if index < 2 :
            pad.cd()
            if plotLabelText : plotLabel.DrawLatexNDC( 0.87, 0.88, plotLabelText )

    if SWeightsType.startswith('simultaneous') and len(KKMassBinBounds) > 2 :
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
                for cat in indexCat.inputCatList() :
                    bins[-1].append( ( cat.GetName(), cat.getIndex(), cat.getLabel() ) )

            pdfs.append( sWeightMassPdf.getPdf( indexCatState.GetName() ) )
            indexCatState = indexCatIter.Next()

        # plot distributions in simultaneous bins
        for ( bin, pdf ) in zip(  bins, pdfs ) :
            massCanvs.append( TCanvas( 'massCanvSigBins%d' % bin[0][1], 'B mass signal range bins' ) )
            plotTitle = observables['mass'].GetTitle() + ' bin %s - signal' % bin[0][2]
            dataCutStr = dict( Cut = ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )
            norm = dataSets['pre'][0].sumEntries( ' && '.join( '%s==%d' % ( c[0], c[1] ) for c in bin[ 1 : ] ) )\
                   / dataSets['pre'][0].sumEntries()

            pad.SetLeftMargin(0.18)
            pad.SetRightMargin(0.05)
            pad.SetBottomMargin(0.18)
            pad.SetTopMargin(0.05)

            binWidth = ( observables['mass'].getMax('Signal') - observables['mass'].getMin('Signal') ) / float(numMassBins[1])
            plot(  massCanvs[-1], observables['mass'], dataSets['pre'][0], pdf
                 , xTitle = 'm(J/#psi K^{+}K^{-}) [MeV/c^{2}]', yTitle = 'Candidates / (%.1f MeV/c^{2})' % binWidth
                 , xTitleOffset = 1.10, yTitleOffset = 1.20
                 #, plotResidHist = 'E3', normalize = True, symmetrize = True
                 , plotResidHist = 'BX', normalize = True, symmetrize = True
                 , frameOpts  = dict( Range = 'Signal', Bins = numMassBins[1], Title = plotTitle )
                 , dataOpts   = dict( MarkerStyle = kFullDotLarge, MarkerSize = 0.4, **dataCutStr )
                 , pdfOpts    = dict( LineColor = kBlue, LineWidth = 3, Normalization = norm )
                 , components = {  'sig*'  : dict( LineColor = kRed,       LineStyle = 7, LineWidth = 3 )
                                 , 'cbkg*' : dict( LineColor = kGreen + 3, LineStyle = 9, LineWidth = 3 )
                                }
                )

            massCanvs[-1].cd()
            binLabel.DrawLatex( 0.18, 0.91, bin[0][2].strip('{}') )

    for it, canv in enumerate(massCanvs) :
        canv.Print( plotsFilePath + '_mass.ps' + ( '(' if it == 0 else ')' if it == len(massCanvs) - 1 else '' ) )

    # write mass plots to ROOT file
    from ROOT import TFile, TObject
    plotsFile = TFile.Open( plotsFilePath + '.root', 'RECREATE' )
    for massPlot in massPlots :
        plotsFile.Append( massPlot[0] )
        plotsFile.Append( massPlot[1] )
    plotsFile.Write( plotsFilePath + '.root', TObject.kOverwrite )
    plotsFile.Close()


###########################################################################################################################################
## compute sWeights ##
######################

weightVars = [ ]
if not simulation :
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: computing sWeights'

    # compute sWeights
    from P2VV.Utilities.SWeights import SData
    sWeightsBuilder = SData( Pdf = sWeightMassPdf, Data = dataSets['pre'][0], Name = 'JpsiKK' )
    dataSets['preS'] = ( sWeightsBuilder.data(), [ ] )

    # import sWeight variables
    for wName in [ 'N_sigMass_sw', 'N_cbkgMass_sw' ] :
        ws.put( dataSets['preS'][0].get().find(wName) )
        weightVars.append( RealVar(wName) )

else :
    dataSets['preS'] = ( dataSets['pre'][0], [ ] )


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
    tagBinsOS = getTagParams( dataSets['preS'][0], estWTagName = wTagOSName, tagCats = tagCatsOS, numSigmas = 2., SameSide = False
                             , WeightVarName = '' if simulation else weightVars[0].GetName() )
    tagBinsSS = getTagParams( dataSets['preS'][0], estWTagName = wTagSSName, tagCats = tagCatsSS, numSigmas = 2., SameSide = True
                             , WeightVarName = '' if simulation else weightVars[0].GetName() )

    # add tagging categories to data sets
    from P2VV.Utilities.DataHandling import addTaggingObservables
    addTaggingObservables( dataSets['preS'][0], 'iTagOS', 'tagCatP2VVOS', tagDecOSName, wTagOSName, tagBinsOS )
    addTaggingObservables( dataSets['preS'][0], 'iTagSS', 'tagCatP2VVSS', tagDecSSName, wTagSSName, tagBinsSS )

    observables['iTagOS']       = Category( ws.put( dataSets['preS'][0].get().find('iTagOS')       ).GetName() )
    observables['iTagSS']       = Category( ws.put( dataSets['preS'][0].get().find('iTagSS')       ).GetName() )
    observables['tagCatP2VVOS'] = Category( ws.put( dataSets['preS'][0].get().find('tagCatP2VVOS') ).GetName() )
    observables['tagCatP2VVSS'] = Category( ws.put( dataSets['preS'][0].get().find('tagCatP2VVSS') ).GetName() )
    obsSetNTuple += [ observables['iTagOS'], observables['iTagSS'], observables['tagCatP2VVOS'], observables['tagCatP2VVSS'] ]

    # add OS+SS tagged/untagged category
    from P2VV.Utilities.DataHandling import addGlobalTagCat
    addGlobalTagCat( dataSets['preS'][0], observables['tagCatP2VVOS'].GetName(), observables['tagCatP2VVSS'].GetName() )
    observables['tagCatP2VV'] = Category( ws.put( dataSets['preS'][0].get().find('tagCatP2VV') ).GetName() )
    obsSetNTuple += [ observables['tagCatP2VV'] ]


###################################################################################################################################
## add track and B momenta to data set ##
#########################################

if addTrackMomenta :
    # add K+, K-, mu+ and mu- momentum magnitudes to data set
    from P2VV.Utilities.DataHandling import addTrackMomenta
    addTrackMomenta( dataSets['preS'][0] )

    observables['Kplus_P']   = RealVar( ws.put( dataSets['preS'][0].get().find('Kplus_P')   ).GetName() )
    observables['Kminus_P']  = RealVar( ws.put( dataSets['preS'][0].get().find('Kminus_P')  ).GetName() )
    observables['muplus_P']  = RealVar( ws.put( dataSets['preS'][0].get().find('muplus_P')  ).GetName() )
    observables['muminus_P'] = RealVar( ws.put( dataSets['preS'][0].get().find('muminus_P') ).GetName() )
    obsSetNTuple += [ observables['Kplus_P'], observables['Kminus_P'], observables['muplus_P'], observables['muminus_P'] ]


###########################################################################################################################################
## create signal and background data sets and print event yields ##
###################################################################

print 120 * '='
print 'P2VV - INFO: createB2CCDataSet: creating data set'

# create data sets with final columns
mainDSList = [ ]
samples = [ ]
obsSetMain = obsSetNTuple + [ var for var in weightVars ]
if splitDataSet :
    # get category states for splitting data set
    states = tuple( [ [ ( cat, observables[cat].GetName(), state.getVal() ) for state in observables[cat] ] for cat in splitDataSet ] )
    from itertools import product
    for state in product(*states) :
        nameStr = '_'.join(  '%s%d'   % ( cState[0], cState[2] ) for cState in state ).replace( '-', 'm' )
        selStr  = '&&'.join( '%s==%d' % ( cState[1], cState[2] ) for cState in state )
        samples.append( ( nameStr, selStr ) )

    # create splitted data sets
    from ROOT import RooDataSet, RooArgSet
    obsSet = RooArgSet( obs._var for obs in obsSetMain )
    mainDS = RooDataSet( 'JpsiKK', 'JpsiKK', obsSet )
    for sample in samples :
        preDS = dataSets['preS'][0].reduce( Cut = sample[1] )
        mainDSList.append( dataTree.buildDataSet( Observables = obsSetMain, Name = 'JpsiKK_' + sample[0], Title = 'JpsiKK'
                                                 , IndexName = 'index', OrigDataSet = preDS ) )
        mainDS.append( mainDSList[-1] )

else :
    # create data set without splitting
    mainDS = dataTree.buildDataSet( Observables = obsSetMain, Name = 'JpsiKK', Title = 'JpsiKK', IndexName = 'index'
                                , OrigDataSet = dataSets['preS'][0] )

dataSets['main'] = ( mainDS, mainDSList )
dataSets.pop('pre')
dataSets.pop('preS')
dataTreeFile.Close()
gROOT.cd('PyROOT:/')

print 'P2VV - INFO: createB2CCDataSet: produced data set:\n' + ' ' * 13,
dataSets['main'][0].Print()
print

if not simulation :
    # create weighted data sets
    from ROOT import RooArgSet
    obsSets = dict( [ ( var.GetName(), RooArgSet( dataSets['main'][0].get() ) ) for var in weightVars ] )
    if not addSWeights :
        for varName, obsSet in obsSets.iteritems() :
            for var in weightVars :
                if var.GetName() == varName : continue
                obsSet.remove( obsSet.find( var.GetName() ) )

    from ROOT import RooDataSet
    dataSets['sigSWeight']  = ( RooDataSet( 'JpsiKK_sigSWeight',  'JpsiKK_sigSWeight', obsSets['N_sigMass_sw']
                                           , Import = dataSets['main'][0], WeightVar = ( weightVars[0].GetName(), True ) ), [ ] )
    dataSets['cbkgSWeight'] = ( RooDataSet( 'JpsiKK_cbkgSWeight', 'JpsiKK_cbkgSWeight', obsSets['N_cbkgMass_sw']
                                           , Import = dataSets['main'][0], WeightVar = ( weightVars[1].GetName(), True ) ), [ ] )
    for sample, data in zip( samples, dataSets['main'][1] ) :
        dataSets['sigSWeight'][1].append(  RooDataSet( 'JpsiKK_sigSWeight_'  + sample[0], 'JpsiKK_sigSWeight', obsSets['N_sigMass_sw']
                                                      , Import = data, WeightVar = ( weightVars[0].GetName(), True ) ) )
        dataSets['cbkgSWeight'][1].append( RooDataSet( 'JpsiKK_cbkgSWeight_' + sample[0], 'JpsiKK_cbkgSWeight', obsSets['N_cbkgMass_sw']
                                                      , Import = data, WeightVar = ( weightVars[1].GetName(), True ) ) )

    print 'P2VV - INFO: createB2CCDataSet: signal data set:\n' + ' ' * 13,
    dataSets['sigSWeight'][0].Print()
    print
    print 'P2VV - INFO: createB2CCDataSet: combinatorial background data set:\n' + ' ' * 13,
    dataSets['cbkgSWeight'][0].Print()
    print


    # print yields
    print 'P2VV - INFO: createB2CCDataSet: event yields:'
    allCats = [ dataSets['main'][0].get().find( obsDict['runPeriod'][0] ) ] if runPeriods else [ ]
    allCats += [  dataSets['main'][0].get().find( obsDict['hlt1ExclB'][0] )
                , dataSets['main'][0].get().find( obsDict['hlt2B'][0] )
               ]
    if len(KKMassBinBounds) > 2 : allCats.append( dataSets['main'][0].get().find( observables['KKMassCat'].GetName() ) )
    allCats = [ cat for cat in allCats if cat ]

    from P2VV.Utilities.DataHandling import printEventYields, printEventYieldsData
    printEventYields(  ParameterSet        = massPdfPars
                     , YieldNames          = yieldNames
                     , SplittingCategories = [ cat for pars in splitPars for cat in pars[1] ]
                    )
    printEventYieldsData(  FullDataSet         = dataSets['main'][0]
                         , WeightedDataSets    = [ dataSets[name][0] for name in [ 'sigSWeight', 'cbkgSWeight' ] ]
                         , DataSetNames        = [ 'Signal', 'Combinatorial background' ]
                         , SplittingCategories = allCats
                        )

if addTaggingObs :
    # print tagging categories distributions for signal and background
    sigDataName = 'main' if simulation else 'sigSWeight'
    from P2VV.RooFitWrappers import ArgSet
    print 'P2VV - INFO: createB2CCDataSet: distribution in opposite side tagging category for signal:'
    dataSets[sigDataName][0].table(  ArgSet( 'sigOSTagSet',  [ observables['tagCatP2VVOS'], observables['iTagOS'] ] ) ).Print('v')
    if not simulation and 'cbkgSWeight' in savedObjects :
        print 'P2VV - INFO: createB2CCDataSet: distribution in opposite side tagging category for combinatorial background:'
        dataSets['cbkgSWeight'][0].table( ArgSet( 'cbkgOSTagSet', [ observables['tagCatP2VVOS'], observables['iTagOS'] ] ) ).Print('v')
    print 'P2VV - INFO: createB2CCDataSet: distribution in same side tagging category for signal:'
    dataSets[sigDataName][0].table(  ArgSet( 'sigSSTagSet',  [ observables['tagCatP2VVSS'], observables['iTagSS'] ] ) ).Print('v')
    if not simulation and 'cbkgSWeight' in savedObjects :
        print 'P2VV - INFO: createB2CCDataSet: distribution in same side tagging category for combinatorial background:'
        dataSets['cbkgSWeight'][0].table( ArgSet( 'cbkgSSTagSet', [ observables['tagCatP2VVSS'], observables['iTagSS'] ] ) ).Print('v')


###########################################################################################################################################
## make sWeight plots ##
#############################

if not simulation and plotsFilePath :
    print 120 * '='
    print 'P2VV - INFO: createB2CCDataSet: plotting sWeights versus J/psiKK mass'

    # plot sWeights
    sWeightCanvs = [ TCanvas('sWeightsAll'), TCanvas('sWeightsSignal') ]
    for canv in sWeightCanvs[ -2 : ] :
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)

    allMassVals = [ ]
    sigMassVals = [ ]
    allSWeights = [ ]
    sigSWeights = [ ]
    for obsSet in dataSets['main'][0] :
        allMassVals.append( obsSet.getRealValue( obsDict['mass'][0] ) )
        allSWeights.append( obsSet.getRealValue( weightVars[0].GetName() ) )
        if obsSet.find( obsDict['mass'][0] ).inRange('Signal') :
            sigMassVals.append( obsSet.getRealValue( obsDict['mass'][0] ) )
            sigSWeights.append( obsSet.getRealValue( weightVars[0].GetName() ) )

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
        graph.GetYaxis().SetTitle('signal weight')
        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset(1.0)

    from ROOT import TLine
    allZeroLine = TLine( allSWeights.GetXaxis().GetXmin(), 0., allSWeights.GetXaxis().GetXmax(), 0. )
    sigZeroLine = TLine( sigSWeights.GetXaxis().GetXmin(), 0., sigSWeights.GetXaxis().GetXmax(), 0. )

    sWeightCanvs[-2].cd()
    allSWeights.Draw('AP')
    allZeroLine.Draw()
    sWeightCanvs[-1].cd()
    sigSWeights.Draw('AP')
    sigZeroLine.Draw()

    for it, canv in enumerate(sWeightCanvs) :
        canv.Print( plotsFilePath + '_sWeights.ps' + ('(' if it == 0 else ')' if it == len(sWeightCanvs) - 1 else '') )


###########################################################################################################################################
## store data sets in ROOT file ##
##################################

if createRangeData :
    # create signal and background data sets with side band ranges
    dataSets['sigRange']  = dataSets['main'][0].reduce( Name = 'JpsiKK_sigRange',  Title = 'JpsiKK_sigRange',  CutRange = 'Signal'       )
    dataSets['cbkgRange'] = dataSets['main'][0].reduce( Name = 'JpsiKK_cbkgRange', Title = 'JpsiKK_cbkgRange', CutRange = 'LeftSideBand' )
    dataSets['cbkgRange'].append( dataSets['main'][0].reduce( CutRange = 'RightSideBand' ) )
    savedObjects += [ 'sigRange', 'cbkgRange' ]

if createNTuple and not simulation :
    # create n-tuple containing signal and background weights
    dataSets['sigTree'] = dataSets['sigSWeight'][0].buildTree( Name = nTupleName, Title = nTupleName, WeightName = sWeightName
                                                              , RooFitFormat = False )
    savedObjects += [ 'sigTree' ]
elif createNTuple :
    # create n-tuple
    dataSets['mainTree'] = dataSets['main'][0].buildTree( Name = nTupleName, Title = nTupleName, RooFitFormat = False )
    savedObjects += [ 'mainTree' ]

# save data sets to file
print 120 * '='
print 'P2VV - INFO: createB2CCDataSet: saving data sets to ROOT file %s:' % dataSetsFilePath
from ROOT import TFile
dataSetsFile = TFile.Open( dataSetsFilePath, 'UPDATE' if appendToFile else 'RECREATE' )

for dataKey, data in dataSets.iteritems() :
    print
    if hasattr( data, 'ClassName' ) and data.ClassName() == 'TTree' :
        # data in TTree
        print 'TTree::%s[%s] = %d entries'\
              % (  data.GetName(), ','.join( br.GetName() for br in data.GetListOfBranches() ), data.GetEntries() )
        if dataKey in savedObjects : dataSetsFile.Append(data)

    else :
        # data in RooDataSet(s)
        if not data[1] :
            data[0].Print()
            if dataKey in savedObjects :
                if data[0].numEntries() > 0 :
                    dataSetsFile.Append(data[0])
                else :
                    print 'P2VV - INFO: createB2CCDataSet: no entries in "%s", not saving data set' % data[0].GetName()
        for subData in data[1] :
            print
            subData.Print()
            if dataKey in savedObjects :
                if subData.numEntries() > 0 :
                    dataSetsFile.Append(subData)
                else :
                    print 'P2VV - INFO: createB2CCDataSet: no entries in "%s", not saving data set' % subData.GetName()

from ROOT import TObject
dataSetsFile.Write( dataSetsFilePath, TObject.kOverwrite )
dataSetsFile.Close()
