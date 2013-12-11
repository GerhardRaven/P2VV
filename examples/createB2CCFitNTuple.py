###########################################################################################################################################
## script settings ##
#####################

nTupleFilePath  = 'data/Reco14/2011_2012_dv33r6p1_s20_201309_tupleB_add.root'
nTupleName       = 'DecayTree'
dataSetsFilePath = 'fitNTuple_2011_2012_Reco14.root'
appendToFile     = False

runPeriods       = [ 2011, 2012 ]
triggerSel       = 'paper2012' # 'noSelection' # 'paper2012' # 'HLT1Unbiased' # 'paper2012'
dataCuts         = 'nominal2011' # 'noSelection' # 'nominal2011'
sWeightName      = 'sigWeight_4KKBins'
splitDataSet     = [ ]
addKKMassCat     = True
KKMassBinBounds  = [ 990., 1020. - 12., 1020., 1020. + 12., 1050. ] # [ 990., 1020. - 12., 1020. - 4., 1020., 1020. + 4., 1020. + 12., 1050. ] # [ 1008., 1020., 1032. ]
eventFracs       = [  dict( N_sigMass = 0.504, N_cbkgMass = None )
                    #, dict( N_sigMass = 0.11,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.70,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.64,  N_cbkgMass = None )
                    #, dict( N_sigMass = 0.18,  N_cbkgMass = None )
                   ]

fitOpts = dict(  NumCPU    = 2
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

obsKeys = [  'runPeriod'
           , 'mass', 'KKMass', 'mumuMass'
           , 'time', 'timeRes'
           , 'ctk', 'ctl', 'phih'
           , 'cpsi', 'cttr', 'phitr'
           , 'wTag', 'tagDec'
           , 'wTagOS', 'tagDecOS'
           , 'wTagSS', 'tagDecSS'
           , 'sel'#, 'selA', 'selB'
           , 'hlt1ExclB', 'hlt2B', 'hlt2UB'#, 'hlt1B', 'hlt1UB'
           , 'trigDecUnb', 'trigDecExclB'
          ]

obsDict = dict(  runPeriod = ( 'runPeriod',            'run period', dict( [ ( 'p%d' % period, period ) for period in runPeriods ] ) )
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
               , wTagOS    = ( 'tagomega_os',          'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , wTagSS    = ( 'tagomega_ss',          'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001     )
               , tagDec    = ( 'tagdecision',          'tag decision',    { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecOS  = ( 'tagdecision_os',       'OS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagDecSS  = ( 'tagdecision_ss',       'SS tag decision', { 'B' : +1, 'Bbar' : -1, 'Untagged' : 0 }         )
               , tagCatOS  = ( 'tagcat_os',            'OS tag category', [ 'unt' ] + [ 'cat%d' % c for c in range(1, 6) ]  )
               , bkgcat    = ( 'bkgcat',               'background cat',  { 'signal' : 0, 'lowMass' : 50 }                  )
               , sel       = ( 'sel',                  'selection',       { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selA      = ( 'selA',                 'selection A',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , selB      = ( 'selB',                 'selection B',     { 'sel'   : 1, 'notSel'   : 0 }                   )
               , hlt1ExclB = ( 'hlt1_excl_biased_dec', 'HLT1 excl. B.',    { 'exclB' : 1, 'notExclB' : 0 }                  )
               , hlt1B     = ( 'hlt1_biased',          'HLT1 B.',          { 'B'     : 1, 'notB'     : 0 }                  )
               , hlt1UB    = ( 'hlt1_unbiased_dec',    'HLT1 UB.',         { 'UB'    : 1, 'notUB'    : 0 }                  )
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
               , muplus_PIDmu               = ( 'muplus_PIDmu',               'muplus_PIDmu',   0.,   -RooInf, +RooInf )
               , muminus_PIDmu              = ( 'muminus_PIDmu',              'muminus_PIDmu',  0.,   -RooInf, +RooInf )
               , Kplus_pidK                 = ( 'Kplus_pidK',                 'Kplus_pidK',     0.,    0.,     +RooInf )
               , Kminus_pidK                = ( 'Kminus_pidK',                'Kminus_pidK',    0.,    0.,     +RooInf )
               , muplus_PX                  = ( 'muplus_PX',                  'muplus_PX',     'MeV/c', 0.,   -RooInf, +RooInf )
               , muplus_PY                  = ( 'muplus_PY',                  'muplus_PY',      'MeV/c', 0.,   -RooInf, +RooInf )
               , muplus_PZ                  = ( 'muplus_PZ',                  'muplus_PZ',      'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_PX                 = ( 'muminus_PX',                 'muminus_PX',     'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_PY                 = ( 'muminus_PY',                 'muminus_PY',     'MeV/c', 0.,   -RooInf, +RooInf )
               , muminus_PZ                 = ( 'muminus_PZ',                 'muminus_PZ',     'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_PX                   = ( 'Kplus_PX',                   'Kplus_PX',       'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_PY                   = ( 'Kplus_PY',                   'Kplus_PY',       'MeV/c', 0.,   -RooInf, +RooInf )
               , Kplus_PZ                   = ( 'Kplus_PZ',                   'Kplus_PZ',       'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_PX                  = ( 'Kminus_PX',                  'Kminus_PX',      'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_PY                  = ( 'Kminus_PY',                  'Kminus_PY',      'MeV/c', 0.,   -RooInf, +RooInf )
               , Kminus_PZ                  = ( 'Kminus_PZ',                  'Kminus_PZ',      'MeV/c', 0.,   -RooInf, +RooInf )
               , B_Pt                       = ( 'B_Pt',                       'B_Pt',           'MeV/c', 0.,    0.,      RooInf )
               , B_P                        = ( 'B_P',                        'B_P',            'MeV/c', 0.,    0.,      RooInf )
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


###########################################################################################################################################
## read data ##
###############

from P2VV.Load import RooFitOutput, LHCbStyle

# create list of required observables
reqObsList = [ 'index', 'mass', 'KKMass', 'tagDecOS', 'tagDecSS', 'wTagOS', 'wTagSS' ]
reqObsList += ['runPeriod'] if runPeriods else [ ]
reqObsList += ['hlt1ExclB'] if triggerSel == 'paper2012' else [ 'hlt1ExclB', 'hlt2B' ] if triggerSel == 'timeEffFit' else [ ]

# create workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# create observables
from P2VV.RooFitWrappers import RealVar, Category
observables  = { }
observables['index'] = RealVar( 'index', Title = 'event index', Observable = True, Value = 0., MinMax = ( -RooInf, +RooInf ) )
obsSetPreDS  = [ observables['index'] ]
obsSetNTuple = [ ]
for obs in obsKeys + reqObsList :
    if obs in observables : continue

    if obs.startswith('sigWeight') :
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
from P2VV.Imports import cutSelStrings, triggerSelStrings
ntupleCuts = cutSelStrings[dataCuts]
ntupleCuts += ( ' && ' if ntupleCuts else '' ) + triggerSelStrings[triggerSel]

# read n-tuple file(s)
print 'P2VV - INFO: createB2CCFitNTuple: reading n-tuple "%s" from file "%s"' % ( nTupleName, nTupleFilePath )
print 'P2VV - INFO: createB2CCFitNTuple: selection cuts:\n    %s' % ntupleCuts
from ROOT import TFile, TTree
dataTreeFile = TFile.Open(nTupleFilePath)
assert dataTreeFile, 'P2VV - ERROR: createB2CCFitNTuple: could not open file "%s"' % nTupleFilePath
dataTree = dataTreeFile.Get(nTupleName)
assert dataTree, 'P2VV - ERROR: createB2CCFitNTuple: could not locate tree "%s" in file "%s"' % ( nTupleName, nTupleFilePath )

# create data set from n-tuple
dataSets = dict( pre = ( dataTree.buildDataSet( Observables = obsSetPreDS, Name = 'JpsiKK', Title = 'JpsiKK', Cuts = ntupleCuts
                        , IndexName = 'index' ), [ ] ) )
print 'P2VV - INFO: createB2CCFitNTuple: data set from n-tuple:\n' + ' ' * 13,
dataSets['pre'][0].Print()
print

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
from P2VV.Parameterizations.MassPDFs import LP2011_Signal_Mass as SignalBMass
signalBMass = SignalBMass( **sigMassArgs )

# build the combinatorial background mass PDF
cbkgMassArgs = dict( Name = 'cbkg_m', mass = observables['mass'] )
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
print 'P2VV - INFO: createB2CCFitNTuple: fitting with mass PDF'
massFitResult = massPdf.fitTo( dataSets['pre'][0], Save = True, **fitOpts )

from P2VV.Imports import parNames
massFitResult.PrintSpecial( text = True, ParNames = parNames )

if runPeriods or triggerSel in ['paper2012', 'timeEffFit'] or len(KKMassBinBounds) > 2 :
    # get splitting categories and parameters
    splitCats = { }
    yieldNames = [ par.GetName() for par in massPdf.Parameters() if par.getAttribute('Yield') ]
    if runPeriods :
        splitCats[observables['runPeriod'].GetName()] = set([ par.GetName() for par in massPdf.Parameters() if not par.isConstant() ])
    if len(KKMassBinBounds) > 2 :
        splitCats[ observables['KKMassCat'].GetName() ] = set(yieldNames)
        for par in backgroundBMass.parameters() :
            if not par.isConstant() : splitCats[ observables['KKMassCat'].GetName() ].add( par.GetName() )
    if triggerSel == 'paper2012' or triggerSel == 'timeEffFit' :
        splitCats[ observables['hlt1ExclB'].GetName() ] = set(yieldNames)
    if triggerSel == 'timeEffFit' :
        splitCats[ observables['hlt2B'].GetName() ] = set(yieldNames)

    splitParsDict = { }
    pdfVars = massPdf.getVariables()
    for cat, params in splitCats.iteritems() :
        assert ws.cat(cat), 'P2VV - ERROR: createB2CCFitNTuple: category "%s" not in workspace' % cat
        for par in params :
            assert ws.var(par), 'P2VV - ERROR: createB2CCFitNTuple: no variable "%s" in workspace' % par
            assert pdfVars.find(par), 'P2VV - ERROR: createB2CCFitNTuple: variable "%s" not in PDF' % par
            if ws[par] not in splitParsDict :
                splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
            else :
                splitParsDict[ ws[par] ].add( ws[cat] )

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

    # build simultaneous mass PDF
    print 'P2VV - INFO: createB2CCFitNTuple: building simultaneous PDF "%s":' % ( massPdf.GetName() + '_simul' )
    print 13 * ' ' + 'splitting categories: [ %s ]' % ' ], [ '.join(', '.join(cat.GetName() for cat in cats) for cats in splitCats)
    print 13 * ' ' + 'split parameters:     [ %s ]' % ' ], [ '.join(', '.join(par.GetName() for par in pars) for pars in splitPars)
    from P2VV.RooFitWrappers import SimultaneousPdf
    sWeightMassPdf = SimultaneousPdf(  massPdf.GetName() + '_simul'
                                     , MasterPdf       = massPdf
                                     , SplitCategories = splitCats
                                     , SplitParameters = splitPars
                                    )

    # set yields for categories
    splitCat      = sWeightMassPdf.indexCat()
    splitCatIter  = splitCat.typeIterator()
    splitCatState = splitCatIter.Next()
    massPdfPars   = sWeightMassPdf.getVariables()
    from P2VV.Utilities.General import getSplitPar
    from math import sqrt
    while splitCatState :
        KKMassState = -1
        if splitCat.isFundamental() :
            selStr = '!(%s-%d)' % ( splitCat.GetName(), splitCatState.getVal() )
            if splitCat.GetName() == observables['KKMassCat'].GetName() : KKMassState = splitCatState.getVal()
        else :
            splitCat.setLabel( splitCatState.GetName() )
            selStr = ' && '.join( '!(%s-%d)' % ( cat.GetName(), cat.getIndex() ) for cat in splitCat.inputCatList() )
            if observables['KKMassCat'] in splitCats[0] :
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

        splitCatState = splitCatIter.Next()

    # determine mass parameters in each subsample with a fit
    print 120 * '='
    print 'P2VV - INFO: createB2CCFitNTuple: fitting with simultaneous mass PDF'
    simMassFitResult = sWeightMassPdf.fitTo( dataSets['pre'][0], Save = True, **fitOpts )

    from P2VV.Imports import parValues
    simMassFitResult.PrintSpecial( text = True, LaTeX = True, normal = True, ParNames = parNames, ParValues = parValues )
    simMassFitResult.covarianceMatrix().Print()
    simMassFitResult.correlationMatrix().Print()

else :
    massPdfPars    = massPdf.getVariables()
    sWeightMassPdf = massPdf


###########################################################################################################################################
## compute sWeights ##
######################

print 'P2VV - INFO: createB2CCFitNTuple: computing sWeights'

# compute sWeights
from P2VV.Utilities.SWeights import SData
sWeightsBuilder = SData( Pdf = sWeightMassPdf, Data = dataSets['pre'][0], Name = 'JpsiKK' )
dataSets['preS'] = ( sWeightsBuilder.data(), [ ] )

# import sWeight variables
weightVars = [ ]
for wName in [ 'N_sigMass_sw', 'N_cbkgMass_sw' ] :
    ws.put( dataSets['preS'][0].get().find(wName) )
    weightVars.append( RealVar(wName) )


###########################################################################################################################################
## create signal and background data sets and print event yields ##
###################################################################

print 120 * '='
print 'P2VV - INFO: createB2CCFitNTuple: creating data set'

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
from ROOT import gROOT
gROOT.cd('PyROOT:/')

print 'P2VV - INFO: createB2CCFitNTuple: produced data set:\n' + ' ' * 13,
dataSets['main'][0].Print()
print

# create weighted data sets
from ROOT import RooArgSet
obsSets = dict( [ ( var.GetName(), RooArgSet( dataSets['main'][0].get() ) ) for var in weightVars ] )
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

print 'P2VV - INFO: createB2CCFitNTuple: signal data set:\n' + ' ' * 13,
dataSets['sigSWeight'][0].Print()
print
print 'P2VV - INFO: createB2CCFitNTuple: combinatorial background data set:\n' + ' ' * 13,
dataSets['cbkgSWeight'][0].Print()
print


# print yields
print 'P2VV - INFO: createB2CCFitNTuple: event yields:'
allCats = [ dataSets['main'][0].get().find( obsDict['runPeriod'][0] ) ] if runPeriods else [ ]
allCats += [  dataSets['main'][0].get().find( obsDict['hlt1ExclB'][0] )
            , dataSets['main'][0].get().find( obsDict['hlt2B'][0] )
           ]
if len(KKMassBinBounds) > 2 : allCats.append( dataSets['main'][0].get().find( observables['KKMassCat'].GetName() ) )
allCats = [ cat for cat in allCats if cat ]

from P2VV.Utilities.DataHandling import printEventYields, printEventYieldsData
printEventYields(  ParameterSet        = massPdfPars
                 , YieldNames          = yieldNames
                 , SplittingCategories = [ cat for catList in splitCats for cat in catList ]
                )
printEventYieldsData(  FullDataSet         = dataSets['main'][0]
                     , WeightedDataSets    = [ dataSets[name][0] for name in [ 'sigSWeight', 'cbkgSWeight' ] ]
                     , DataSetNames        = [ 'Signal', 'Combinatorial background' ]
                     , SplittingCategories = allCats
                    )


###########################################################################################################################################
## store n-tuple in ROOT file ##
################################

print 120 * '='

# create n-tuple containing signal and background weights
sigTree = dataSets['sigSWeight'][0].buildTree( Name = nTupleName, Title = nTupleName, WeightName = sWeightName, RooFitFormat = False )
print 'P2VV - INFO: createB2CCFitNTuple: n-tuple:'
print 'TTree::%s[%s] = %d entries'\
      % (  sigTree.GetName(), ','.join( br.GetName() for br in sigTree.GetListOfBranches() ), sigTree.GetEntries() )

# save n-tuple to file
print 'P2VV - INFO: createB2CCFitNTuple: saving n-tuple to ROOT file %s' % dataSetsFilePath
from ROOT import TFile
dataSetsFile = TFile.Open( dataSetsFilePath, 'UPDATE' if appendToFile else 'RECREATE' )
dataSetsFile.Append(sigTree)
dataSetsFile.Write()
dataSetsFile.Close()
