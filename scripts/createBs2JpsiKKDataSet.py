###########################################################################################################################################
## set script parameters ##
###########################

# parse command-line options
import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '--jobName', '-N', default = 'Bs2JpsiKK' )
parser.add_argument( '--model', '-m', default = 'ipatia2' )  # 'ipatia2', 'doublegauss'
parser.add_argument( '--weightIn', '-x', default = 'pbkgWeight' )
parser.add_argument( '--biasedOnly', '-b', default = True )
parser.add_argument( '--minimalData', '-r', default = True )
parser.add_argument( '--addMomenta', '-a', default = False )
parser.add_argument( '--printYields', '-y', default = True )
parser.add_argument( '--numCPU', '-c', type = int, default = 2 )
parser.add_argument( '--runHesse', '-e', default = False )
parser.add_argument( '--runMinos', '-s', default = '' )
parser.add_argument( '--minosPars', '-p', default = None )
parser.add_argument( '--dataPath', '-d', default = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14' )
parser.add_argument( '--workPath', '-w', default = '/project/bfys/jleerdam/softDevel/P2VV2/test' )
parser.add_argument( '--nTupleFileIn', '-f', default = 'nTupleC_hope_the_last_one_add_20140309.root' )
parser.add_argument( '--nTupleNameIn', '-n', default = 'DecayTree' )
parser.add_argument( '--dataFormatOut', '-j', default = 'dataset' )  # 'dataset', 'ntuple'
parser.add_argument( '--dataSetFileOut', '-g' )
parser.add_argument( '--dataSetNameOut', '-q' )
parser.add_argument( '--parFileIn', '-i' )
parser.add_argument( '--parFileOut', '-o' )

args = parser.parse_args()
model = args.model.lower()
assert model in [ 'ipatia2', 'doublegauss' ]
biasedOnly = False if not args.biasedOnly or str( args.biasedOnly ).lower() in [ 'false', '0' ] else True
minimalData = False if not args.minimalData or str( args.minimalData ).lower() in [ 'false', '0' ] else True
addMomenta = False if ( not args.addMomenta and minimalData ) or str( args.addMomenta ).lower() in [ 'false', '0' ] else True
printYields = False if not args.printYields or str( args.printYields ).lower() in [ 'false', '0' ] else True
assert type(args.numCPU) == int and args.numCPU > 0 and args.numCPU < 20
runHesse = False if not args.runHesse or str( args.runHesse ).lower() in [ 'false', '0' ] else True
minosPars = args.minosPars.split(',') if args.minosPars and str(args.minosPars) != 'None' else [ ]
runMinos = False if ( args.runMinos == '' and not minosPars ) or str( args.runMinos ).lower() in [ 'false', '0' ] else True
dataPath = args.dataPath
workPath = args.workPath
if dataPath and dataPath[-1] != '/' : dataPath += '/'
if workPath and workPath[-1] != '/' : workPath += '/'
assert args.nTupleFileIn
nTupleFileIn = dataPath + args.nTupleFileIn
assert args.nTupleNameIn
dataFormatOut = args.dataFormatOut.lower()
dataSetFileOut = args.dataSetFileOut
if dataSetFileOut :
    dataSetFileOut = workPath + dataSetFileOut
else :
    dataSetFileOut = workPath + '%sdata.root' % ( args.jobName + '_' if args.jobName else '' )
dataSetNameOut = args.dataSetNameOut
if not dataSetNameOut :
    dataSetNameOut = 'JpsiKK_sigSWeight' if dataFormatOut == 'dataset' else 'DecayTree'
assert dataFormatOut in [ 'dataset', 'ntuple' ]
parFileIn = args.parFileIn
if parFileIn == None :
    parFileIn = workPath + 'eventYields.par'
elif parFileIn :
    parFileIn = workPath + parFileIn
parFileOut = args.parFileOut
if parFileOut :
    parFileOut = workPath + parFileOut

# print script settings
print 'job parameters:'
if args.jobName :
    print '  job name: %s' % args.jobName
print '  model: %s' % model
print '  input event weight: %s' % args.weightIn
print '  only HLT2-biased data: %s' % ( 'true' if biasedOnly else 'false' )
print '  minimal output data set: %s' % ( 'true' if minimalData else 'false' )
print '  add track momenta: %s' % addMomenta
print '  print event yields: %s' % printYields
print '  number of cores: %d' % args.numCPU
print '  run Hesse: %s' % ( 'true' if runHesse else 'false' )
print '  run Minos: %s' % ( 'true' if runMinos else 'false' )
print '  Minos parameters:', ', '.join( '"%s"' % par for par in minosPars )
print '  data path: %s' % dataPath
print '  work path: %s' % workPath
print '  n-tuple file in: %s' % nTupleFileIn
print '  n-tuple name in: %s' % args.nTupleNameIn
print '  data format out: %s' % dataFormatOut
print '  data set file out: %s' % dataSetFileOut
print '  data set name out: %s' % dataSetNameOut
print '  input parameter file: %s' % parFileIn
print '  output parameter file: %s' % parFileOut

# clear command-line options
import sys
sys.argv = sys.argv[ 0 : 1 ]

# additional settings
obsKeys = [  'runPeriod', 'hlt1ExclB', 'hlt2B', 'hlt2UB'
           , 'mass', 'KKMass', 'KKMassCat', 'mumuMass'
           , 'time', 'timeRes', 'ctk', 'ctl', 'phih'
           , 'wTagOS', 'wTagSS'
          ]
if not minimalData :
    obsKeys += [ 'wMC', 'sel', 'selClTail'
                , 'firstData', 'hlt2_prescale', 'polarity'
                , 'trigDecUnb', 'trigDecExclB'
                , 'tagDecOS', 'tagDecSS'
               ]
if addMomenta :
    obsKeys += [  'B_P', 'B_Pt', 'B_eta', 'B_phi'
                , 'Kplus_PX',   'Kplus_PY',   'Kplus_PZ',   'Kplus_LOKI_ETA'
                , 'Kminus_PX',  'Kminus_PY',  'Kminus_PZ',  'Kminus_LOKI_ETA'
                , 'muplus_PX',  'muplus_PY',  'muplus_PZ',  'muplus_LOKI_ETA'
                , 'muminus_PX', 'muminus_PY', 'muminus_PZ', 'muminus_LOKI_ETA'
               ]

from math import pi
from ROOT import RooNumber
RooInf  = RooNumber.infinity()
obsDict = dict(  runPeriod        = ( 'runPeriod',                 'run period',      dict( p2011 = 2011, p2012 = 2012      )       )
               , firstData        = ( 'firstData',                 'first data',      dict( first = +1, later = 0           )       )
               , hlt2_prescale    = ( 'hlt2_prescale',             'HLT2 prescale',   dict( presc = +1, noPresc = 0         )       )
               , polarity         = ( 'polarity',                  'magnet polarity', dict( up = +1, down = -1, noPol = 0   )       )
               , tagDecOS         = ( 'tagdecision_os_cb',         'OS tag decision', dict( B = +1, Bbar = -1, Untagged = 0 )       )
               , tagDecSS         = ( 'tagdecision_ss_nn',         'SS tag decision', dict( B = +1, Bbar = -1, Untagged = 0 )       )
               , sel              = ( 'sel',                       'selection',       dict( sel   = 1, notSel   = 0         )       )
               , selClTail        = ( 'sel_cleantail',             'clean tail sel.', dict( sel   = 1, notSel   = 0         )       )
               , hlt1ExclB        = ( 'hlt1_excl_biased',          'HLT1 excl. B.',   dict( exclB = 1, notExclB = 0         )       )
               , hlt2B            = ( 'hlt2_biased',               'HLT2 B.',         dict( B     = 1, notB     = 0         )       )
               , hlt2UB           = ( 'hlt2_unbiased',             'HLT2 UB.',        dict( UB    = 1, notUB    = 0         )       )
               , trigDecUnb       = ( 'triggerDecisionUnbiased',   'trigDecUnb',      dict( UB    = 1, notUB    = 0         )       )
               , trigDecExclB     = ( 'triggerDecisionBiasedExcl', 'trigDecExclB',    dict( exclB = 1, notExclB = 0         )       )
               , wMC              = ( 'wMC',                       'pbkgWeight',              '',          0.,    -RooInf, +RooInf  )
               , mass             = ( 'mass',                      'm(J/#psi K^{+}K^{-})',    'MeV/c^{2}', 5368.,  5200.,   5550.   )
               , mumuMass         = ( 'mdau1',                     'm(#mu^{+}#mu^{-})',       'MeV/c^{2}', 3096.,  3030.,   3150.   )
               , KKMass           = ( 'mdau2',                     'm(K^{+}K^{-})',           'MeV/c^{2}', 1020.,  990.,    1050.   )
               , time             = ( 'time',                      'Decay time',              'ps',        1.5,    0.3,     14.     )
               , timeRes          = ( 'sigmat',                    '#sigma(t)',               'ps',        0.01,   0.0001,  0.12    )
               , ctk              = ( 'helcosthetaK',              'cos(#theta_{K})',         '',          0.,    -1.,     +1.      )
               , ctl              = ( 'helcosthetaL',              'cos(#theta_{#mu})',       '',          0.,    -1.,     +1.      )
               , phih             = ( 'helphi',                    '#phi_{h}',                'rad',       0.,    -pi,     +pi      )
               , wTagOS           = ( 'tagomega_os_cb',            'OS est. wrong-tag prob.', '',          0.25,   0.,      0.50001 )
               , wTagSS           = ( 'tagomega_ss_nn',            'SS est. wrong-tag prob.', '',          0.25,   0.,      0.50001 )
               , muplus_PX        = ( 'muplus_PX',                 'muplus_PX',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , muplus_PY        = ( 'muplus_PY',                 'muplus_PY',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , muplus_PZ        = ( 'muplus_PZ',                 'muplus_PZ',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , muplus_LOKI_ETA  = ( 'muplus_LOKI_ETA',           'muplus_LOKI_ETA',         '',          0.,    -RooInf, +RooInf  )
               , muminus_PX       = ( 'muminus_PX',                'muminus_PX',              'MeV/c',     0.,    -RooInf, +RooInf  )
               , muminus_PY       = ( 'muminus_PY',                'muminus_PY',              'MeV/c',     0.,    -RooInf, +RooInf  )
               , muminus_PZ       = ( 'muminus_PZ',                'muminus_PZ',              'MeV/c',     0.,    -RooInf, +RooInf  )
               , muminus_LOKI_ETA = ( 'muminus_LOKI_ETA',          'muminus_LOKI_ETA',        '',          0.,    -RooInf, +RooInf  )
               , Kplus_PX         = ( 'Kplus_PX',                  'Kplus_PX',                'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kplus_PY         = ( 'Kplus_PY',                  'Kplus_PY',                'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kplus_PZ         = ( 'Kplus_PZ',                  'Kplus_PZ',                'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kplus_LOKI_ETA   = ( 'Kplus_LOKI_ETA',            'Kplus_LOKI_ETA',          '',          0.,    -RooInf, +RooInf  )
               , Kminus_PX        = ( 'Kminus_PX',                 'Kminus_PX',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kminus_PY        = ( 'Kminus_PY',                 'Kminus_PY',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kminus_PZ        = ( 'Kminus_PZ',                 'Kminus_PZ',               'MeV/c',     0.,    -RooInf, +RooInf  )
               , Kminus_LOKI_ETA  = ( 'Kminus_LOKI_ETA',           'Kminus_LOKI_ETA',         '',          0.,    -RooInf, +RooInf  )
               , B_Pt             = ( 'B_Pt',                      'B_Pt',                    'MeV/c',     0.,     0.,      RooInf  )
               , B_P              = ( 'B_P',                       'B_P',                     'MeV/c',     0.,     0.,      RooInf  )
               , B_eta            = ( 'B_eta',                     'B_eta',                   '',          0.,    -RooInf, +RooInf  )
               , B_phi            = ( 'B_phi',                     'B_phi',                   'rad',       0.,    -pi,     +pi      )
              )

constSplitPars = [ 'm_sig_alpha_1', 'm_sig_alpha_2', 'm_sig_n_1', 'm_sig_n_2', 'm_sig_lambda', 'm_sig_zeta' ]
constSplitVals = {  'bin0' : [ 3.1265598805552086, 2.4478555551857504, 0.01,                0.4626575371693939,  -1.940905888712173, 0.0010225614571429675  ]
                  , 'bin1' : [ 2.235667390313107,  2.296350006990024,  1.1725685774531813,  1.696945683056572,   -2.590902278024503, 0.2967019395179116     ]
                  , 'bin2' : [ 2.8232366360188683, 2.570677068540935,  2.046626071513164,   2.6043347350995187,  -2.76448170618875,  0.00019473384424262274 ]
                  , 'bin3' : [ 3.1064191306717484, 2.096424131956603,  1.69939291075719,    2.7310345435871,     -2.903782724519466, 0.0007918143518637066  ]
                  , 'bin4' : [ 2.5929196863701067, 2.1195375110354044, 1.2511053016584157,  1.432994401017839,   -2.626462929014679, 0.19046609538221904    ]
                  , 'bin5' : [ 3.0236236294665466, 3.033593188175826,  0.39183360830444464, 0.43646654437298116, -2.216215361029956, 0.00019353774759423015 ]
                 }


###########################################################################################################################################
## read data ##
###############

from P2VV.Load import RooFitOutput, LHCbStyle

# create list of required observables
reqObsList  = [ 'index', 'runPeriod', 'mass', 'KKMass', 'tagDecOS', 'tagDecSS', 'wTagOS', 'wTagSS', 'hlt1ExclB', 'hlt2B' ]
reqObsList += [ args.weightIn ] if args.weightIn else [ ]
reqObsList += [ '%s_P%s' % ( part, comp ) for part in [ 'Kplus', 'Kminus', 'muplus', 'muminus' ] for comp in ( 'X', 'Y', 'Z' ) ]\
              if addMomenta else [ ]

# create workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject(workspace = 'JpsiKKWS').ws()

# create observables
from P2VV.RooFitWrappers import RealVar, Category
observables  = { }
observables['index'] = RealVar( 'index', Title = 'event index', Observable = True, Value = 0., MinMax = ( -RooInf, +RooInf ) )
obsSetPreDS  = [ observables['index'] ]
if args.weightIn :
    observables[args.weightIn] = RealVar( args.weightIn, Title = 'weight', Observable = True, Value = 0., MinMax = ( -RooInf, +RooInf ) )
    obsSetPreDS.append( observables[args.weightIn] )
obsSetNTuple = [ ]
for obs in obsKeys + reqObsList :
    if obs in observables : continue
    if obs == 'KKMassCat' : continue

    if type( obsDict[obs][2] ) == dict or type( obsDict[obs][2] ) == list :
        observables[obs] = Category( obsDict[obs][0], Title = obsDict[obs][1], Observable = True, States = obsDict[obs][2] )

    else :
        observables[obs] = RealVar( obsDict[obs][0], Title = obsDict[obs][1], Unit = obsDict[obs][2], Observable = True
                                   , Value = obsDict[obs][3], MinMax = ( obsDict[obs][4], obsDict[obs][5] ) )

    if obs in obsKeys    : obsSetNTuple.append( observables[obs] )
    if obs in reqObsList : obsSetPreDS.append(  observables[obs] )

# add mass ranges
observables['mass'].setRanges( dict(  LeftSideBand  = ( 5200., 5300. )
                                    , Signal        = ( 5300., 5440. )
                                    , RightSideBand = ( 5440., 5550. )
                                    , PeakBkg       = ( 5380., 5440. )
                                   )
                             )

# read n-tuple file(s)
from P2VV.Imports import cutSelStrings, triggerSelStrings
ntupleCuts = '%s && %s' % ( cutSelStrings['nominal2011'], triggerSelStrings[ 'paper2012_tos' if biasedOnly else 'timeEffFit_tos' ] )
print 'P2VV - INFO: createBs2JpsiKKDataSet: reading n-tuple "%s" from file "%s"' % ( args.nTupleNameIn, nTupleFileIn )
print 'P2VV - INFO: createBs2JpsiKKDataSet: selection cuts:\n    %s' % ntupleCuts
from ROOT import TFile, TTree
dataTreeFile = TFile.Open(nTupleFileIn)
assert dataTreeFile, 'P2VV - ERROR: createBs2JpsiKKDataSet: could not open file "%s"' % nTupleFileIn
dataTree = dataTreeFile.Get(args.nTupleNameIn)
assert dataTree, 'P2VV - ERROR: createBs2JpsiKKDataSet: could not locate tree "%s" in file "%s"' % ( args.nTupleNameIn, nTupleFileIn )

# create data set from n-tuple
preDataSet = dataTree.buildDataSet( Observables = obsSetPreDS, Name = 'JpsiKK', Title = 'JpsiKK', Cuts = ntupleCuts
                                   , IndexName = 'index', WeightName = args.weightIn if args.weightIn else '' )
print 'P2VV - INFO: createBs2JpsiKKDataSet: data set from n-tuple:\n' + ' ' * 13,
preDataSet.Print()
print

# create KK mass binning
from array import array
KKMassBinsArray = array( 'd', [ 990., 1008., 1016, 1020., 1024., 1032., 1050. ] )
from ROOT import RooBinning
KKMassBinning = RooBinning( len(KKMassBinsArray) - 1, KKMassBinsArray, 'KKMassBinning' )
observables['KKMass'].setBinning( KKMassBinning, 'KKMassBinning' )

# add KK mass split category to data
from P2VV.RooFitWrappers import BinningCategory
observables['KKMassCat'] = BinningCategory( 'KKMassCat', Observable = observables['KKMass'], Binning = KKMassBinning, Fundamental = True
                                           , Data = [preDataSet], CatTypeName = 'bin' )
if 'KKMassCat' in obsKeys : obsSetNTuple.append( observables['KKMassCat'] )


###########################################################################################################################################
## build J/psiKK mass PDFs ##
#############################

# initialize PDF components
sWeightNames = [ 'N_sigMass_sw', 'N_cbkgMass_sw' ]
from P2VV.RooFitWrappers import Component
sigMassComps  = Component( 'sigMass',  [ ], Yield = ( 0.4 * preDataSet.sumEntries(), 0., preDataSet.sumEntries() ) )
cbkgMassComps = Component( 'cbkgMass', [ ], Yield = ( 0.6 * preDataSet.sumEntries(), 0., preDataSet.sumEntries() ) )
massComps = [ sigMassComps, cbkgMassComps ]

# build the signal mass PDF
sigMassArgs = dict( Name = 'sig_m', mass = observables['mass'] )
if model == 'ipatia2' :
    from P2VV.Parameterizations.MassPDFs import Ipatia2_Signal_Mass as SignalBMass
    sigMassArgs['m_sig_lambda']  = -2.5
    sigMassArgs['m_sig_zeta']    =  0.01
    sigMassArgs['m_sig_alpha_1'] =  3.0
    sigMassArgs['m_sig_alpha_2'] =  2.5
    sigMassArgs['m_sig_n_1']     =  1.5
    sigMassArgs['m_sig_n_2']     =  2.0
else :
    from P2VV.Parameterizations.MassPDFs import DoubleGauss_Signal_Mass as SignalBMass
signalBMass = SignalBMass( **sigMassArgs )
sigMassComps += signalBMass.pdf()

# build the combinatorial background mass PDF
from P2VV.Parameterizations.MassPDFs import LP2011_Background_Mass as BackgroundBMass
backgroundBMass = BackgroundBMass( Name = 'cbkg_m', mass = observables['mass'] )
cbkgMassComps += backgroundBMass.pdf()

# build total mass PDF
from P2VV.RooFitWrappers import buildPdf
massPdf = buildPdf( massComps, Observables = [ observables['mass'] ], Name = 'JpsiKKMass' )

# get splitting categories and parameters
yieldNames  = [ par.GetName() for par in massPdf.Parameters() if par.getAttribute('Yield') ]
allParNames = [ par.GetName() for par in massPdf.Parameters() if not par.getAttribute('Yield') if not par.isConstant() ]
sigParNames = [ par.GetName() for par in signalBMass.parameters() if not par.isConstant() ]
bkgParNames = [ par.GetName() for par in backgroundBMass.parameters() if not par.isConstant() ]
splitCats = {  observables['runPeriod'].GetName() : set( yieldNames + bkgParNames )
             , observables['KKMassCat'].GetName() : set( yieldNames               )
            }
if not biasedOnly :
    splitCats[ observables['hlt2B'].GetName() ] = set( yieldNames + bkgParNames )
if model == 'ipatia2' :
    # split Ipatia 2 signal parameters
    hlt2BName = observables['hlt2B'].GetName()
    KKCatName = observables['KKMassCat'].GetName()
    for parName in sigParNames :
        if not biasedOnly and 'sigma' in parName : splitCats[hlt2BName].add(parName)
        splitCats[KKCatName].add(parName)
    for parName in constSplitPars :
        splitCats[KKCatName].add(parName)

splitParsDict = { }
pdfVars = massPdf.getVariables()
for cat, params in splitCats.iteritems() :
    assert ws.cat(cat), 'P2VV - ERROR: createBs2JpsiKKDataSet: category "%s" not in workspace' % cat
    for par in params :
        assert ws.var(par), 'P2VV - ERROR: createBs2JpsiKKDataSet: no variable "%s" in workspace' % par
        assert pdfVars.find(par), 'P2VV - ERROR: createBs2JpsiKKDataSet: variable "%s" not in PDF' % par
        if ws[par] not in splitParsDict :
            splitParsDict[ ws[par] ] = set( [ ws[cat] ] )
        else :
            splitParsDict[ ws[par] ].add( ws[cat] )

# create list of split parameters and categories
from P2VV.Utilities.General import createSplitParsList
splitPars = createSplitParsList(splitParsDict)

# build simultaneous mass PDF
print 'P2VV - INFO: createBs2JpsiKKDataSet: building simultaneous PDF "%s":' % ( massPdf.GetName() + '_simul' )
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
massPdfObs  = sWeightMassPdf.getObservables(preDataSet)
massPdfPars = sWeightMassPdf.getParameters(preDataSet)

if model == 'ipatia2' :
    # set values of constant Ipatia 2 parameters
    from P2VV.Utilities.General import getSplitPar
    for splitCat in sWeightMassPdf.indexCat() :
        sWeightMassPdf.indexCat().setIndex( splitCat.getVal() )
        KKState = observables['KKMassCat'].getLabel()
        for it, parName in enumerate(constSplitPars) :
            par = getSplitPar( parName, KKState, massPdfPars )
            par.setVal( constSplitVals[KKState][it] )

if parFileIn or parFileOut :
    # create PDF configuration object
    from P2VV.Parameterizations.FullPDFs import PdfConfiguration as PdfConfig
    pdfConfig = PdfConfig()

print 120 * '='
if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(sWeightMassPdf)
print 'P2VV - INFO: createBs2JpsiKKDataSet: observables in mass PDF:'
massPdfObs.Print('v')
print 'P2VV - INFO: createBs2JpsiKKDataSet: parameters in mass PDF:'
massPdfPars.Print('v')

# determine mass parameters with a fit
print 'P2VV - INFO: createBs2JpsiKKDataSet: fitting mass PDF'
simMassFitResult = sWeightMassPdf.fitTo( preDataSet, Save = True, SumW2Error = False, NumCPU = args.numCPU, Optimize = 1, Timer = True
                                        , Minimizer = 'Minuit2', Strategy = 1, Offset = True, Hesse = runHesse, Minos = minosPars )
simMassFitResult.PrintSpecial( text = True )

if parFileOut :
    # write parameters to file
    pdfConfig.getParametersFromPdf( sWeightMassPdf, preDataSet )
    pdfConfig.writeParametersToFile( filePath = parFileOut )


###########################################################################################################################################
## compute sWeights ##
######################

print 120 * '='
print 'P2VV - INFO: createBs2JpsiKKDataSet: computing sWeights'

# compute sWeights
from P2VV.Utilities.SWeights import SData
sWeightsBuilder = SData( Pdf = sWeightMassPdf, Data = preDataSet, Name = 'JpsiKK' )
preSDataSet = sWeightsBuilder.data()

# import sWeight variables
weightVars = [ ]
for wName in sWeightNames :
    ws.put( preSDataSet.get().find(wName) )
    weightVars.append( RealVar(wName) )


###################################################################################################################################
## add tagging observables to data sets ##
##########################################

print 120 * '='
print 'P2VV - INFO: createBs2JpsiKKDataSet: building tagging categories'

# tagging observable names
wTagOSName   = obsDict['wTagOS'][0]
wTagSSName   = obsDict['wTagSS'][0]
tagDecOSName = obsDict['tagDecOS'][0]
tagDecSSName = obsDict['tagDecSS'][0]

# get tagging category bins
from P2VV.Parameterizations.FlavourTagging import getTagCatParamsFromData as getTagParams
tagCats = [ ( 'Untagged', 0, 0.500000001 ), ( 'Tagged',   1, 0.499999999 ) ]
tagBinsOS = getTagParams( preSDataSet, estWTagName = wTagOSName, tagCats = tagCats, numSigmas = 1., SameSide = False
                         , WeightVarName = weightVars[0].GetName() )
tagBinsSS = getTagParams( preSDataSet, estWTagName = wTagSSName, tagCats = tagCats, numSigmas = 1., SameSide = True
                         , WeightVarName = weightVars[0].GetName() )

# add tagging categories to data sets
from P2VV.Utilities.DataHandling import addTaggingObservables
addTaggingObservables( preSDataSet, 'iTagOS', 'tagCatP2VVOS', tagDecOSName, wTagOSName, tagBinsOS )
addTaggingObservables( preSDataSet, 'iTagSS', 'tagCatP2VVSS', tagDecSSName, wTagSSName, tagBinsSS )

observables['iTagOS']       = Category( ws.put( preSDataSet.get().find('iTagOS')       ).GetName() )
observables['iTagSS']       = Category( ws.put( preSDataSet.get().find('iTagSS')       ).GetName() )
observables['tagCatP2VVOS'] = Category( ws.put( preSDataSet.get().find('tagCatP2VVOS') ).GetName() )
observables['tagCatP2VVSS'] = Category( ws.put( preSDataSet.get().find('tagCatP2VVSS') ).GetName() )
obsSetNTuple += [ observables['iTagOS'], observables['iTagSS'], observables['tagCatP2VVOS'], observables['tagCatP2VVSS'] ]


###################################################################################################################################
## add track and B momenta to data set ##
#########################################

if addMomenta :
    # add K+, K-, mu+ and mu- momentum magnitudes to data set
    from P2VV.Utilities.DataHandling import addTrackMomenta
    addTrackMomenta(preSDataSet)

    observables['Kplus_P']   = RealVar( ws.put( preSDataSet.get().find('Kplus_P')   ).GetName() )
    observables['Kminus_P']  = RealVar( ws.put( preSDataSet.get().find('Kminus_P')  ).GetName() )
    observables['muplus_P']  = RealVar( ws.put( preSDataSet.get().find('muplus_P')  ).GetName() )
    observables['muminus_P'] = RealVar( ws.put( preSDataSet.get().find('muminus_P') ).GetName() )
    obsSetNTuple += [ observables['Kplus_P'], observables['Kminus_P'], observables['muplus_P'], observables['muminus_P'] ]


###########################################################################################################################################
## create signal and background data sets and print event yields ##
###################################################################

print 120 * '='
print 'P2VV - INFO: createBs2JpsiKKDataSet: creating data set'

# create data sets with final columns
obsSetMain = obsSetNTuple + [ var for var in weightVars ]
mainDataSet = dataTree.buildDataSet( Observables = obsSetMain, Name = 'JpsiKK', Title = 'JpsiKK', IndexName = 'index'
                                    , OrigDataSet = preSDataSet )
del preSDataSet
del preDataSet
dataTreeFile.Close()
from ROOT import gROOT
gROOT.cd('PyROOT:/')

from ROOT import RooArgSet, RooDataSet
dataSets = { }
for it, var in enumerate(weightVars) :
    argSet = RooArgSet( mainDataSet.get() )
    for remVar in weightVars :
        if remVar != var : argSet.remove( argSet.find( remVar.GetName() ) )
    dataName = dataSetNameOut if it == 0 else 'JpsiKK_%d' % it
    dataSets[var] = RooDataSet( dataName, dataName, argSet, Import = mainDataSet, WeightVar = ( var.GetName(), True ) )
    if not printYields : break

print 'P2VV - INFO: createBs2JpsiKKDataSet: produced data sets:'
mainDataSet.Print()
for var in weightVars :
    dataSets[var].Print()
    if not printYields : break
print

if printYields :
    # print event yields
    print 'P2VV - INFO: createBs2JpsiKKDataSet: event yields:'
    allCats = [  mainDataSet.get().find( obsDict['runPeriod'][0] )
               , mainDataSet.get().find( obsDict['hlt1ExclB'][0] )
               , mainDataSet.get().find( obsDict['hlt2B'][0] )
               , mainDataSet.get().find( observables['KKMassCat'].GetName() )
              ]
    from P2VV.Utilities.DataHandling import printEventYields, printEventYieldsData
    printEventYields(  ParameterSet        = massPdfPars
                     , YieldNames          = yieldNames
                     , SplittingCategories = [ cat for pars in splitPars for cat in pars[1] ]
                    )
    printEventYieldsData(  FullDataSet         = mainDataSet
                         , WeightedDataSets    = [ dataSets[var] for var in weightVars ]
                         , DataSetNames        = [ 'Signal', 'Combinatorial background' ]
                         , SplittingCategories = allCats
                        )

print 120 * '='
print 'P2VV - INFO: createBs2JpsiKKDataSet: saving data sets to ROOT file %s' % dataSetFileOut
outFile = TFile.Open( dataSetFileOut, 'RECREATE' )
if dataFormatOut == 'dataset' :
    outFile.Add( dataSets[ weightVars[0] ] )
else :
    wName = 'sWeights_ipatia' if model == 'ipatia2' else 'sWeights_dgauss'
    nTuple = dataSets[ weightVars[0] ].buildTree( Name = dataSetNameOut, Title = dataSetNameOut, WeightName = wName, RooFitFormat = False )
    print 'TTree::%s[%s] = %d entries'\
          % (  nTuple.GetName(), ','.join( br.GetName() for br in nTuple.GetListOfBranches() ), nTuple.GetEntries() )
from ROOT import TObject
outFile.Write( dataSetFileOut, TObject.kOverwrite )
outFile.Close()
