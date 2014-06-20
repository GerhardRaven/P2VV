import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-r', '--runPeriod',  type=int )
parser.add_argument( '-s', '--kaonSign',   type=int )
parser.add_argument( '-w', '--sigWeight',  type=str )
options = parser.parse_args()
for val in options.__dict__.itervalues():
    if val==None: assert False, 'P2VV - ERROR: Specify job parameters (runPeriod (-r) / kaonSign (-s) / sWeight (-w) ) with command arguments.'

# data set paths
inputPath      = '/data/bfys/vsyropou/Bs2JpsiKst/DiegosUnrefinedTuples/'
inputTreeName  = 'DecayTree'
protoTreePaths = {} 
protoTreePaths['2011negKaons'] = inputPath + '2011n.root'
protoTreePaths['2011posKaons'] = inputPath + '2011p.root'
protoTreePaths['2012negKaons'] = inputPath + '2012n.root'
protoTreePaths['2012posKaons'] = inputPath + '2012p.root'

# configure job
saveDataSetToLocalFolder = True   # protect official output location
delIntermidiatTree       = False  # delete temporary file

prodDate           = '120614'
runPeriods         = [2011,2012]
createFitNtuple    = True
createAngAccNtuple = True
addKpiMassCategory = True
addKaonSignInfo    = True
addRunPeriodInfo   = True
calculateHelAngles = True
saveOldHelAngles   = True 
applySelection     = False 
KpiMassBranchName  = 'Kst_892_0_M'
KpiMassBinBounds   = [826, 825+35, 892, 965-35, 966]
KpiMassInds        = range( len(KpiMassBinBounds) + 1 )
bdtgSelCuts        = {2011:0.2, 2012:-0.24}
runPeriod          = options.runPeriod
kaonSign           = options.kaonSign
kaonSignStr        = 'neg' if options.kaonSign < 0 else 'pos'
dataSetKey         = str(runPeriod) + kaonSignStr + 'Kaons'
protoTreePath      = protoTreePaths[dataSetKey]
wghtArg            = options.sigWeight
if   's' in wghtArg or 'S' in wghtArg: weightName, weightVarName = 'Bs', 'Bs_sWeight' 
elif 'd' in wghtArg or 'D' in wghtArg: weightName, weightVarName = 'Bd', 'Bd_sWeight'
else: assert False, 'P2VV - ERROR: Wrong sWeight speicifier. Choose either Bs or Bd.'

# output paths
outDataSetPaths = {}
outputPath      = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/' if not saveDataSetToLocalFolder else './'
outDataSetPaths['2011negKaons'] = outputPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_negKaons'
outDataSetPaths['2011posKaons'] = outputPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_posKaons'
outDataSetPaths['2012negKaons'] = outputPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_negKaons'
outDataSetPaths['2012posKaons'] = outputPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_posKaons'

# select and rename branches
from ROOT import RooNumber
from math import pi
RooInf  = RooNumber.infinity()
KKMin   = KpiMassBinBounds[0]
KKMax   = KpiMassBinBounds[-1]

obsDict = dict(   runPeriod         = ( 'runPeriod',       'run period', { 'p%s'%runPeriod : runPeriod for runPeriod in runPeriods }     )
                 , KpiMassCat       = ( 'KpiMassCat',      'KpiMassCat', {       'bin%s'%b : b  for b in range(1,len(KpiMassBinBounds))} )               
                 , kaonSign         = ( 'kaonSign',        'kaon sign',  {        'neg': -1, 'pos': 1}                                   )
                 , mass             = ( 'mass',            'm(J/#psi K^{+}#pi^{-})', 'MeV/c^{2}', 5368.,  5110.,   5690.        ) # 5110.,    5690.       )
                 , Mjpsik           = ( 'mJpsiK',          'm(J/#psi K^{+}',         'MeV/c^{2}', 4594,   -RooInf,  RooInf      )
                 , Kst_892_0_MM     = ( 'mdau2',           'm(K^{+}#pi^{-})',        'MeV/c^{2}', 892.,   KKMin,    KKMax       )
                 , J_psi_1S_MM      = ( 'mdau1',           'm(mu^{+}#mu^{-})',       'MeV/c^{2}', 3094,  -RooInf,    RooInf    )
                 #, sWeights_Bs      = ( 'sWeights_Bs',     'Bs sWeights',            '',           0.,    -RooInf, RooInf       )
                 #, sWeights_Bd      = ( 'sWeights_Bd',     'Bd sWeights',            '',           0.,    -RooInf, RooInf       )
                 , cor_sWeights_Bs  = ( 'Bs_sWeight', 'corrected Bs sWeights',       '',           0.,    -RooInf, RooInf       )
                 , cor_sWeights_Bd  = ( 'Bd_sWeight', 'corrected Bd sWeights',       '',           0.,    -RooInf, RooInf       )
                 , helcosthetaK     = ( 'helcosthetaK',    'cos(#theta_{K})',        '',           0.,      -1.,     +1.        )
                 , helcosthetaL     = ( 'helcosthetaL',    'cos(#theta_{#mu})',      '',           0.,      -1.,     +1.        )
                 , helphi           = ( 'helphi',          '#phi_{h}',               'rad',        0.,      -pi,     +pi        )
                 )

# add branches for angEff or helicity angles calculation
if createAngAccNtuple or calculateHelAngles:
    for name in [ '%s_P%s' % ( part, comp ) for part in [ 'muplus', 'muminus', 'Kplus', 'piminus' ] for comp in ( 'X', 'Y', 'Z' ) ]:
        obsDict[name] = ( name, name, 'MeV/c', 0.,   -RooInf, +RooInf )
#    if calculateHelAngles: 
#        for key in ['helcosthetaK','helcosthetaL','B0_Phi']: obsDict.pop(key)

# remove unecesessary weights 
if weightName == 'Bs': obsDict.pop('cor_sWeights_Bd')
if weightName == 'Bd': obsDict.pop('cor_sWeights_Bs')

# minimal dataset for fitting
if createFitNtuple: 
    for key in ['Mjpsik', 'J_psi_1S_MM']: obsDict.pop(key)

#####################################################################################################################
## refine all ntuple ##
#######################

# load p2vv library
from P2VV.Load import P2VVLibrary

# read input file
from ROOT import TFile
protoFile   = TFile.Open(protoTreePath)
protoTree   = protoFile.Get(inputTreeName)
initNumEvts = protoTree.GetEntries()
print 'P2VV - INFO: Processing tree (%s initial entries): %s'%(initNumEvts,protoTreePath)

# save old helicity angles in case you calculate new ones
if saveOldHelAngles:
    print 'P2VV - INFO: Saving old Helicity angles'
    protoTree.SetBranchStatus('*',0)
    for branch in ['helcosthetaK','helcosthetaL','B0_Phi']:protoTree.SetBranchStatus(branch,1)

    from ROOT import TFile
    oldHelAngFile = TFile.Open('oldHelAnglesFor_%s.root'%dataSetKey,'recreate')
    oldHelAngTree = protoTree.CloneTree()

    from ROOT import copyFloatInTree 
    for oldBranchName, newBranchName in {'helcosthetaK':'helcosthetaK_old', 'helcosthetaL':'helcosthetaL_old', 'B0_Phi':'helphi_old'}.iteritems():
        print 'P2VV - INFO: Copying and renaming branch: %s --> %s'%(oldBranchName,newBranchName)
        copyFloatInTree( oldHelAngTree, oldBranchName, newBranchName[0] )

    oldHelAngTree.Write()
    oldHelAngFile.Close()
    del oldHelAngFile

    protoTree.SetBranchStatus('*',1)
    for branch in ['helcosthetaK','helcosthetaL','B0_Phi']:protoTree.SetBranchStatus(branch,0)
    print 'P2VV - INFO: Wrote old helicity angles to file %s'%'oldHelAnglesFor_%s.root'%dataSetKey

# create temporary intermediate file
intermediateFileName = 'temp_FileFor_%s.root'%dataSetKey
intermediateFile = TFile.Open(intermediateFileName, 'recreate')
intermediateTree = protoTree.CloneTree()

protoFile.Close()
del protoFile

# create new branches
from ROOT import addIntegerToTree 
if addRunPeriodInfo:   
    print 'P2VV - INFO: Adding run period (%s) info to ntuple'%runPeriod
    addIntegerToTree(intermediateTree, int(runPeriod), 'runPeriod' )

if addKaonSignInfo:
    print 'P2VV - INFO: Adding kaon sign (%+i) info to ntuple'%kaonSign
    addIntegerToTree(intermediateTree, int(kaonSign), 'kaonSign' )

if calculateHelAngles: 
    from ROOT import TDatabasePDG, addHelicityAnglesToTree
    MeV = 1000 # TDatabasePDG is in GeV
    PDG = TDatabasePDG()
    Mmu = PDG.GetParticle('mu-').Mass()*MeV
    Mk  = PDG.GetParticle('K-').Mass()*MeV
    Mpi = PDG.GetParticle('pi-').Mass()*MeV
    
    print 'P2VV - INFO: Mass of daughter particles used for helicity angles calculation: m(mu) = %s, m(K) = %s, m(pi) = %s. Units MUST be in GeV!!'%(Mmu,Mk,Mpi)
    ## TODO: make this more flexible in h1<->h2 and remove old hel angles and add yours.
    addHelicityAnglesToTree(intermediateTree, 'Kplus', 'piminus', 'muplus', 'muminus',
                            Mk, Mpi, Mmu,Mmu, 
                            'helcosthetaK_test', 'helcosthetaL_test', 'helphi_test',
                            )
if addKpiMassCategory: 
    from ROOT import addCategoryToTree, std
    print 'P2VV - INFO: Adding Kpi-mass category with indices "%s" and KK-mass boundaries "%s" to n-tuple' % ( KpiMassInds, KpiMassBinBounds )
    bounds = std.vector('Double_t')()
    inds   = std.vector('Int_t')()
    for bound in KpiMassBinBounds : bounds.push_back(bound)
    for ind in KpiMassInds : inds.push_back(ind)
    addCategoryToTree( intermediateTree, KpiMassBranchName, 'KpiMassCat', bounds, inds )

# copy-rename branches
from ROOT import copyFloatInTree
for oldBranchName, newBranchName in obsDict.iteritems():
    if oldBranchName != newBranchName[0]:
        print 'P2VV - INFO: Copying and renaming branch: %s --> %s'%(oldBranchName,newBranchName[0])
        copyFloatInTree( intermediateTree, oldBranchName, newBranchName[0] )

# close intermediate files
intermediateTree.Write()
intermediateFile.Close()
del intermediateFile
print 'P2VV - INFO: Wrote intermediate tree to file: %s'%intermediateFileName
print 'P2VV - INFO: Finished refining tree. Continuing to create RooDataSet'        

##########################################################################################
## create RooDataSet ##
##########################################################################################

# create workspace
from P2VV.RooFitWrappers import RooObject, RealVar, Category
ws = RooObject(workspace = 'JpsiphiWorkspace').ws()

# create observables
observables  = { }
for obs in obsDict.keys():
    if type( obsDict[obs][2] ) == dict or type( obsDict[obs][2] ) == list :
        observables[obs] = Category( obsDict[obs][0], Title = obsDict[obs][1], Observable = True, States = obsDict[obs][2] )
    else :
        observables[obs] = RealVar( obsDict[obs][0], Title = obsDict[obs][1], Unit = obsDict[obs][2], Observable = True
                                   , Value = obsDict[obs][3], MinMax = ( obsDict[obs][4], obsDict[obs][5] ) )

# ntuple selection string
dataSetArgs = {}
if applySelection:
    from P2VV.Imports import cutSelStringsJpsiKst as selection
    selection.pop('noSelection')
    if runPeriod == 2011: selection.pop('bdtg_2012')
    if runPeriod == 2012: selection.pop('bdtg_2011')
    dataSetArgs['ntupleCuts'] =  ' && '.join( selection.itervalues() )

# read tree into RooDataSet
from P2VV.Utilities.DataHandling import readData 
outDataSet = readData( intermediateFileName, inputTreeName, 
                       NTuple       = True,
                       observables  = [ obs for obs in observables.itervalues()], 
                       WeightVar    = ( weightVarName,True ),            
                       ImportIntoWS = False,
                       **dataSetArgs
                       )
outDataSet.Print()
if initNumEvts != outDataSet.numEntries(): print 'P2VV -INFO: You lost %s entries due to selection cuts:'% (initNumEvts-outDataSet.numEntries()) 

# create final output file
outFileName = outDataSetPaths[dataSetKey]
if createFitNtuple: outFileName += '_fitNtuple'
if createAngAccNtuple or calculateHelAngles: outFileName += '_angEff'
outFileName += '_%s_%s.root'%(prodDate,weightName + '_weighted')

outFile = TFile.Open(outFileName, 'recreate')
outDataSet.Write()
outFile.Close()
del outFile

if delIntermidiatTree:
    import os
    os.remove(intermediateFileName)

