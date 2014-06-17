import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-r', '--runPeriod',  type=int )
parser.add_argument( '-s', '--kaonSign',   type=int )
parser.add_argument( '-w', '--sigWeight',  type=str )
options = parser.parse_args()
for val in options.__dict__.itervalues():
    if val==None: assert False, 'P2VV - ERROR: Specify job parameters (runPeriod (-r) / kaonSign (-s) / sWeight (-w) ) with command arguments.'

# data set paths
inputPath      = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/DiegosUnrefinedTuples/'
inputTreeName  = 'DecayTree'
protoTreePaths = {} 
protoTreePaths['2011negKaons'] = inputPath + '2011n.root'
protoTreePaths['2011posKaons'] = inputPath + '2011p.root'
protoTreePaths['2012negKaons'] = inputPath + '2012n.root'
protoTreePaths['2012posKaons'] = inputPath + '2012p.root'

outputPath      = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/'
outDataSetPaths = {}
outDataSetPaths['2011negKaons'] = outputPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_negKaons'
outDataSetPaths['2011posKaons'] = outputPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_posKaons'
outDataSetPaths['2012negKaons'] = outputPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_negKaons'
outDataSetPaths['2012posKaons'] = outputPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_posKaons'

# configure job
prodDate           = '120614'
runPeriods         = [2011,2012]
createFitNtuple    = True
createAngAccNtuple = True
addKpiMassCategory = True
addKaonSignInfo    = True
addRunPeriodInfo   = True
calculateHelAngles = False
applySelection     = False
KpiMassBranchName  = 'Kst_892_0_M'
KpiMassBinBounds   = [826, 825+35, 892, 965-35, 966]
KpiMassInds        = range( len(KpiMassBinBounds) + 1 )
bdtgSelCuts        = {2011:0.2, 2012:-0.24}
runPeriod          = options.runPeriod
kaonSign           = options.kaonSign
kaonSignStr        = 'neg' if options.kaonSign < 0 else 'pos'
dataSetKey         = str(runPeriod) + kaonSignStr + 'Kaons'
path               = protoTreePaths[dataSetKey]
wghtArg            = options.sigWeight
if   's' in wghtArg or 'S' in wghtArg: weightName, weightVarName = 'Bs', 'Bs_sWeight' 
elif 'd' in wghtArg or 'D' in wghtArg: weightName, weightVarName = 'Bd', 'Bd_sWeight'
else: assert False, 'P2VV - ERROR: Wrong sWeight speicifier. Choose either Bs or Bd.'

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
if createAngAccNtuple:
    for name in [ '%s_P%s' % ( part, comp ) for part in [ 'muplus', 'muminus', 'Kplus', 'piminus' ] for comp in ( 'X', 'Y', 'Z' ) ]:
        obsDict[name] = ( name, name, 'MeV/c', 0.,   -RooInf, +RooInf )

# remove unecesessary weights 
if weightName == 'Bs': obsDict.pop('cor_sWeights_Bd')
if weightName == 'Bd': obsDict.pop('cor_sWeights_Bs')

#####################################################################################################################
## refine all ntuple ##
#######################

# read input file
from ROOT import TFile
protoFile   = TFile.Open(path)
protoTree   = protoFile.Get(inputTreeName)

# create temporary intermediate file
intermediateFileName = 'temp_FileFor_%s.root'%dataSetKey
intermediateFile = TFile.Open(intermediateFileName, 'recreate')
intermediateTree = protoTree.CloneTree()

protoFile.Close()
del protoFile

# create new branches
from P2VV.Load import P2VVLibrary
from ROOT import addIntegerToTree 

print 'P2VV - INFO: Processing tree (%s initial entries): %s'%(intermediateTree.GetEntries(),path)
if addRunPeriodInfo:   
    print 'P2VV - INFO: Adding run period (%s) info to ntuple'%runPeriod
    addIntegerToTree(intermediateTree, int(runPeriod), 'runPeriod' )

if addKaonSignInfo:
    print 'P2VV - INFO: Adding kaon sign (%+i) info to ntuple'%kaonSign
    addIntegerToTree(intermediateTree, int(kaonSign), 'kaonSign' )

if calculateHelAngles: calculateHelicityAngles(intermediateTree)

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

# create final output file
outFileName = outDataSetPaths[dataSetKey]
if createFitNtuple: outFileName += '_fitNtuple'
if createAngAccNtuple: outFileName += '_angEff'
outFileName += '_%s_%s.root'%(prodDate,weightName + '_weighted')

outFile = TFile.Open(outFileName, 'recreate')
outDataSet.Write()
outFile.Close()
del outFile





# define helping functions
def addKpiMassCategory(tree, BinBounds, KpiMassBranchName):
    # create binning
    nBins = len(BinBounds) -1
    Bins = {}
    for i_thBin in range(nBins): Bins['bin%s'%i_thBin] = (BinBounds[i_thBin],BinBounds[i_thBin+1])
    
    def getBinNumber(KpiMass):
        binIndx = 0
        for bin, binRange in Bins.iteritems():
            if KpiMass > binRange[0] and KpiMass <= binRange[1]: 
                binIndx = bin[-1]
                break
        return int(binIndx)

    # create branch for Kpi mass category
    from array import array
    br_addr = array('i',[0])
    branch  = tree.Branch('KpiMassCat', br_addr, 'KpiMassCat/I')    
    
    # loop over treee and fill
    for event in tree:
        #KpiMass = getattr(tree, KpiMassBranchName)
        br_addr[0] = 1 #getBinNumber(KpiMass)
        branch.Fill()
