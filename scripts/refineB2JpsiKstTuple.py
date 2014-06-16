import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-r', '--runPeriod', type=int )
parser.add_argument( '-s', '--pionSign',  type=int )
options = parser.parse_args()

# load p2vv library
from ROOT import gSystem, TFile
gSystem.Load('libP2VV')

# data set paths
treeName = 'DecayTree'
protoTreePaths = {} 
protoTreePaths['2011negPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/temp/2011n.root'
protoTreePaths['2011posPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/temp/2011p.root'
protoTreePaths['2012negPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/temp/2012n.root'
protoTreePaths['2012posPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/temp/2012p.root'

outDataSetPaths = {}
outDataSetPaths['2011negPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/P2VVDataSet_2011Reco14_Bs2JpsiKst_negPions_'
outDataSetPaths['2011posPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/P2VVDataSet_2011Reco14_Bs2JpsiKst_posPions_'
outDataSetPaths['2012negPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/P2VVDataSet_2012Reco14_Bs2JpsiKst_negPions_'
outDataSetPaths['2012posPions'] = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/P2VVDataSet_2012Reco14_Bs2JpsiKst_posPions_'

# configure job
prodDate           = '_120614'
createFitNtuple    = True
createAngAccNtuple = True
addKpiMassCat      = True
addPionSignInfo    = True
addRunPeriodInfo   = True
calculateHelAngles = False
bdtgSelCuts        = {2011:0.2, 2012:-0.24}
KpiMassBinBounds   = [825, 825+40, 892, 965-40, 965]
KpiMassBranchName  = 'Kst_892_0_M'

# select and rename branches
obsDict = dict(
    Mjpsik          = 'mass',
    Kst_892_0_M     = 'mdau2',
    J_psi_1S_M      = 'mdau1',
    BDTG            = 'bdtg_sel',
    sWeights_Bd     = 'sWeights_Bd',
    cor_sWeights_Bd = 'sWeights_Bd_corr',
    sWeights_Bs     = 'sWeights_Bs',
    cor_sWeights_Bs = 'sWeights_Bs_corr',
    B0_Phi          = 'helphi',
    helcosthetaK    = 'helcosthetaK',
    helcosthetaL    = 'helcosthetaL'
    )

if createAngAccNtuple:
    for name in [ '%s_P%s' % ( part, comp ) for part in [ 'muplus', 'muminus' ] for comp in ( 'X', 'Y', 'Z' ) ]:
        obsDict[name] = name

# auto configure
runPeriod   = options.runPeriod
pionSign    = options.pionSign
pionSignStr = 'neg' if options.pionSign < 0 else 'pos'
dataSetKey = str(runPeriod) + pionSignStr + 'Pions'
path = protoTreePaths[dataSetKey]

# define helping functions
def addRunPeriodInfo(tree, period):
    from array import array
    br_addr = array('i',[int(period)])
    branch  = outTree.Branch('runPeriod', br_addr, 'runPeriod/I')
    for event in outTree: branch.Fill()

def addPionSignInfo(tree, sign):
    from array import array
    br_addr = array('i',[int(sign)])
    branch  = outTree.Branch('pionSign', br_addr, 'pionSign/I')
    for event in outTree: branch.Fill()

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

def calculateHelicityAngles(tree, h1name, h2name, l1name, l2name):
    a = 1
    ## implement this late, once you solve all other sbtlties


# loop to refine all trees
# for key, path in protoTreePaths.iteritems():

protoFile = TFile.Open(path,'read')
protoTree = protoFile.Get(treeName)
print 'P2VV - INFO: Processing tree (%s): %s'%(protoTree.GetEntries(),path)

if createAngAccNtuple:
    for name in [ '%s_P%s' % ( part, comp ) for part in [ 'Kplus', 'piminus' ] for comp in ( 'X', 'Y', 'Z' ) ]:
        obsDict[name] = name
        
# actually select and rename branches
protoTree.SetBranchStatus('*',0)
for branch in protoTree.GetListOfBranches(): 
    oldBrName = branch.GetName()
    if oldBrName in obsDict.keys():
        protoTree.SetBranchStatus(oldBrName,1)
        protoTree.GetBranch(oldBrName).SetName(obsDict[oldBrName])

# create output file        
outFileName = outDataSetPaths[dataSetKey]
if createFitNtuple: outFileName += 'fitNtuple' + prodDate
if createAngAccNtuple: outFileName += 'anfEff' + prodDate
outFileName += '.root'
outFile = TFile.Open(outFileName, 'recreate')

# apply selection
selection = '%s > %s'%(obsDict['BDTG'], bdtgSelCuts[runPeriod])
outTree = protoTree.CopyTree(selection)
print 'P2VV - INFO: Applying selection to initial tree, %s entries after selection.'%(outTree.GetEntries())

# create new branches 
if addRunPeriodInfo: addRunPeriodInfo(outTree, runPeriod)
if addPionSignInfo:  addPionSignInfo(outTree, pionSign)
if calculateHelAngles: calculateHelicityAngles(tree)      
if addKpiMassCategory:  addKpiMassCategory(outTree, KpiMassBinBounds, KpiMassBranchName )


# outTree.Show()


# from ROOT import TCanvas
# c = TCanvas('vd','csd')
# outTree.Draw('KpiMassCat')
# c.Print(str(pionSign) + runPeriod + '.pdf')

outFile.cd()
outTree.Write()
outFile.Close()
protoFile.Close()




# check the name conventions again
# try set entries
