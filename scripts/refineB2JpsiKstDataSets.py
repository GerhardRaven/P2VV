import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-w', '--sigWeight',  type=str )
options = parser.parse_args()

# configure job
wghtArg            = options.sigWeight
if   's' in wghtArg or 'S' in wghtArg: weightName = 'Bs'
elif 'd' in wghtArg or 'D' in wghtArg: weightName = 'Bd'
else: assert False, 'P2VV - ERROR: Wrong sWeight speicifier. Choose either Bs or Bd.'

# input dataset paths
inPath       = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/'
dataSetPaths = {}
dataSetName  = 'Bs2JpsiKst'
dataSetPaths['2011negKaons'] = inPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_negKaons_fitNtuple_120614_%s_weighted.root'%weightName
dataSetPaths['2011posKaons'] = inPath + 'P2VVDataSet_2011Reco14_Bs2JpsiKst_posKaons_fitNtuple_120614_%s_weighted.root'%weightName
dataSetPaths['2012negKaons'] = inPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_negKaons_fitNtuple_120614_%s_weighted.root'%weightName
dataSetPaths['2012posKaons'] = inPath + 'P2VVDataSet_2012Reco14_Bs2JpsiKst_posKaons_fitNtuple_120614_%s_weighted.root'%weightName


outPath = '/project/bfys/vsyropou/data/Bs2JpsiKst/RealData/'
outFileName = 'P2VVDataSet_20112012Reco14_Bs2JpsiKst_allKaons_fitNtuple_120614_%s_weighted.root'%weightName
outDataName = 'Bs2JpsiKst'

# read datasets into a combined dataset
from ROOT import TFile, RooDataSet, RooFit
dataSets = {}
for key, path in dataSetPaths.iteritems(): dataSets[key] = TFile.Open(path).Get(dataSetName)

combinedDataSet = RooDataSet( outDataName, outDataName, dataSets['2012posKaons'].get() ) 
for dataSet in dataSets.itervalues(): combinedDataSet.append(dataSet)
combinedDataSet.Print()

outFile = TFile.Open(outPath + outFileName,'recreate')
combinedDataSet.Write()
outFile.Close()
del outFile

print 'P2VV - INFO: Wrote dataset to file: %s'%outPath + outFileName
