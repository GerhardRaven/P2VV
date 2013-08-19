

# What this script does:
# It takes as imput a MC sample and an sWeighted dataset.
# And does the folowing:
# (1) Constructs a MC pdf                                                                         ( DONE )
# (2) Reweights the MC sample with a weight w = pdf(t,Omega;pars_sData) / pdf(t,Omega;pars_MC)    ( DONE )
#  in order to match the physics between data and MC                                              ( DONE )
# (3) Writes the output reweighted MC sample in 'physReweightOutputFile'                          ( DONE )
# (4) Takes the physics reweighted MC sample and reweights it again to mach the K momentum        ( DONE )
# (5) Takes the physics + K momentum matched MC sample and calculates the angular acceptance      ( Under Devel ) 
# (6) Performs the sFit on data using the new accptance.                                          ( Under Devel )
# At this stage a loop takes you to step (1)                                                      ( Under Devel )       


# NOTE: VERY IMPORTANT: The KK  matching reweight cannot handle weighted datasets at this stage. 
#  in orderr to prceed the weighted datasets are mimiced (see method matchWeightedDistributions.MimicWeightedDistribution ).
#  As a result the output datsets have less entries. It has not been checked if this increses the statistical error of the
#  acceptance function. 
#----------------------------------------------------------------------------------------------------------------------------





#Import stuff
# MUST SetUp Urania FIRST and then P2VV
import sys
import os    # exit
sys.path.append(os.environ["BS2MUMUROOT"] +"/python/Bs2MuMu/")

# Global output folder and default tree name
# Modify this according to where you want to put your output datasets.
globalOutputFolder = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/output/'
# Default tree name
treeName = 'DecayTree'

# Specify input data
# Mc input File 
mcTuplePath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628.root'
# sData input file
sDatPath = '/project/bfys/vsyropou/PhD/macros/iterativeAngAcc/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_sWeights.root'









####################################### Begin iterative procedure ###########################################################
# Output file name 
physReweightOutputFile = 'Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628_physWeights_1stIteration.root'

#Parameters obtained from sFit to data, on KK mass binning, wide KK mass window. 
dataParameters = dict( prefix       = 'data_'
                       ,AperpMag2   =  0.246729
                       ,AperpPhase  =  3.08416
                       ,A0Mag2      =  0.523498
                       ,A0Phase     =  0 # cosntrained
                       ,AparPhase   =  3.2097
                       ,ASOddPhase  =-0.0988961
                       ,f_S         = 0.0471
                       ,C_SP        = 0.326
                       ,dM          = 17.688
                       ,dGamma      = 0.102546
                       ,Gamma       = 0.672938
                       ,phiCP       = 0.0899423
                       ,lambdaCPSq  = 0.901019**2
                           )
#### Build angular pdf, calculate and write physics matcing weights
from P2VV.GeneralUtils import matchMCphysics2Data
matchPhysics = matchMCphysics2Data( mcTuplePath )
matchPhysics.buildMonteCarloPdf()
matchPhysics.calculateWeights( dataParameters )
matchPhysics.writeWeightsToFile( globalOutputFolder + physReweightOutputFile )

# Comment in to make compare B momenta to see if physics match between data and MC.
from P2VV.GeneralUtils import compareWeightedDistributions
from ROOT import TFile
#t_mc = TFile.Open(globalOutputFolder + physReweightOutputFile).Get(treeName)
#t_sD = TFile.Open(sDatPath).Get(treeName)
#c,a = compareWeightedDistributions(t_mc, t_sD, 'B_P', weight='weightPhys', sWeight='sWeight', assymPlot=True )



#### Reweight Kaon momenta of MC to match the Kaon momenta of sData.
momMatcingInputFile   = globalOutputFolder + physReweightOutputFile 
momReweightOutputFile = globalOutputFolder +  'Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20130628_physWeights_momWeights_1stIteration'
from P2VV.GeneralUtils import matchWeightedDistributions
reweightArgs = dict( sDInfo     = dict(path=sDatPath,            name=treeName, weight='sWeight'    )
                     ,mcInfo    = dict(path=momMatcingInputFile, name=treeName, weight='weightPhys' )
                     ,whichVars = 'KaonMomenta'
                     ,nBins     = 500  ### WATCH OUT IN THE END YOU WILL SET THIS TO 1000  
                     )
matchMC2Data = matchWeightedDistributions( momReweightOutputFile, **reweightArgs )
matchMC2Data.mimicWeights()
matchMC2Data.reweightMC()

# Check if K_minus before and after reweighting DO match.
t_mc = TFile.Open(momReweightOutputFile + '.root').Get('T')
t_sD = TFile.Open(sDatPath).Get(treeName)
c,a = compareWeightedDistributions(t_mc, t_sD, 'kminus_pmod', sVar='Kminus_P', weight='weightPhys', sWeight='sWeight', rangeX=[0,10e4], assymPlot=True )








































########## Stuff for xChecking #######################################################################

# Check physics reweighting
dataInfdict = dict(arg=True, path=sDatPath, name='DecayTree',sweight='sWeight' )
c, a= matchPhysics.comparisionPlots(kinematicVars['Kpl_P'], compWithData=dataInfdict, Xrange=[0,8e3] )
c, a= matchPhysics.comparisionPlots(kinematicVars['Kpl_P'], Xrange=[0,8e3] )



# Check mimicing of the weighted dataset
from ROOT import TFile
t = TFile.Open('output/Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010_weights.root').Get('DecayTree')
t.Draw('Kplus_PX>>hm(100,-14421,12580)','weightPhys')
from ROOT import gPad
hm = gPad.GetPrimitive('hm')
from Utilities import getSumOfWeights
print hm.GetSumOfWeights()

from ROOT import TH1F
test = TH1F('test','test',100,-14421,12580)
for i in xrange( len(matchMC2Data._mimicedVars['mc']['Kplus_PX']) ):
    test.Fill( matchMC2Data._mimicedVars['mc']['Kplus_PX'][i] )
print len(matchMC2Data._mimicedVars['mc']['Kplus_PX'])






# Check KK momomnta reweighitg
from ROOT import TFile
tmc = TFile.Open('Bs2JpsiPhiPrescaled_MC11a_ntupleB_for_fitting_20121010_Reweighted_1.root').Get('T')
tmc.Draw('kplus_pmod','weightPhys')

tsD = TFile.Open('Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp_sWeights.root').Get('DecayTree')
K_mom = dict(Kpl_P = 'TMath::Sqrt(Kplus_PX**2 + Kplus_PY**2 + Kplus_PZ**2)' )
tsD.Draw(K_mom['Kpl_P'],'sWeight')

from Utilities import getSumOfWeights
tmc.getSumOfWeights(tmc,'weightPhys')
tsD.getSumOfWeights(tsD,'sWeight')
#scale = 











## A method to make comparision plots
   # def comparisonPlots( self, whichVar, **kwargs):
    #     Xrange    = kwargs.pop( 'rangeX',        []                                         )
    #     sDataInfo = kwargs.pop( 'compWith',  dict(arg=False,path=None,name=None,sweight='') )
    #     AssymPlot = kwargs.pop( 'assymPlot', False                                          )
    #     # arg flag==True: Compares sData with Mc after physics reweighting
    #     # arg flag==False: Compares MC before and after physics reweighting

    #     if sDataInfo['arg']: 
    #         arg,p,n,s=sDataInfo['arg'],sDataInfo['path'],sDataInfo['name'],sDataInfo['sweight']
    #     else: arg=False

    #     from ROOT import TFile
    #     if arg: 
    #         file_sdata = TFile.Open(p)        
    #         tree_sdata = file_sdata.Get(n)        
        
    #     file_mc = TFile.Open(self._weightedNtuplePath)
    #     tree_mc = file_mc.Get(self._weightedNtupleName)
        
    #     if arg:
    #         print 'arg=True'
    #         print 'P2VV - INFO: Comparing ' + whichVar + ' between MC and sWeighted data after MC physics reweighting.'
    #         if AssymPlot:
    #               c,d = compareWeightedDistributions(tree_mc, tree_sdata, whichVar,                   \
    #                     weight=self._physWeightsName, sWeight=s, rangeX=Xrange, assymPlot=AssymPlot )
    #         else: c = compareWeightedDistributions(tree_mc, tree_sdata, whichVar,                     \
    #                     weight=self._physWeightsName, sWeight=s, rangeX=Xrange                      )       
    #     else:
    #         print 'arg=True'
    #         print 'P2VV - INFO: Comparing ' + whichVar + ' in MC before and after MC physics reweighting.'
    #         if AssymPlot:
    #             c,d = compareWeightedDistributions(tree_mc, tree_mc, whichVar,                        \
    #                     sWeight=self._physWeightsName, rangeX=Xrange, assymPlot=AssymPlot           )
    #         else: c = compareWeightedDistributions(tree_mc, tree_mc, whichVar,                        \
    #                     sWeight=self._physWeightsName, rangeX=Xrange                                )
        
    #     #file_sdata.Close()
    #     #file_mc.Close()

    #     if AssymPlot: return c,d
    #     else:         return c
