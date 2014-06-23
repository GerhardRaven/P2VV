import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-H', '--flipHadrons',  default = False,  action = 'store_true' )
parser.add_argument( '-L', '--flipLeptons',  default = False,  action = 'store_true' )
parser.add_argument( '-c', '--compare',      default = False,  action = 'store_true' )
options = parser.parse_args()

#######################################################
## configure ##
#######################################################
dataSetPath = '/data/bfys/vsyropou/Bs2JpsiKst/DiegosUnrefinedTuples/2011n.root'
dataSetName = 'DecayTree'

daughterPartNames  = dict( posHad='Kplus'    if not options.flipHadrons else 'piminus',  
                           negHad='piminus'  if not options.flipHadrons else 'Kplus',  
                           posLep='muplus'   if not options.flipLeptons else 'muminus',  
                           negLep='muminus'  if not options.flipLeptons else 'muplus'
                          )

caseSpecifier = '_%s_%s_%s_%s'%( daughterPartNames['posHad'][:3],daughterPartNames['negHad'][:3],\
                                 daughterPartNames['posLep'][:3],daughterPartNames['negLep'][:3] )

# calcualted angle names
angleNames = {}
angleNames['helcosthetaK'] = 'helcosthetaK_%s' %caseSpecifier
angleNames['helcosthetaL'] = 'helcosthetaL_%s'%caseSpecifier
angleNames['helphi']       = 'helphi_%s'%caseSpecifier

oldangleNames['helcosthetaK'] = 'helcosthetaK'
oldangleNames['helcosthetaL'] = 'helcosthetaL'
oldangleNames['helphi']       = 'B0_Phi'

##################################
## calculate helicity angles ##
####################################

# open input file
f = TFile.Open(dataSetPath)
t = f.Get(dataSetName)
t.SetBranchStatus('*',0)
for name in ['helcosthetaL','helcosthetaK','B0_Phi']: 
    t.SetBranchStatus(name,1)
for name in [ '%s_P%s' % ( part, comp ) for part in [ 'muplus', 'muminus', 'Kplus', 'piminus' ] for comp in ( 'X', 'Y', 'Z' ) ]: 
    t.SetBranchStatus(name,1)

# intermediate file
tempFile = TFile.Open(dataSetPath[-10:-5] + '_%s.root'%caseSpecifier,'recreate')
tree = t.CloneTree()

# close initial file
f.Close()
del f

# masses
from P2VV.Load import P2VVLibrary
from ROOT import TDatabasePDG, TFile, addHelicityAnglesToTree, TCanvas, TH1D
from math import pi

MeV = 1000 # TDatabasePDG is in GeV
PDG = TDatabasePDG()
Mmu = PDG.GetParticle('mu-').Mass()*MeV
Mk  = PDG.GetParticle('K-').Mass()*MeV
Mpi = PDG.GetParticle('pi-').Mass()*MeV

posHadMass = Mk  if 'K' in daughterPartNames['posHad'] else Mpi
negHadMass = Mpi if 'pi' in daughterPartNames['negHad'] else Mk
lepMass    = Mmu

print ' P2VV - INFO: The following associations will be made:\n '\
    ' Positive hadron name: %s, mass=%s \n  Negative hadron name: %s, mass=%s \n '\
    ' Positive lepton name: %s, mass=%s \n  Negative lepton name: %s, mass=%s \n'\
    ' Units MUST be in GeV. Check!!!'\
    %(daughterPartNames['posHad'], posHadMass, daughterPartNames['negHad'], negHadMass,\
      daughterPartNames['posLep'], lepMass,    daughterPartNames['negLep'], lepMass)
print ' Helicity angles names:\n helcosthetaK = %s \n helcosthetaL = %s \n helphi = %s'\
    %( angleNames['helcosthetaK'], angleNames['helcosthetaL'], angleNames['helphi'] )

addHelicityAnglesToTree(tree, 
                        daughterPartNames['posHad'], daughterPartNames['negHad'], 
                        daughterPartNames['posLep'], daughterPartNames['negLep'],
                        posHadMass, negHadMass, lepMass, lepMass,
                        angleNames['helcosthetaK'], angleNames['helcosthetaL'], angleNames['helphi']
                        )

# close outfile
tempFile.cd()
tree.Write()
tree.Show()
tempFile.Close()
del tempFile, tree
print 'P2VV - INFO: Wrote file: %s'%dataSetPath[-10:-5] + '_%s.root'%caseSpecifier

# re-open outfile
file_ = TFile.Open(dataSetPath[-10:-5] + '_%s.root'%caseSpecifier)
tree = file_.Get('DecayTree')

h_cthK = TH1D('cthK','cthK',100,-1,1)
h_cthL = TH1D('cthL','cthL',100,-1,1)
h_phi  = TH1D('phi','phi',100,-pi-.5,pi+.5)

h_my_cthK = TH1D('mycthK','mycthK',100,-1,1)
h_my_cthL = TH1D('mycthL','mycthL',100,-1,1)
h_my_phi  = TH1D('myphi','myphi',100,-pi-.5,pi+.5)

h_pull_cthK = TH1D('cThKpull','cThKpull', 100, -1e-6, 1e-6)
h_pull_cthL = TH1D('cThLpull','cThLpull', 100, -1e-6, 1e-6)
h_pull_phi  = TH1D('phipull','phipull', 100, -5e-6, 5e-6)

# new comparision
for entry in tree:
    h_my_cthK.Fill(getattr(entry,angleNames['helcosthetaK']))
    h_my_cthL.Fill(getattr(entry,angleNames['helcosthetaL']))
    h_my_phi.Fill( getattr(entry,angleNames['helphi']      ))       
    
    h_cthK.Fill(getattr(entry,oldangleNames['helcosthetaK'])) 
    h_cthL.Fill(getattr(entry,oldangleNames['helcosthetaL']))
    h_phi.Fill(getattr(entry,oldangleNames['helphi']))
    
    h_pull_cthK.Fill(getattr(entry,angleNames['helcosthetaK'])    - getattr(entry,oldangleNames['helcosthetaK']) )
    h_pull_cthL.Fill(getattr(entry,angleNames['helcosthetaL'])    - getattr(entry,oldangleNames['helcosthetaL']) )
    h_pull_phi.Fill (getattr(entry,angleNames['helphi'])          - getattr(entry,oldangleNames['helphi']) )
    
# compare
if options.compare:
    c = TCanvas('calculated angles','calculated angles')
    c.Divide(3,2)
    
    c.cd(1)
    h_cthK.Draw()
    
    c.cd(2)
    h_cthL.Draw()
    
    c.cd(3)
    h_phi.Draw()
    
    c.cd(4)
    h_my_cthK.Draw()
    
    c.cd(5)
    h_my_cthL.Draw()
    
    c.cd(6)
    h_my_phi.Draw()
    
    
    # plot compare
    c2 = TCanvas('pulls','pulls')
    c2.Divide(2,2)
    
    c2.cd(1)
    h_pull_cthK.Draw()
    
    c2.cd(2)
    h_pull_cthL.Draw()
    
    c2.cd(3)
    h_pull_phi.Draw()
