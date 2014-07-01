import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-H', '--flipHadrons',      default = False,  action = 'store_true' )
parser.add_argument( '-L', '--flipLeptons',      default = False,  action = 'store_true' )
parser.add_argument( '-c', '--compare',          default = False,  action = 'store_true' )
parser.add_argument( '-m', '--minimalTuples',    default = False,  action = 'store_true' )
parser.add_argument( '-a', '--addAllConfigrtns', default = False,  action = 'store_true' )
parser.add_argument( '-f', '--KstFlavorCut',     default = ''                            )
options = parser.parse_args()
if options.addAllConfigrtns:
    options.flipHadrons = False
    options.flipLeptons = False

#######################################################
## configure ##
#######################################################
from P2VV.Load import P2VVLibrary
dataSetPath = '/home/vilain/Desktop/workDir/data/storage/data/Bd2JpsiKstar_DTT_after_yuehongs_script_20120203.root'
# '~/Desktop/workDir/data/storage/data/2011p.root' # '~/Desktop/workDir/data/storage/data/Bd2JpsiKstar_DTT_after_yuehongs_script_20120203.root'
dataSetName = 'DecayTree'

daughterPartNames  = dict( posHad='Kplus'    if not options.flipHadrons else 'piminus',  
                           negHad='piminus'  if not options.flipHadrons else 'Kplus',  
                           posLep='muplus'   if not options.flipLeptons else 'muminus',  
                           negLep='muminus'  if not options.flipLeptons else 'muplus'
                          )

caseSpecifier = '_%s_%s_%s_%s'%( daughterPartNames['posHad'][:3],daughterPartNames['negHad'][:3],\
                                 daughterPartNames['posLep'][:3],daughterPartNames['negLep'][:3] )

# calcualted angle names
angleNames, oldangleNames = {},{}
angleNames['helcosthetaK'] = 'helcosthetaK_%s' %caseSpecifier
angleNames['helcosthetaL'] = 'helcosthetaL_%s'%caseSpecifier
angleNames['helphi']       = 'helphi_%s'%caseSpecifier

oldangleNames['helcosthetaK'] = 'B0_ThetaK'
oldangleNames['helcosthetaL'] = 'B0_ThetaL'
oldangleNames['helphi']       = 'B0_Phi'

minimalNtuple = options.minimalTuples
if options.KstFlavorCut=='neg':
    print 'P2VV - INFO: Selecting negative Kaons only.'
    selectionString = 'Kplus_ID == -321' 
elif options.KstFlavorCut=='pos':
    selectionString = 'Kplus_ID == 321' 
    print 'P2VV - INFO: Selecting positive Kaons only.'
else: selectionString = ''
# else: assert False, 'P2VV - INFO: Provide valid argument (-f) for Kaon sign {pos/neg}.'

##################################
## calculate helicity angles ##
####################################

# open input file
from ROOT import TFile
f = TFile.Open(dataSetPath)
t = f.Get(dataSetName)

# create intermediate file
tempFile = TFile.Open(dataSetPath[-10:-5] + '_%s.root'%caseSpecifier,'recreate')
if selectionString:
    print 'P2VV - INFO: Applying the following cuts %s. Initial entries: %s'%(selectionString,t.GetEntries())
    tree = t.CopyTree(selectionString)
    print 'P2VV - INFO: Entries after cuts: %s'%tree.GetEntries()
else:tree = t.CloneTree() 

# switch off unncessessary branches
if minimalNtuple:
    tree.SetBranchStatus('*',0)
    for name in [oldangleNames['helcosthetaK'], oldangleNames['helcosthetaL'], oldangleNames['helphi']]: 
        tree.SetBranchStatus(name,1)
    for name in [ '%s_P%s' % ( part, comp ) for part in [ 'muplus', 'muminus', 'Kplus', 'piminus' ] for comp in ( 'X', 'Y', 'Z' ) ]: 
        tree.SetBranchStatus(name,1)

# close initial file
f.Close()
del f

# import stuff
from ROOT import TDatabasePDG, addHelicityAnglesToTree
from math import pi

# masses
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
    ' Units MUST be in MeV. Check!!!'\
    %(daughterPartNames['posHad'], posHadMass, daughterPartNames['negHad'], negHadMass,\
      daughterPartNames['posLep'], lepMass,    daughterPartNames['negLep'], lepMass)
print ' Helicity angles names:\n helcosthetaK = %s \n helcosthetaL = %s \n helphi = %s'\
    %( angleNames['helcosthetaK'], angleNames['helcosthetaL'], angleNames['helphi'] )

addHelicityAnglesToTree(tree, 
                        daughterPartNames['posHad'], daughterPartNames['negHad'], 
                        daughterPartNames['posLep'], daughterPartNames['negLep'],
                        posHadMass, negHadMass, lepMass, lepMass,
                        angleNames['helcosthetaK'], angleNames['helcosthetaL'], angleNames['helphi'],
                        'Float_t')


if options.addAllConfigrtns:
    print 'P2VV - INFO: Fliping hadrons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['negHad'], daughterPartNames['posHad'], 
                            daughterPartNames['posLep'], daughterPartNames['negLep'],
                            negHadMass, posHadMass, lepMass, lepMass,
                            'helcosthetaK__pim_Kpl_mup_mum','helcosthetaL__pim_Kpl_mup_mum','helphi__pim_Kpl_mup_mum',
                            'Float_t')

    print 'P2VV - INFO: Fliping leptons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['posHad'], daughterPartNames['negHad'] , 
                            daughterPartNames['negLep'], daughterPartNames['posLep'], 
                            posHadMass, negHadMass, lepMass, lepMass,
                            'helcosthetaK__Kpl_pim_mum_mup','helcosthetaL__Kpl_pim_mum_mup','helphi__Kpl_pim_mum_mup',
                            'Float_t')

    print 'P2VV - INFO: Fliping hadrons and leptons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['negHad'], daughterPartNames['posHad'], 
                            daughterPartNames['negLep'], daughterPartNames['posLep'], 
                            negHadMass, posHadMass, lepMass, lepMass,
                            'helcosthetaK__pim_Kpl_mum_mup','helcosthetaL__pim_Kpl_mum_mup','helphi__pim_Kpl_mum_mup',
                            'Float_t')

# close outfile
tempFile.cd()
tree.Write()
tree.Show()
tempFile.Close()
del tempFile, tree
print 'P2VV - INFO: Wrote file: %s'%dataSetPath[-10:-5] + '_%s.root'%caseSpecifier

if options.compare:
    # re-open outfile
    file_ = TFile.Open(dataSetPath[-10:-5] + '_%s.root'%caseSpecifier)
    tree = file_.Get('DecayTree')
    
    from math import cos 
    from ROOT import TH1D, TH2D

    h_cthK = TH1D('cthK','cthK',100,-1,1)
    h_cthL = TH1D('cthL','cthL',100,-1,1)
    h_phi  = TH1D('phi','phi',100,-pi-.5,pi+.5)
    
    h_my_cthK = TH1D('mycthK','mycthK',100,-1,1)
    h_my_cthL = TH1D('mycthL','mycthL',100,-1,1)
    h_my_phi  = TH1D('myphi','myphi',100,-pi-.5,pi+.5)
    
    # h_pull_cthK = TH1D('cThKpull','cThKpull', 100, -1e-6, 1e-6)
    # h_pull_cthL = TH1D('cThLpull','cThLpull', 100, -1e-6, 1e-6)
    # h_pull_phi  = TH1D('phipull','phipull',   100, -1e-6, 1e-6)

    h_scat_cthK = TH2D('cThKscat','cThKscat', 100, -1, 1,100, -1, 1)
    h_scat_cthL = TH2D('cThLscat','cThLscat', 100, -1, 1,100, -1, 1)
    h_scat_phi  = TH2D('phiscat','phiscat',   100, -pi, pi, 100, -pi, pi)

    # new comparision
    for entry in tree:
        h_my_cthK.Fill(getattr(entry,angleNames['helcosthetaK']))
        h_my_cthL.Fill(getattr(entry,angleNames['helcosthetaL']))
        h_my_phi.Fill( getattr(entry,angleNames['helphi']      ))       

        h_cthK.Fill(cos(getattr(entry,oldangleNames['helcosthetaK']))) 
        h_cthL.Fill(cos(getattr(entry,oldangleNames['helcosthetaL'])))
        h_phi.Fill(getattr(entry,oldangleNames['helphi']))

        h_scat_cthK.Fill(getattr(entry,angleNames['helcosthetaK']),   cos(getattr(entry,oldangleNames['helcosthetaK'])))
        h_scat_cthL.Fill(getattr(entry,angleNames['helcosthetaL']),   cos(getattr(entry,oldangleNames['helcosthetaL'])))
        h_scat_phi. Fill(getattr(entry,angleNames['helphi']      ),       getattr(entry,oldangleNames['helphi']       ))
        
        # if getattr(entry,oldangleNames['helphi']) < 0: 
        #            h_scat_phi. Fill(angleNames['helphi'], getattr(entry,oldangleNames['helphi']) - pi)
        # else: h_scat_phi. Fill(angleNames['helphi'], getattr(entry,oldangleNames['helphi']) + pi)
        

           
        # h_cthK.Fill(getattr(entry,oldangleNames['helcosthetaK'])) 
        # h_cthL.Fill(getattr(entry,oldangleNames['helcosthetaL']))
        # h_phi.Fill(getattr(entry,oldangleNames['helphi']))
        
        # h_pull_cthK.Fill(getattr(entry,angleNames['helcosthetaK'])    - getattr(entry,oldangleNames['helcosthetaK']) )
        # h_pull_cthL.Fill(getattr(entry,angleNames['helcosthetaL'])    - getattr(entry,oldangleNames['helcosthetaL']) )
        # h_pull_phi.Fill (getattr(entry,angleNames['helphi'])          - getattr(entry,oldangleNames['helphi']) )
        
# compare
if options.compare:
    from ROOT import TCanvas
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
    
    # scater plots
    c3 = TCanvas('scatter','scatter')
    c3.Divide(2,2)
    c3.cd(1)
    tree.Draw('%s:%s'%(oldangleNames['helcosthetaK'],angleNames['helcosthetaK']))
    c3.cd(2)
    tree.Draw('%s:%s'%(oldangleNames['helcosthetaL'],angleNames['helcosthetaL']))    
    c3.cd(3)
    tree.Draw('%s:%s'%(oldangleNames['helphi'],angleNames['helphi']))
    
    
    # pulls
    c2 = TCanvas('pulls','pulls')
    c2.Divide(2,2)
    
    c2.cd(1)
    h_pull_cthK.Draw()
    
    c2.cd(2)
    h_pull_cthL.Draw()
    
    c2.cd(3)
    h_pull_phi.Draw()

    c.Print('angles_%s%s.pdf'%(dataSetPath[-10:-5],caseSpecifier))
    c2.Print('angles_pulls_%s%s.pdf'%(dataSetPath[-10:-5],caseSpecifier))
    c3.Print('angles_scatters_%s%s.pdf'%(dataSetPath[-10:-5],caseSpecifier))
