import argparse
parser = argparse.ArgumentParser()
parser.add_argument( '-H', '--flipHadrons',      default = False,  action = 'store_true' )
parser.add_argument( '-L', '--flipLeptons',      default = False,  action = 'store_true' )
parser.add_argument( '-c', '--compare',          default = False,  action = 'store_true' )
parser.add_argument( '-m', '--minimalTuples',    default = False,  action = 'store_true' )
parser.add_argument( '-a', '--addAllConfigrtns', default = False,  action = 'store_true' )
parser.add_argument( '-r', '--reduced',          default = False,  action = 'store_true' )
parser.add_argument( '-f', '--ManualKstFlavor',  default = ''                            )
options = parser.parse_args()
if options.addAllConfigrtns:
    options.flipHadrons = False
    options.flipLeptons = False

#######################################################
## configure ##
#######################################################
from P2VV.Load import P2VVLibrary
from ROOT import std
dataSetPath   = '/data/bfys/vsyropou/Bs2JpsiKst/DiegosUnrefinedTuples/MC/' 
              # '/data/bfys/vsyropou/Bs2JpsiKst/edinburgTuples/'  
dataSetFile   = 'Bd_MCT_2014_fuckingAngles3_costhetaKrevenge.root'
              # 'Bd2JpsiKstar_DTT_after_yuehongs_script_20120203.root' 
dataSetName   = 'T' #'DecayTree'
minimalNtuple = options.minimalTuples # switch of unused branches

# branch names as they appear in the input tree
daughterPartNames = dict(  posHad = 'k1'  if not options.flipHadrons else 'pi1' , # 'Kplus'   if not options.flipHadrons else 'piminus', # 'k1'  if not options.flipHadrons else 'p1' ,  
                           negHad = 'pi1' if not options.flipHadrons else 'k1' , # 'piminus' if not options.flipHadrons else 'Kplus',   # 'p1'  if not options.flipHadrons else 'k1' , 
                           posLep = 'mu1' if not options.flipLeptons else 'mu2', # 'muplus'  if not options.flipLeptons else 'muminus', # 'mu1' if not options.flipLeptons else 'mu2', 
                           negLep = 'mu2' if not options.flipLeptons else 'mu1', # 'muminus' if not options.flipLeptons else 'muplus'   # 'mu2' if not options.flipLeptons else 'mu1'  
                          )

# sufixes for the components of the base momenta name 
daughterPartNameSufixes = dict(x='p1', y='p2', z='p3') # dict(x='_PX', y='_PY', z='_PZ')  # 
sufixes = std.vector('TString')()
for component in ['x','y','z']: sufixes.push_back(daughterPartNameSufixes[component])

# masses of daughters
daughterPartMasses = dict( posHad = 'kaon' if not options.flipHadrons else 'pion',
                           negHad = 'pion' if not options.flipHadrons else 'kaon',
                           posLep = 'muon', 
                           negLep = 'muon' 
                         )
# kaon id branch name
kaonIDname = 'k1ID' # 'Kplus_ID'

# automatic bookeeping specifier
caseSpecifier = '_%s_%s_%s_%s'%( daughterPartNames['posHad'][:3],daughterPartNames['negHad'][:3],\
                                 daughterPartNames['posLep'][:3],daughterPartNames['negLep'][:3] )

# new angle names
angleNames, oldangleNames = {},{}
angleNames['helcosthetaK'] = 'helcosthetaK_%s' %caseSpecifier
angleNames['helcosthetaL'] = 'helcosthetaL_%s'%caseSpecifier
angleNames['helphi']       = 'helphi_%s'%caseSpecifier

# old angle names
oldangleNames['helcosthetaK'] =  'helcosthetaK' # 'B0_ThetaK' #'cK'     
oldangleNames['helcosthetaL'] =  'helcosthetaL' # 'B0_ThetaL' #'cL'     
oldangleNames['helphi']       =  'helphi'       # 'B0_Phi'    #'ph_old' 
rawAngles = False # if true, it applies cos(theta{K,L}) when filling comaprision histograms

# manually select Kst flavour, if datasets are not already splited.
if options.ManualKstFlavor=='neg':
    print 'P2VV - INFO: Selecting negative Kaons only.'
    selectionString = kaonIDname + ' == -321' 
elif options.ManualKstFlavor=='pos':
    selectionString = kaonIDname + '  == 321' 
    print 'P2VV - INFO: Selecting positive Kaons only.'
else: 
    selectionString = ''
    print 'P2VV - INFO: Assuming that dataset is already split for pos/neg Kaon. No cut applied.'
if   options.reduced and selectionString:     selectionString += ' && runNumber < 92e3'
elif options.reduced and not selectionString: selectionString += 'runNumber < 92e3'


##################################
## calculate helicity angles ##
####################################

# open input file
from ROOT import TFile
f = TFile.Open(dataSetPath + dataSetFile)
t = f.Get(dataSetName)

# create intermediate file
tempFile = TFile.Open(dataSetFile[:-5] + '%s_%s.root'%(caseSpecifier,options.ManualKstFlavor + 'Kaons'),'recreate')
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

if   daughterPartMasses['posHad'] == 'kaon': posHadMass = Mk  
elif daughterPartMasses['posHad'] == 'pion': posHadMass = Mpi
else: assert False, 'P2VV - ERROR: Cannot assign mass to positive hadron.'
if   daughterPartMasses['negHad'] == 'kaon': negHadMass = Mk 
elif daughterPartMasses['negHad'] == 'pion': negHadMass = Mpi 
else: assert False, 'P2VV - ERROR: Cannot assign mass to negative hadron.'
if daughterPartMasses['posLep'] == 'muon': lepMass = Mmu
else: assert False, 'P2VV - ERROR: Cannot assign mass to possitve hadron.'
if daughterPartMasses['negLep'] == 'muon': lepMass = Mmu
else: assert False, 'P2VV - ERROR: Cannot assign mass to negative hadron.'

# add helicity angles to tree
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
                        sufixes)

# add all the possible configurations as well  
if options.addAllConfigrtns:
    print 'P2VV - INFO: Fliping hadrons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['negHad'], daughterPartNames['posHad'], 
                            daughterPartNames['posLep'], daughterPartNames['negLep'],
                            negHadMass, posHadMass, lepMass, lepMass,
                            'helcosthetaK__pim_Kpl_mup_mum','helcosthetaL__pim_Kpl_mup_mum','helphi__pim_Kpl_mup_mum',
                            sufixes)

    print 'P2VV - INFO: Fliping leptons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['posHad'], daughterPartNames['negHad'] , 
                            daughterPartNames['negLep'], daughterPartNames['posLep'], 
                            posHadMass, negHadMass, lepMass, lepMass,
                            'helcosthetaK__Kpl_pim_mum_mup','helcosthetaL__Kpl_pim_mum_mup','helphi__Kpl_pim_mum_mup',
                            sufixes)

    print 'P2VV - INFO: Fliping hadrons and leptons and recalculating helicity angles'
    addHelicityAnglesToTree(tree, 
                            daughterPartNames['negHad'], daughterPartNames['posHad'], 
                            daughterPartNames['negLep'], daughterPartNames['posLep'],
                            negHadMass, posHadMass, lepMass, lepMass,
                            'helcosthetaK__pim_Kpl_mum_mup','helcosthetaL__pim_Kpl_mum_mup','helphi__pim_Kpl_mum_mup',
                            sufixes)

# close outfile
tempFile.cd()
tree.Write()
tree.Show()
tempFile.Close()
del tempFile, tree
print 'P2VV - INFO: Wrote file: %s'%dataSetFile[:-5] + '%s.root'%caseSpecifier

if options.compare:
    # re-open outfile
    file_ = TFile.Open(dataSetFile[:-5] + '_%s.root'%caseSpecifier)
    tree = file_.Get(dataSetName)
    
    from math import cos 
    from ROOT import TH1D, TH2D, TCanvas

    h_cthK = TH1D('cthK','cthK',100,-1,1)
    h_cthL = TH1D('cthL','cthL',100,-1,1)
    h_phi  = TH1D('phi','phi',100,-pi-.5,pi+.5)
    
    h_my_cthK = TH1D('mycthK','mycthK',100,-1,1)
    h_my_cthL = TH1D('mycthL','mycthL',100,-1,1)
    h_my_phi  = TH1D('myphi','myphi',100,-pi-.5,pi+.5)

    h_scat_cthK = TH2D('helcosthetaK', 'helcosthetaK', 100, -1, 1,   100, -1, 1)
    h_scat_cthL = TH2D('helcosthetaL', 'helcosthetaL', 100, -1, 1,   100, -1, 1)
    h_scat_phi  = TH2D('helphi',       'helphi',       100, -pi, pi, 100, -pi, pi)

    # value getter function to choose between raw angle and cos(rawAngles)
    _val = lambda ent,name,flag: cos(getattr(ent,name)) if flag else getattr(ent,name)

    # new comparision
    for entry in tree:
        h_my_cthK.Fill(getattr(entry,angleNames['helcosthetaK']))
        h_my_cthL.Fill(getattr(entry,angleNames['helcosthetaL']))
        h_my_phi.Fill( getattr(entry,angleNames['helphi']      ))       

        h_cthK.Fill( _val(entry,oldangleNames['helcosthetaK'],rawAngles) ) 
        h_cthL.Fill( _val(entry,oldangleNames['helcosthetaL'],rawAngles) )
        h_phi.Fill(getattr(entry,oldangleNames['helphi']))

        h_scat_cthK.Fill( _val(entry,oldangleNames['helcosthetaK'],rawAngles), getattr(entry,angleNames['helcosthetaK'])  )
        h_scat_cthL.Fill( _val(entry,oldangleNames['helcosthetaL'],rawAngles), getattr(entry,angleNames['helcosthetaL'])  )
        h_scat_phi. Fill(      getattr(entry,oldangleNames['helphi'] )       , getattr(entry,angleNames['helphi']      ) )      
     
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
    for hist in [h_scat_cthK,h_scat_cthL,h_scat_phi]:
        hist.SetStats(0)
        hist.SetXTitle('#varphi')
        hist.SetYTitle('#varphi^{vasilis}')
        hist.GetXaxis().SetTitleOffset(0.4)
        hist.GetXaxis().SetTitleSize(0.08)
        hist.GetYaxis().SetTitleOffset(0.5)
        hist.GetYaxis().SetTitleSize(0.08)


    c3 = TCanvas('scatter','scatter')
    c3.Divide(3,1)
    c3.cd(1)
    h_scat_cthK.Draw()
    # tree.Draw('TMath::Cos(%s):%s'%(oldangleNames['helcosthetaK'],angleNames['helcosthetaK']))
    c3.cd(2)
    h_scat_cthL.Draw()
    # tree.Draw('TMath::Cos(%s):%s'%(oldangleNames['helcosthetaL'],angleNames['helcosthetaL']))    
    c3.cd(3)
    h_scat_phi.Draw()
    #tree.Draw('%s:%s'%(oldangleNames['helphi'],angleNames['helphi']))
    
    canvNameSufix = (caseSpecifier) + '_' + options.ManualKstFlavor + 'Kaons'
    c.Print('angles_%s%s.pdf'%(dataSetFile[:-5],canvNameSufix))
    c3.Print('angles_scatters_%s%s.pdf'%(dataSetFile[:-5],canvNameSufix))
