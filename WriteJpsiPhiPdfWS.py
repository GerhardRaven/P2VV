##############################################################
##############################################################
###                                                        ###
### Basic script to make Jpsiphi-pdf,                      ###
### generate toy data, save data, fit and plot             ###
###                                                        ###
### RooDataSet generated here can be used in other scripts ###
###                                                        ###
### Tristan du Pree                                        ###
###                                                        ###
##############################################################
##############################################################

from ROOT import *
import GaudiPython
P2VV = GaudiPython.gbl.P2VV
#to load functions (made with namespace function) like makePVVPdf:
GaudiPython.loaddict('P2VVDict')
import math
from math import pi

#############################################
### Here we go...your favourite settings  ###
#############################################

read      = False
filename  = "testdata.root"
_nJpsiphi = 800

##########################
### physics parameters ###
##########################

#_Phi_s  = -0.04
_Phi_s  = 0.

#_tau_s  =  1.47
#That defines gamma:
_G_s = 0.68

_dG_s   =  0.075

_dm_s   = 17.7
#_dGG_s  =  0.1
#_dGG_s  = _dG_s*_tau_s

_Azero_jpsiphi  = 0.524
_Aperp_jpsiphi  = 0.231

_phi_perp_jpsiphi = 2.91
_phi_par_jpsiphi  =  -2.93
_phi_zero_jpsiphi =  0.

########################
### Define workspace ###
########################

ws = RooWorkspace("ws")

###################
### Observables ###
###################

#t      = RooRealVar("t",     "time",        -2., 20. )
t      = RooRealVar("t",     "time",        0.3, 12. )

#todo: make these variables in the constructor of the xAngles
#ctheta = RooRealVar("ctheta","cos(#theta)", -1.0,  1.0)
#cpsi   = RooRealVar("cpsi",  "cos(#psi)",   -1.0,  1.0)
#phi    = RooRealVar("phi",   "#phi",        -pi ,  pi )

#The names of these RooRealVars are inspired by the common names in the NTuple....
ctheta = RooRealVar("trcostheta","cos(#theta)", -1.0,  1.0)
cpsi   = RooRealVar("trcospsi",  "cos(#psi)",   -1.0,  1.0)
phi    = RooRealVar("trphi",   "#phi",        -pi ,  pi )

myThreeAngles = P2VV.ThreeAngles(ctheta,cpsi,phi)

############ organizing the flavour tagging categories ############ 

### to generate
# bsjpsiphi (when you generate Jpsiphi, note: name = "state")
#The names of this RooCategory are inspired by the common names in the NTuple....
jpsiphiState = RooCategory("tagdecision","state Jpsiphi")
jpsiphiState.defineType('Bs_Jpsiphi',1)
jpsiphiState.defineType('Bsbar_Jpsiphi',-1)
jpsiphiState.defineType('untagged',0)

#jpsiphiState = RooCategory("jpsiphiState","state Jpsiphi")
#jpsiphiState.defineType("Bs_Jpsiphi")
#jpsiphiState.defineType("Bsbar_Jpsiphi")

#### and now map the states to the real categories in the fit pdf ####

#Flavor initial state
tagFlavCatJpsiphi = RooMappedCategory("tagFlavCatJpsiphi","tagFlavCatJpsiphi",jpsiphiState)
tagFlavCatJpsiphi.map("Bs_*",    "Bstag",   +1)
tagFlavCatJpsiphi.map("Bsbar_*", "Bsbartag",-1)

#Flavor Final state
recFlavCatJpsiphi = RooMappedCategory("recFlavCatJpsiphi","recFlavCatJpsiphi",jpsiphiState,"*_Jpsiphi",0)

#### finally also define some ranges (for practical purposes like plotting) ####
jpsiphiState.setRange("BsJpsiphi",   "Bs_Jpsiphi")
jpsiphiState.setRange("BsbarJpsiphi","Bsbar_Jpsiphi")

###########################
### Sets of Observables ###
###########################

#to generate data 
dataObsSetJpsiphi = RooArgSet(t,
                              myThreeAngles.ctheta(),
                              myThreeAngles.cpsi(),
                              myThreeAngles.phi(),
                              jpsiphiState)


####################
### Ang eff corr ###
####################

#todo: make these variables in the constructor of the list
phi1_jpsiphi = RooRealVar("phi1_jpsiphi", "phi1_jpsiphi", 32./3.*pi/3.)
phi2_jpsiphi = RooRealVar("phi2_jpsiphi", "phi2_jpsiphi", 32./3.*pi/3.)
phi3_jpsiphi = RooRealVar("phi3_jpsiphi", "phi3_jpsiphi", 32./3.*pi/3.)
phi4_jpsiphi = RooRealVar("phi4_jpsiphi", "phi4_jpsiphi",           0.)
phi5_jpsiphi = RooRealVar("phi5_jpsiphi", "phi5_jpsiphi",           0.)
phi6_jpsiphi = RooRealVar("phi6_jpsiphi", "phi6_jpsiphi",           0.)

myAngularAcceptanceListJpsiphi = P2VV.THREEANGLES.FFSS.AngularAcceptanceCorrectionList(phi1_jpsiphi,
                                                                                       phi2_jpsiphi,
                                                                                       phi3_jpsiphi,
                                                                                       phi4_jpsiphi,
                                                                                       phi5_jpsiphi,
                                                                                       phi6_jpsiphi)

#####################
### amps & phases ###
#####################

#you can (should) fix the sum to one
Azero_sq_jpsiphi = RooRealVar("Azero_sq_jpsiphi",  "Azero_sq_jpsiphi", _Azero_jpsiphi, 0.3, 0.6)
Aperp_sq_jpsiphi = RooRealVar("Aperp_sq_jpsiphi",  "Aperp_sq_jpsiphi", _Aperp_jpsiphi, 0.2, 0.5)
Apar_sq_jpsiphi  = RooFormulaVar("Apar_sq_jpsiphi","Apar_sq_jpsiphi", "1-@0-@1",RooArgList(Azero_sq_jpsiphi,Aperp_sq_jpsiphi))
#Aperp_sq = RooRealVar("Aperp_sq","Aperp_sq", 0.2, 0., 1.)

#you can constrain and e.g. fix one of the phases (e.g. phi_zero) to zero
#phi_perp_jpsiphi = RooRealVar("phi_perp_jpsiphi","phi_perp_jpsiphi", _phi_perp_jpsiphi, 0., 2.*pi)
#phi_par_jpsiphi  = RooRealVar("phi_par_jpsiphi", "phi_par_jpsiphi",  _phi_par_jpsiphi,  0., 2.*pi)
#phi_zero_jpsiphi = RooRealVar("phi_zero_jpsiphi","phi_zero_jpsiphi", _phi_zero_jpsiphi,  0., 2.*pi)

phi_perp_jpsiphi = RooRealVar("phi_perp_jpsiphi","phi_perp_jpsiphi", _phi_perp_jpsiphi)
phi_par_jpsiphi  = RooRealVar("phi_par_jpsiphi", "phi_par_jpsiphi",  _phi_par_jpsiphi)
phi_zero_jpsiphi = RooRealVar("phi_zero_jpsiphi","phi_zero_jpsiphi", _phi_zero_jpsiphi)

###################################
### Strong phases or parameters ###
###################################

Aperp_Azero_cos_jpsiphi = RooFormulaVar("Aperp_Azero_cos_jpsiphi","sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Aperp_sq_jpsiphi,phi_perp_jpsiphi,phi_zero_jpsiphi))
Aperp_Azero_sin_jpsiphi = RooFormulaVar("Aperp_Azero_sin_jpsiphi","sqrt(@0)*sqrt(@1)*sin(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Aperp_sq_jpsiphi,phi_perp_jpsiphi,phi_zero_jpsiphi))
Azero_Apar_cos_jpsiphi  = RooFormulaVar("Azero_Apar_cos_jpsiphi", "sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Azero_sq_jpsiphi, Apar_sq_jpsiphi, phi_par_jpsiphi ,phi_zero_jpsiphi))
Aperp_Apar_cos_jpsiphi  = RooFormulaVar("Aperp_Apar_cos_jpsiphi", "sqrt(@0)*sqrt(@1)*cos(@2-@3)",
                                        RooArgList(Apar_sq_jpsiphi, Aperp_sq_jpsiphi, phi_perp_jpsiphi,phi_par_jpsiphi))
Aperp_Apar_sin_jpsiphi  = RooFormulaVar("Aperp_Apar_sin_jpsiphi", "sqrt(@0)*sqrt(@1)*sin(@2-@3)",
                                        RooArgList(Apar_sq_jpsiphi, Aperp_sq_jpsiphi, phi_perp_jpsiphi,phi_par_jpsiphi)) 
    
###################################
### Amplitudes and phases lists ###
###################################

myAmplitudesAndPhasesListJpsiphi = P2VV.AmplitudesAndPhasesList(Azero_sq_jpsiphi,
                                                                Apar_sq_jpsiphi,
                                                                Aperp_sq_jpsiphi,
                                                                Aperp_Azero_cos_jpsiphi,
                                                                Aperp_Azero_sin_jpsiphi,
                                                                Azero_Apar_cos_jpsiphi,
                                                                Aperp_Apar_cos_jpsiphi,
                                                                Aperp_Apar_sin_jpsiphi)

##########################
### Physics parameters ###
##########################

G_s = RooRealVar('G_s','G_s',_G_s, 0.4, 0.9)

tau_s  = RooFormulaVar("tau_s","#tau_{s}","1/@0",RooArgList(G_s))

dG_s = RooRealVar('dG_s','dG_s',_dG_s, -0.3,0.3)

dGG_s  = RooFormulaVar("dGG_s","#Delta#Gamma_{s}/#Gamma_{s}","@0/@1",RooArgList(dG_s,G_s))

dm_s   = RooRealVar("dm_s", "#Delta m_{s}",
                    _dm_s)#, 15., 25.)

#Phi_s = RooRealVar("Phi_s","Phi_s",_Phi_s, -pi, pi)
Phi_s = RooRealVar("Phi_s","Phi_s",_Phi_s)

Sf_s   =  RooFormulaVar("Sf_s", "sin(Phi_s)", "sin(@0)", RooArgList(Phi_s))
Df_s   =  RooFormulaVar("Df_s", "cos(Phi_s)", "cos(@0)", RooArgList(Phi_s))

################
### Blinding ###
################

## G_s_blindstring = 'BsCalvin'
## dG_s_blindstring = 'BsHobbes'
## Phi_s_blindstring = 'BsBsGoofy'
## Sf_s_blindstring = 'BsMickeyMouse'
## Df_s_blindstring = 'BsMiniMouse'


## G_s_blind = RooUnblindUniform('G_s_blind','G_s_blind',G_s_blindstring,0.4,G_s)
## dG_s_blind = RooUnblindUniform('dG_s_blind','dG_s_blind',dG_s_blindstring,0.2,dG_s)
## #Phi_s_blind = RooUnblindUniform('Phi_s_blind','Phi_s_blind',Phi_s_blindstring,3,Phi_s)
## Sf_s_blind = RooUnblindUniform('Sf_s_blind','Sf_s_blind',Sf_s_blindstring,2,Sf_s)
## Df_s_blind = RooUnblindUniform('Df_s_blind','Df_s_blind',Df_s_blindstring,2,Df_s)

################################
### physics parameters lists ###
################################

myPhysicsParametersListJpsiphi = P2VV.PhysParamList(tau_s,
                                                    dGG_s,
                                                    dm_s,
                                                    Sf_s,
                                                    Df_s)

###############################
### Experimental parameters ###
###############################

resol_t_mean  = RooRealVar("resol_t_mean",   "mean(t resol)",  0.)
resol_t_sigma = RooRealVar("resol_t_sigma",  "sigma(t resol)", 0.050)#, 0.01, 0.1 )
tres  = RooGaussModel("tres","prop time resolution", t, resol_t_mean, resol_t_sigma)

#tres  = RooTruthModel("tres","tres",t)

wmistag = RooRealVar("wmistag","mistag rate", 0.5)

##################################### building the PDFs ###################################

threeDEff = RooFormulaVar("eff_een","eff_een","1.",RooArgList())

myJpsiphiPdf_withEff = P2VV.Functions.makePVVPdf("myJpsiphiPdf_withEff",
                                                 "myJpsiphiPdf_withEff",
                                                 myThreeAngles,
                                                 t,
                                                 tagFlavCatJpsiphi,
                                                 recFlavCatJpsiphi,
                                                 myAmplitudesAndPhasesListJpsiphi,
                                                 myPhysicsParametersListJpsiphi,
                                                 tres,
                                                 wmistag,
                                                 threeDEff,
                                                 "Jpsiphi")

myJpsiphiPdf_noEff = P2VV.Functions.makePVVPdf("myJpsiphiPdf_noEff",
                                               "myJpsiphiPdf_noEff",
                                               myThreeAngles,
                                               t,
                                               tagFlavCatJpsiphi,
                                               recFlavCatJpsiphi,
                                               myAmplitudesAndPhasesListJpsiphi,
                                               myPhysicsParametersListJpsiphi,
                                               tres,
                                               wmistag,
                                               "Jpsiphi")

myJpsiphiPdf_withWeights = P2VV.Functions.makePVVPdf("myJpsiphiPdf_withWeights",
                                                     "myJpsiphiPdf_withWeights",
                                                     myThreeAngles,
                                                     t,
                                                     tagFlavCatJpsiphi,
                                                     recFlavCatJpsiphi,
                                                     myAmplitudesAndPhasesListJpsiphi,
                                                     myPhysicsParametersListJpsiphi,
                                                     tres,
                                                     wmistag,
                                                     myAngularAcceptanceListJpsiphi,
                                                     "Jpsiphi")
############################
### Writing to workspace ###
############################

getattr(ws,'import')(dataObsSetJpsiphi)
getattr(ws,'import')(tres)
getattr(ws,'import')(myJpsiphiPdf_noEff)
#getattr(ws,'import')(myJpsiphiPdf_withWeights)

file = TFile("WSJpsiPhiPdf.root","RECREATE")
ws.Write("workspace")
file.Close()


