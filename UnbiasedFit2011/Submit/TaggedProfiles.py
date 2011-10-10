########################################
### Author: Daan van Eijk
### Updated on: Jun 5 11
### Description: This script reads the root file with the workspace for the tagged fit and fits
###              The MakeProfile function is implemented in RooFitDecorators.py, but it makes a DLL profile in one go.
###              For big grids this becomes too time-consuming. So there are now scripts to make ganga jobs for profile production:
###              These are TaggedProfiles.py and SubmitTaggedProfiles.py
########################################

from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-j", type="int", dest="step")
parser.add_option("-n", type="int", dest="npoints")
(options, args) = parser.parse_args()

_step = options.step
_npoints = options.npoints
print 'step =', _step
print 'npoints =', _npoints

#import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
#myStyle = rootStyle.plainstyle()
#gROOT.SetStyle(myStyle.GetName())
#gROOT.ForceStyle()
#gStyle.UseCurrentStyle()

angcorr = False
taggingsyst = True
blinded = True
phisparam = True
splitmKK = False

name = 'BsJpsiPhi2011_fit_phis'

wsfile = TFile('TaggedWS.root')
ws = wsfile.Get('ws')

if taggingsyst:
    if not angcorr:
        pdfbeforeblinding = ws.pdf('pdf_ext_inc_tag_syst')
    if angcorr:
        pdfbeforeblinding = ws.pdf('pdf_ext_angcorr_inc_tag_syst')
    
else:
    pdfbeforeblinding = ws.pdf('pdf_ext_angcorr')

################
### Blinding ###
################

if blinded:
    #Building blinded parameters
    if phisparam:
        ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
        ws.factory("RooUnblindUniform::phis_blind('BsCustard',3.,phis)")    

        #calling RooCustomizer
        customizer = RooCustomizer(pdfbeforeblinding,'blinded')
        customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
        customizer.replaceArg( ws['phis'], ws['phis_blind'] )

    else:
        ws.factory("RooUnblindUniform::t_sig_dG_blind('BsRooBarb',0.2,t_sig_dG)")
        ws.factory("RooUnblindUniform::S_blind('BsMickeyMouseLP2011',2.,S)")
        #ws.factory("RooUnblindUniform::D_blind('BsMiniMouseLP2011',2.,D)")
        #ws.factory("RooUnblindUniform::C_blind('BsDonaldDuckLP2011',2.,C)")

        #calling RooCustomizer
        customizer = RooCustomizer(pdfbeforeblinding,'blinded')
        customizer.replaceArg( ws['t_sig_dG'], ws['t_sig_dG_blind'] )
        customizer.replaceArg( ws['S'], ws['S_blind'] )
        #customizer.replaceArg( ws['D'], ws['D_blind'] )
        #customizer.replaceArg( ws['C'], ws['C_blind'] )

    blindedpdf = customizer.build()
    ws.put(blindedpdf)

if blinded:
    pdf = blindedpdf
else:
    pdf = pdfbeforeblinding

#ws['rs2'].setVal(0.)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.77],t_sig_dm_sigma[0.12])")

pdf = pdf
data = ws['unbiaseddata']

################
### Profiles ###
################

#setting back values
ws['#Gamma'].setVal(0.68)
ws['#Gamma'].setMax(2.)
ws['t_sig_dG'].setVal(0.060)
ws['t_sig_dG'].setMin(-0.7)
ws['t_sig_dG'].setMax(0.7)

ws['phis'].setVal(-pi)
ws['phis'].setMin(-2*pi)
ws['phis'].setMax(0)

phis = ws['phis']
deltaGamma = ws['t_sig_dG']

#From the old function
BaseNLL = RooRealVar('BaseNLL','BaseNLL',0.)
GridPointNLL = RooRealVar("GridPointNLL","GridPointNLL",0.)

stepphis = RooRealVar("stepphis","stepphis",0.)
valuephis = RooRealVar("valuephis","valuephis",0.)

stepdeltaGamma = RooRealVar("stepdeltaGamma","stepdeltaGamma",0.)
valuedeltaGamma = RooRealVar("valuedeltaGamma","valuedeltaGamma",.0)

dlogLLArgSet = RooArgSet(BaseNLL,GridPointNLL,stepphis,valuephis,stepdeltaGamma,valuedeltaGamma)
dlogLLDataSet = RooDataSet('dlogLLDataSet','dlogLLDataSet',dlogLLArgSet)

phis_min = -2*pi
phis_max = 0
phis_int = phis_max-phis_min

#get range for deltaGamma
deltaGamma_min = -0.7
deltaGamma_max = 0.7
deltaGamma_int = deltaGamma_max - deltaGamma_min

nll = pdf.createNLL(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))

#Check some inital stuff
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(True),RooFit.Minos(False),RooFit.Save(True),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))
BaseNLL.setVal(nll.getVal())
print 'BaseNLL =',BaseNLL.getVal()

for i in range(_npoints):  
    phis.setVal((phis_min+(phis_int/(2.*_npoints)))+i*(phis_int/_npoints))
    deltaGamma.setVal((deltaGamma_min+(deltaGamma_int/(2.*_npoints)))+_step*(deltaGamma_int/_npoints))
    print '***************************************************************************'
    print 'At gridpoint i = %i from %i and j = %i from %i'%(i+1,_npoints,_step+1,_npoints)
    print '%s min ='%(phis.GetName()), phis.getMin()
    print '%s max ='%(phis.GetName()), phis.getMax()
    print '%s at current gridpoint ='%(phis.GetName()), phis.getVal()
    print '%s min ='%(deltaGamma.GetName()), deltaGamma.getMin()
    print '%s max ='%(deltaGamma.GetName()), deltaGamma.getMax()
    print '%s at current gridpoint ='%(deltaGamma.GetName()), deltaGamma.getVal()
    print '***************************************************************************'
    stepphis.setVal(i+1)
    valuephis.setVal(phis.getVal())
    stepdeltaGamma.setVal(_step+1)
    valuedeltaGamma.setVal(deltaGamma.getVal())
    
    phis.setConstant(kTRUE)
    deltaGamma.setConstant(kTRUE)
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true),RooFit.ExternalConstraints(RooArgSet(ws.pdf('dmsconstraint'))))
    GridPointNLL.setVal(nll.getVal())
    print 'BaseNLL =', BaseNLL.getVal()
    print 'GridPointNLL =',GridPointNLL.getVal()

    #dlogLL.setVal(profile.getVal())
    
    dlogLLDataSet.add(dlogLLArgSet)
            
tfile = TFile('profile.root','RECREATE')
dlogLLDataSet.Write() 
tfile.Close()

