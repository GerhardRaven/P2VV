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
parser.add_option("-x", type="int", dest="phisstep")
parser.add_option("-n", type="int", dest="npoints")
(options, args) = parser.parse_args()

_phisstep = options.phisstep
_npoints = options.npoints
print 'phisstep =', _phisstep
print 'npoints = ', _npoints

angcorr = True
taggingsyst = True
blinded = True
phisparam = True

wsfile = TFile('UnbiasedWS_flatbkg.root')
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

pdf = pdf
data = ws['unbiaseddata']

assert False
####SWAVE
#ws['rs2'].setVal(0.)
#ws['rs2'].setConstant()
#ws['deltas'].setConstant()

ws['rs2'].setVal(0.03)
ws['rs2'].setMin(0.)
ws['rs2'].setMax(1.)

ws['deltas'].setVal(3.0)
ws['deltas'].setMin(-2.*pi)
ws['deltas'].setMax(2.*pi)
###########

ws['rz2'].setVal(0.485)
ws['rz2'].setMin(0.)
ws['rz2'].setMax(1.)

ws['rperp2'].setVal(0.259)
ws['rperp2'].setMin(0.)
ws['rperp2'].setMax(1.)

ws['deltaperp'].setVal(2.79)
ws['deltaperp'].setMin(-2.*pi)
ws['deltaperp'].setMax(2.*pi)

ws['deltapar'].setVal(3.24)
ws['deltapar'].setMin(-2.*pi)
ws['deltapar'].setMax(2.*pi)

ws['phis'].setVal(-1.1)
ws['phis'].setMin(-2.)
ws['phis'].setMax(0.)


################
### Profiles ###
################

#setting back values
ws['#Gamma'].setVal(0.68)
ws['#Gamma'].setMax(2.)

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

phis_min = -2.
phis_max = 0.
phis_int = float(phis_max-phis_min)

ws['phis'].setVal(-pi/2)
ws['phis'].setMin(phis_min)
ws['phis'].setMax(phis_max)

#get range for deltaGamma
deltaGamma_min = -0.1
deltaGamma_max = 0.2
deltaGamma_int = deltaGamma_max - deltaGamma_min

ws['t_sig_dG'].setVal(0.060)
ws['t_sig_dG'].setMin(deltaGamma_min)
ws['t_sig_dG'].setMax(deltaGamma_max)

nll = pdf.createNLL(data,RooFit.Extended(true),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))

#Check some inital stuff
result = pdf.fitTo(data,RooFit.Extended(true),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))
BaseNLL.setVal(nll.getVal())
print 'BaseNLL =',BaseNLL.getVal()

for i in range(_npoints):  
    phis.setVal((phis_min+(phis_int/(2.*_npoints)))+_phisstep*(phis_int/_npoints))
    deltaGamma.setVal((deltaGamma_min+(deltaGamma_int/(2.*_npoints)))+i*(deltaGamma_int/_npoints))
    print '***************************************************************************'
    print 'At gridpoint i = %i from %i and j = %i from %i'%(_phisstep+1,_npoints,i+1,_npoints)
    print '%s at current gridpoint ='%(phis.GetName()), phis.getVal()
    print '%s at current gridpoint ='%(deltaGamma.GetName()), deltaGamma.getVal()
    print '***************************************************************************'
    stepphis.setVal(_phisstep+1)
    valuephis.setVal(phis.getVal())
    stepdeltaGamma.setVal(i+1)
    valuedeltaGamma.setVal(deltaGamma.getVal())
    
    phis.setConstant(kTRUE)
    deltaGamma.setConstant(kTRUE)
    result = pdf.fitTo(data,RooFit.Extended(true),RooFit.ExternalConstraints(RooArgSet(ws['p0'],ws['p1'],ws['dmsconstraint'],ws['tres_SFconstraint'])))
    GridPointNLL.setVal(nll.getVal())
    print 'BaseNLL =', BaseNLL.getVal()
    print 'GridPointNLL =',GridPointNLL.getVal()

    dlogLLDataSet.add(dlogLLArgSet)
            
tfile = TFile('profile.root','RECREATE')
dlogLLDataSet.Write() 
tfile.Close()
