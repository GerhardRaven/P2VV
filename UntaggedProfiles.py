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

wsfile = TFile('UntaggedWS.root')
ws = wsfile.Get('ws')

angcorrpdf = ws.pdf('pdf_ext_angcorrpdf')
data = ws.data('data')

reblinded = True

if reblinded:
    ws.factory("RooUnblindUniform::phis_blind('BsOlivier',3,phis)")

    customizer = RooCustomizer(angcorrpdf,'reblinded')

    phis_unblind = ws.var('phis')
    phis_blind = ws.function('phis_blind')

    customizer.replaceArg( phis_unblind, phis_blind )

    reblindedpdf = customizer.build()
    ws.put(reblindedpdf)

if reblinded:
    pdf = ws.pdf('pdf_ext_angcorrpdf_reblinded')
else:
    pdf = ws.pdf('pdf_ext_angcorrpdf')

#setting back values
ws.var('wtag').setVal(0.5)
ws.var('#Gamma').setVal(0.68)
ws.var('#Gamma').setMax(2.)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(0.)
ws.var('phis').setConstant(kFALSE)
#With phis unconstrained we now also have sensitivity to deltaperp!!!!
ws.var('deltaperp').setMin(-2*pi)
ws.var('deltaperp').setMax(2*pi)
ws.var('deltaperp').setConstant(kFALSE)

phis = ws.var('phis')
deltaGamma = ws.var('t_sig_dG')

#From the old function
BaseNLL = RooRealVar('BaseNLL','BaseNLL',0.)
GridPointNLL = RooRealVar("GridPointNLL","GridPointNLL",0.)

stepphis = RooRealVar("stepphis","stepphis",0.)
valuephis = RooRealVar("valuephis","valuephis",0.)

stepdeltaGamma = RooRealVar("stepdeltaGamma","stepdeltaGamma",0.)
valuedeltaGamma = RooRealVar("valuedeltaGamma","valuedeltaGamma",.0)

dlogLLArgSet = RooArgSet(BaseNLL,GridPointNLL,stepphis,valuephis,stepdeltaGamma,valuedeltaGamma)
dlogLLDataSet = RooDataSet('dlogLLDataSet','dlogLLDataSet',dlogLLArgSet)

phis.setMin(0)
phis.setMax(2*pi)

phis_min = 0
phis_max = 2*pi
phis_int = phis_max-phis_min

#get range for deltaGamma
deltaGamma_min = -0.7
deltaGamma_max = 0.7
deltaGamma_int = deltaGamma_max - deltaGamma_min

nll = pdf.createNLL(data,RooFit.NumCPU(8),RooFit.Extended(true))
#profile = nll.createProfile(RooArgSet(phis,deltaGamma))

#Check some inital stuff
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
BaseNLL.setVal(nll.getVal())
print 'BaseNLL =',BaseNLL.getVal()

for i in range(_npoints):  
    phis.setVal((phis_min+(phis_int/(2.*_npoints)))+i*(phis_int/_npoints))
    deltaGamma.setVal((deltaGamma_min+(deltaGamma_int/(2.*_npoints)))+_step*(deltaGamma_int/_npoints))
    print '***************************************************************************'
    print 'At gridpoint i = %i from %i and j = %i from %i'%(i+1,_npoints,_step+1,_npoints)
    print '%s at current gridpoint ='%(phis.GetName()), phis.getVal()
    print '%s at current gridpoint ='%(deltaGamma.GetName()), deltaGamma.getVal()
    print '***************************************************************************'
    stepphis.setVal(i+1)
    valuephis.setVal(phis.getVal())
    stepdeltaGamma.setVal(_step+1)
    valuedeltaGamma.setVal(deltaGamma.getVal())
    
    phis.setConstant(kTRUE)
    deltaGamma.setConstant(kTRUE)
    result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))
    GridPointNLL.setVal(nll.getVal())
    print 'BaseNLL =', BaseNLL.getVal()
    print 'GridPointNLL =',GridPointNLL.getVal()

    #dlogLL.setVal(profile.getVal())
    
    dlogLLDataSet.add(dlogLLArgSet)
            
tfile = TFile('profile.root','RECREATE')
dlogLLDataSet.Write() 
tfile.Close()

