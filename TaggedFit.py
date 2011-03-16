from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *

#import rootStyle
#from ROOT import (gROOT,gStyle,TStyle)
#myStyle = rootStyle.plainstyle()
#gROOT.SetStyle(myStyle.GetName())
#gROOT.ForceStyle()
#gStyle.UseCurrentStyle()

taggingsyst = True
blinded = False
reblinded = True

wsfile = TFile('TaggedWS.root')
ws = wsfile.Get('ws')

if taggingsyst:
    pdfbeforeblinding = ws.pdf('pdf_ext_angcorrpdf_inc_tag_syst')
else:
    pdfbeforeblinding = ws.pdf('pdf_ext_angcorrpdf')

data = ws.data('data')

#For tagged fit with phis set to zero:
#reblinded = False
#blinded = False
#ws.var('phis').setVal(0.0)
#ws.var('phis').setConstant(kTRUE)

################
### Blinding ###
################

if blinded:
    #Building blinded parameters
    ws.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
    ws.factory("expr::t_sig_tau_blind('1/@0',{#Gamma_unblind})")

    ws.factory("RooUnblindUniform::t_sig_dG_blind('BsHobbes',0.2,t_sig_dG)")

    ws.factory("RooUnblindUniform::phis_blind('BsGoofy',3,phis)")

    #calling RooCustomizer
    customizer = RooCustomizer(pdfbeforeblinding,'blinded')

    tau_unblind = ws.function('t_sig_tau')
    tau_blind = ws.function('t_sig_tau_blind')

    dG_unblind = ws.var('t_sig_dG')
    dG_blind = ws.function('t_sig_dG_blind')

    #Replacing
    customizer.replaceArg( tau_unblind, tau_blind )
    customizer.replaceArg( dG_unblind, dG_blind )

    blindedpdf = customizer.build()
    ws.put(blindedpdf)


if reblinded:
    ws.factory("RooUnblindUniform::phis_blind('BsOlivier',3,phis)")

    customizer = RooCustomizer(pdfbeforeblinding,'reblinded')

    phis_unblind = ws.var('phis')
    phis_blind = ws.function('phis_blind')

    customizer.replaceArg( phis_unblind, phis_blind )

    reblindedpdf = customizer.build()
    ws.put(reblindedpdf)

if reblinded:
    pdf = reblindedpdf
elif blinded:
    pdf = blindedpdf
else:
    pdf = pdfbeforeblinding

#Fix deltams
#ws.var('t_sig_dm').setVal(17.8)
#ws.var('t_sig_dm').setConstant(kTRUE)
#result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

ws.var('#Gamma').setVal(0.68)
ws.var('#Gamma').setMax(2.)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(pi)
ws.var('phis').setMin(0.)
ws.var('phis').setMax(2*pi)

#Constrain deltams
ws.factory("Gaussian::dmsconstraint(t_sig_dm,t_sig_dm_mean[17.77],t_sig_dm_sigma[0.12])")
result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true),RooFit.ExternalConstraints(RooArgSet(ws.pdf('p0'),ws.pdf('p1'),ws.pdf('dmsconstraint'))))

paramlist = pdf.getParameters(data)
writeFitParamsLatex(paramlist)

dict = writeCorrMatrixLatex(result)

################
### Profiles ###
################

#setting back values
ws.var('#Gamma').setVal(0.68)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(0.0)

phis = ws.var('phis')
deltaGamma = ws.var('t_sig_dG')

#MakeProfile('ProfiledGamma_phis_tagged',data,pdf,15,phis,0,2*pi,deltaGamma,-0.7,0.7)
