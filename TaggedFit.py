from ROOT import *
gSystem.Load("libp2vv")
from math import sqrt,pi
from array import array

from RooFitDecorators import *
import rootStyle
from ModelBuilders import _buildAngularFunction
#from ROOT import (gROOT,gStyle,TStyle)
myStyle = rootStyle.plainstyle()
gROOT.SetStyle(myStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

wsfile = TFile('TaggedWS.root')
ws = wsfile.Get('ws')

angcorrpdf = ws.pdf('pdf_ext_angcorrpdf')
data = ws.data('data')

#For tagged fit with phis set to zero:
#reblinded = False
#blinded = False
#ws.var('phis').setVal(0.0)
#ws.var('phis').setConstant(kTRUE)

################
### Blinding ###
################

blinded = False
reblinded = True

if blinded:
    #Building blinded parameters
    ws.factory("RooUnblindUniform::#Gamma_unblind('BsCalvin',0.4,#Gamma)")
    ws.factory("expr::t_sig_tau_blind('1/@0',{#Gamma_unblind})")

    ws.factory("RooUnblindUniform::t_sig_dG_blind('BsHobbes',0.2,t_sig_dG)")

    ws.factory("RooUnblindUniform::phis_blind('BsGoofy',3,phis)")

    #calling RooCustomizer
    customizer = RooCustomizer(angcorrpdf,'blinded')

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

    customizer = RooCustomizer(angcorrpdf,'reblinded')

    phis_unblind = ws.var('phis')
    phis_blind = ws.function('phis_blind')

    customizer.replaceArg( phis_unblind, phis_blind )

    reblindedpdf = customizer.build()
    ws.put(reblindedpdf)

if reblinded:
    pdf = ws.pdf('pdf_ext_angcorrpdf_reblinded')
elif blinded:
    pdf = ws.pdf('pdf_ext_angcorrpdf_blinded')
else:
    pdf = ws.pdf('pdf_ext_angcorrpdf')



result = pdf.fitTo(data,RooFit.NumCPU(8),RooFit.Extended(true),RooFit.Minos(false),RooFit.Save(true))

paramlist = pdf.getParameters(data)
writeFitParamsLatex(paramlist)

dict = writeCorrMatrixLatex(result)

assert False

#Make PDF conditional on WA deltams
ws.factory("Gaussian::dmGauss(t_sig_dm,t_sig_dm_mean[17.77],t_sig_dm_sigma[0.12])")
ws.factory("PROD::condpdf(pdf_ext_fourier_eff_full|t_sig_dm,dmGauss)")
condpdf = ws.pdf('condpdf')

frame = ws.var('t_sig_dm').frame()
ws.pdf('dmGauss').plotOn(frame)
frame.Draw()


################
### Profiles ###
################

#setting back values
ws.var('#Gamma').setVal(0.68)
ws.var('t_sig_dG').setVal(0.060)
ws.var('phis').setVal(0.0)

phis = ws.var('phis')
deltaGamma = ws.var('t_sig_dG')

MakeProfile('ProfiledGamma_phis_tagged',data,angcorrpdf,15,phis,-pi,pi,deltaGamma,-1,1)
