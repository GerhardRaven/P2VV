from ROOT import *
gSystem.Load("libp2vv")
import math
gStyle.SetPalette(1)
from math import pi

def abasis(w,label,i,j,k,l,c):
    ctheta = w.var('ctheta')
    cpsi = w.var('cpsi')
    phi = w.var('phi')
    string = "%s_%d_%d"%(label,i,j)
    if l<0:
        name = "%s_%d_m%d"%(string,k,l)
    else:
        name = "%s_%d_%d"%(string,k,l)
    b = w.function(name)
    basisfunc = RooP2VVAngleBasis(name,name,cpsi,ctheta,phi,i,j,k,l,c)
    if not b:
        im = getattr(w,'import')(basisfunc)
    return basisfunc

def moment(data,basis,pdf_marginal,allObs,factor):
    m0 = 0
    m1 = 0
    m2 = 0
    
    for i in range(0,data.numEntries()):
        allObs.assignValueOnly(data.get(i))
        #print 'basisval = ', basis.getVal()
        #print 'pdfval = ', pdf_marginal.getVal(allObs)
        #div = basis.getVal()/pdf_marginal.getVal(allObs)
        div = basis.getVal()
        #print 'div = ', div
        m0+=1
        m1+=div
        m2+=div*div

    mu = factor*m1/m0
    sig2 = m2/m0 - mu*mu
    muerr = sqrt(sig2/(m0-1))
    musigni = mu/sqrt(sig2/m0)
    print  "moment(", basis.GetName(),"_", pdf.GetName(), ") = ",mu," +- ", muerr ," significance: ", musigni
    print 'm0 =',m0
    print 'm1 =',m1
    print 'm2 =',m2
    print 'mu =',mu
    print 'merr =',muerr
    print 'musigni =',musigni

    ret = {'mu':mu,'muerr':muerr,'musigni':musigni}
    return ret


tfile = TFile("p2vv_3.root")
#tfile = TFile("p2vv_onlytag1.root")
w = tfile.Get('w')

ctheta = w.var('ctheta')
cpsi = w.var('cpsi')
phi = w.var('phi')


w.factory("RooUniform::uni({ctheta,cpsi,phi})")

sigpdf = w.pdf('pdf')

pdf = w.pdf('uni')#uniform for background description!
pdf.Print()

data = w.data('pdfData') 

#create moments -- for this, we need the PDF used to generate the data
allObs = pdf.getObservables( data.get() )
marginalObs = pdf.getObservables( data.get() )
marginalObs.remove( w.argSet("cpsi,ctheta,phi") )

#marginalize pdf over 'the rest' so we get the normalization of the moments right...
pdf_marginal = pdf.createProjection(marginalObs)

mu_vars = []
coeffList = RooArgList('coeffList')

muerr_vars = []
coefferrList = RooArgList('coefferrList')

musigni_vars = []
coeffsigniList = RooArgList('coeffsigniList')

basislist = []

for i in range(0,3):
    for l in range(0,3):
        for m in range(-l,l+1):
            basis = abasis(w,'mom',i,0,l,m,float(2*i+1)/2. )
            factor = 2./float(2*i+1)
            
            basislist.append(basis)
            
            momentdict = moment(data,basis,pdf_marginal,allObs,factor)
            
            mu_vars.append(RooRealVar('C_%s_0_%s_%s'%(i,l,m),'C_%s_0_%s_%s'%(i,l,m),momentdict['mu']))
            coeffList.add(mu_vars[-1])
            
            muerr_vars.append(RooRealVar('Cerr_%s_0_%s_%s'%(i,l,m),'Cerr_%s_0_%s_%s'%(i,l,m),momentdict['muerr']))
            coefferrList.add(muerr_vars[-1])
            
            musigni_vars.append(RooRealVar('Csigni_%s_0_%s_%s'%(i,l,m),'Csigni_%s_0_%s_%s'%(i,l,m),momentdict['musigni']))
            coeffsigniList.add(musigni_vars[-1])

basisList = RooArgList('basisList')
for i in basislist:
    basisList.add(i)

basisList.Print()
coeffList.Print()

#basis0000  = abasis(w,'mom',0,0,0,0,1.)
#basislist.append(basis0000)
## basis0020  = abasis(w,'mom',0,0,2,0,1.)
## basislist.append(basis0020)
## basis0022  = abasis(w,'mom',0,0,2,2,1.)
## basislist.append(basis0022)
## basis2000  = abasis(w,'mom',2,0,0,0,1.)
## basislist.append(basis2000)
## basis2020  = abasis(w,'mom',2,0,2,0,1.)
## basislist.append(basis2020)
## basis2022  = abasis(w,'mom',2,0,2,2,1.)
## basislist.append(basis2022)
## basis2200  = abasis(w,'mom',2,2,0,0,1.)
## basislist.append(basis2200)
## basis2220  = abasis(w,'mom',2,2,2,0,1.)
## basislist.append(basis2220)
## basis2222  = abasis(w,'mom',2,2,2,2,1.)
## basislist.append(basis2222)
## basis222m1 = abasis(w,'mom',2,2,2,-1,1.)
## basislist.append(basis222m1)
## basis2121  = abasis(w,'mom',2,1,2,1,1.)
## basislist.append(basis2121)
## basis212m2 = abasis(w,'mom',2,1,2,-2,1.)
## basislist.append(basis212m2)

#This will only be interesting for jpiskstar where the time and angle parts separate, then you can indeed just put the coefficients from the table....

#ShouldBecoeffList = RooArgList('ShouldBecoeffList')
#shouldbevars = []

#for b in basislist:
#    momentdict = moment(data,b,pdf_marginal,allObs)
#    mu_vars.append(RooRealVar('C_%s'%(b.getTitle()),'C_%s'%(b.getTitle()),momentdict['mu']))

#    coeffList.add(mu_vars[-1])
#    muerr_vars.append(RooRealVar('C_err%s'%(b.getTitle()),'C_err%s'%(b.getTitle()),momentdict['muerr']))
#    coefferrList.add(muerr_vars[-1])
    
#    musigni_vars.append(RooRealVar('C_signi%s'%(b.getTitle()),'C_signi%s'%(b.getTitle()),momentdict['musigni']))
#    coeffsigniList.add(musigni_vars[-1])

## for coeff in range(0,coeffList.getSize()) :
##     mean = (coeffList[coeff]).getVal()
##     error = (coefferrList[coeff]).getVal()/math.sqrt(data.numEntries())
##     if (coeff==0 or abs(mean)<3.*error):
##         (coeffList[coeff]).setConstant(kTRUE)
##         print "setting this constant"
##     else :
##         min =mean-3*error
##         max =mean+3*error
##         (coeffList[coeff]).setError(error)
##         (coeffList[coeff]).setMin(min)
##         (coeffList[coeff]).setMax(max)
##         (coeffList[coeff]).setConstant(kFALSE)
##         print "floating this parameter; with min = ", min , " and max = ", max

signifcoeffList = RooArgList('signifcoeffList')
signifbasisList = RooArgList('signifbasisList')
for coeff in range(0,coeffList.getSize()):
    signif = coeffsigniList[coeff].getVal()
    if abs(signif)>500:
        signifcoeffList.add(coeffList[coeff])
        signifbasisList.add(basisList[coeff])

#make the roorealsumpdf:
#bkgPDF = RooRealSumPdf("bkgPDF","bkgPDF",signifbasisList,signifcoeffList)
bkgPDF = RooRealSumPdf("bkgPDF","bkgPDF",basisList,coeffList)

#ShouldbebkgPDF = RooRealSumPdf("bkgPDF","bkgPDF",basisList,ShouldBecoeffList)

#bkgResult = bkgPDF.fitTo(data,RooFit.Save(),RooFit.NumCPU(4));

#w.Print()

cthetaframe = w.var('ctheta').frame()
cpsiframe = w.var('cpsi').frame()
phiframe = w.var('phi').frame()

#cutdata = data.reduce('phi>%d'%pi)
#cutdata = data.reduce('cpsi>0')
                      
data.plotOn(cthetaframe)
bkgPDF.plotOn(cthetaframe)
#sigpdf.plotOn(cthetaframe,RooFit.LineColor(RooFit.kRed))
#ShouldbebkgPDF.plotOn(cthetaframe,RooFit.LineColor(RooFit.kRed))

data.plotOn(cpsiframe)
bkgPDF.plotOn(cpsiframe)
#sigpdf.plotOn(cpsiframe,RooFit.LineColor(RooFit.kRed))
#ShouldbebkgPDF.plotOn(cpsiframe,RooFit.LineColor(RooFit.kRed))

data.plotOn(phiframe)
bkgPDF.plotOn(phiframe)
#sigpdf.plotOn(phiframe,RooFit.LineColor(RooFit.kRed))
#ShouldbebkgPDF.plotOn(phiframe,RooFit.LineColor(RooFit.kRed))


c =TCanvas('c','c')
c.Divide(3)

c.cd(1)
cpsiframe.Draw()

c.cd(2)
cthetaframe.Draw()

c.cd(3)
phiframe.Draw()


#Nice plots

## #try:
## cthetaframe = w.var('ctheta').frame()
## cpsiframe = w.var('cpsi').frame()
## phiframe = w.var('phi').frame()

## leg00 = RooLegendre('leg00','leg00',w.var('cpsi'),0,0)
## leg20 = RooLegendre('leg20','leg20',w.var('cpsi'),2,0)
## leg21 = RooLegendre('leg21','leg21',w.var('cpsi'),2,1)
## leg22 = RooLegendre('leg22','leg22',w.var('cpsi'),2,2)

## cleg = TCanvas('cleg','cleg')
## leg00.plotOn(cpsiframe)
## leg20.plotOn(cpsiframe)
## leg21.plotOn(cpsiframe)
## leg22.plotOn(cpsiframe)
## cpsiframe.Draw()

## csph = TCanvas('csph','csph',1300,400)
## csph.Divide(6)
## ctheta = w.var('ctheta')
## phi = w.var('phi')

## sph00 = RooSpHarmonic('sph00','sph00',w.var('ctheta'),w.var('phi'),0,0)
## sph20 = RooSpHarmonic('sph20','sph20',w.var('ctheta'),w.var('phi'),2,0)
## sph22 = RooSpHarmonic('sph22','sph22',w.var('ctheta'),w.var('phi'),2,2)
## sph21 = RooSpHarmonic('sph21','sph21',w.var('ctheta'),w.var('phi'),2,1)
## sph2m1 = RooSpHarmonic('sph2m1','sph2m1',w.var('ctheta'),w.var('phi'),2,-1)
## sph2m2 = RooSpHarmonic('sph2m2','sph2m2',w.var('ctheta'),w.var('phi'),2,-2)

## csph.cd(1)
## sph00_2D = sph00.createHistogram("ctheta,phi",50,50)
## sph00_2D.Draw('surf2')

## csph.cd(2)
## sph20_2D = sph20.createHistogram("ctheta,phi",50,50)
## sph20_2D.Draw('surf2')

## csph.cd(3)
## sph22_2D = sph22.createHistogram("ctheta,phi",50,50)
## sph22_2D.Draw('surf2')

## csph.cd(4)
## sph21_2D = sph21.createHistogram("ctheta,phi",50,50)
## sph21_2D.Draw('surf2')

## csph.cd(5)
## sph2m1_2D = sph2m1.createHistogram("ctheta,phi",50,50)
## sph2m1_2D.Draw('surf2')

## csph.cd(6)
## sph2m2_2D = sph2m2.createHistogram("ctheta,phi",50,50)
## sph2m2_2D.Draw('surf2')
