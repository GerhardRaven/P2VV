from itertools import count
from math import pi
from ModelBuilders import buildJpsiphi,buildJpsikstar,declareObservables,definePolarAngularAmplitudes
from RooFitDecorators import *
from ROOT import *

ws = RooWorkspace('ws')
declareObservables(ws)
definePolarAngularAmplitudes(ws)

ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")
ws.factory("{wmistag[0.0]}")

useTransversityAngles = False
pdf = buildJpsiphi(ws,'jpsiphipdf',useTransversityAngles)  ## for now we rely quite a bit on a naming convention -- 
                                     ## in future we should pass more information into the builder
                                     ## maybe a dictionary of what's what...


### let's make some nice plots to show what this PDF looks like...
if useTransversityAngles :
    obs = ws.argSet('trcospsi,trcostheta,trphi,t,tagdecision')
else:
    obs = ws.argSet('helcosthetaK,helcosthetaL,helphi,t,tagdecision')

canvas = TCanvas('canvas','canvas')
canvas.Divide(5,4)

p = ws.argSet('rz,rpar,rperp')
for i in  range(len(p)+1) :
    if i>0 : 
        for (j,x) in enumerate(p) : x.setVal( 1 if i==j+1 else 0 )
    for l in p.nameList() + ['deltaz','deltapar','deltaperp'] : ws.var(l).Print()

    data = None  # pdf.generate( obs, 100000)
    if data : 
        if i==0 : 
            ws.put(data)
            file = TFile("p2vv_9.root","RECREATE")
            ws.Write("w")
            file.Close()
        pdf.fitTo(data, RooFit.NumCPU(2))

    angt = RooArgSet( obs )
    angt.remove( ws['tagdecision'] )
    for (j,k) in zip( angt,count(5*i+1)) :
       canvas.cd(k)
       f = j.frame() 
       if data : 
            data.plotOn(f)
            pdf.plotOn(f)
       else :
            proj = RooArgSet( obs )
            proj.remove( j )
            pdf.plotOn(f,RooFit.Project(proj))
       f.Draw()

    tagAsym = RooFit.Asymmetry(ws["tagdecision"])
    canvas.cd(5*i+5)
    f = ws['t'].frame(RooFit.Range(-.5,3.5))
    if data :
        data.plotOn(f,tagAsym)
        pdf.plotOn(f,tagAsym)
    else :
        proj = RooArgSet( obs )
        proj.remove( ws.argSet('t,tagdecision' ) )
        pdf.plotOn(f,tagAsym, RooFit.Project(proj))
    f.Draw()

canvas.Flush()
