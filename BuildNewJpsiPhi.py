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

#######################################################################################################################################

pdf = buildJpsiphi(ws,'jpsiphipdf')  ## for now we rely quite a bit on a naming convention -- 
                                     ## in future we should pass more information into the builder
                                     ## maybe a dictionary of what's what...


### let's make some nice plots to show what this PDF looks like...

obsNames =[ 'trcospsi','trcostheta','trphi','t','tagdecision' ]
obs = ws.argSet(','.join(obsNames))

canvas = TCanvas('canvas','canvas')
canvas.Divide(5,4)

p = [ 'rz','rpar','rperp' ]
for i in  range(len(p)+1) :
    if i>0 : 
        for j in range(len(p)) : 
            x = ws.var(p[j])
            x.setVal( 1 if i==j+1 else 0 )
    for l in p + ['deltaz','deltapar','deltaperp' ] : ws.var(l).Print()

    data = pdf.generate( obs, 100000)
    if data : pdf.fitTo(data, RooFit.NumCPU(2))

    for (j,k) in zip(['trcospsi','trcostheta','trphi','t'],count(5*i+1)) :
       canvas.cd(k)
       f = ws[j].frame() 
       if data : 
            data.plotOn(f)
            pdf.plotOn(f)
       else :
            proj = RooArgSet( obs )
            proj.remove( ws[j] )
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
        proj.remove( RooArgSet( ws['t'], ws['tagdecision'] ) )
        pdf.plotOn(f,tagAsym, RooFit.Project(proj))

    f.Draw()

canvas.Flush()
