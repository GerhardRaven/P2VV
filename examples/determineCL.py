from itertools import count
from math import pi
from ModelBuilders import buildJpsiphi,buildJpsikstar,declareObservables,definePolarAngularAmplitudes
from P2VV.RooFitDecorators import *
from ROOT import *

if False :
    ws = RooWorkspace('ws')
    declareObservables(ws)
    definePolarAngularAmplitudes(ws)

    ws.factory("RooGaussModel::tres_sig(t,mu[0],sigma[0.05])")
    ws.factory("{wmistag[0.0]}")

    useTransversityAngles = False
    pdf = buildJpsiphi(ws,'jpsiphipdf',useTransversityAngles)  

    obs = ws.set('transversityangles' if useTransversityAngles else 'helicityangles')
    obs.add( ws.argSet('t,tagdecision') )
    data = pdf.generate( obs, 10000)
    ws.put(data)
    file = TFile("p2vv_datapdf.root","RECREATE")
    ws.Write("w")
    file.Close()
else :
    file = TFile("p2vv_datapdf.root")
    ws = file.Get("w")
    pdf = ws['jpsiphipdf']
    data = ws['jpsiphipdfData']


phis = ws['phis']
dG   = ws['dG']
# let phi_s vary!
phis.setMin(-pi)
phis.setMax(pi)
phis.setConstant(False)

if False :
    iplc = RooStats.ProfileLikelihoodCalculator(data,pdf, RooArgSet(phis) )
    iplc.SetConfidenceLevel(0.95) # 95% interval
    plot = RooStats.LikelihoodIntervalPlot(  iplc.GetInterval() )
    plot.Draw();
            
if True :
    opt = ''
    for (lvl,col) in zip( [ 0.05, 0.10, 0.32 ],[ kYellow, kOrange, kRed ]) :
        plc = RooStats.ProfileLikelihoodCalculator(data, pdf, RooArgSet(phis,dG ) )
        plc.SetTestSize(lvl) 
        contourPlot = RooStats.LikelihoodIntervalPlot( plc.GetInterval() )
        contourPlot.SetNPoints(7)
        contourPlot.SetLineColor(kBlack)
        contourPlot.SetContourColor(col)
        contourPlot.Draw(','.join(['nominuit',opt]))
        if not opt : opt = 'same'
