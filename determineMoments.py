from ROOT import *
gSystem.Load('libp2vv')
gStyle.SetOptStat(0);
from ModelBuilders import buildMomentPDF,abasis,declareObservables,buildMassPDFs
from itertools import count,product
from math import pi


def doit(name, angles, tree, irange, lrange, mrange ) :
    # TODO: veto signal mass window... or use sweight to veto signal
    # TODO: use sweights to split J/psi from mumu combinatoric
    # TODO: veto insignifcant moments iff i,l,abs(m)>2 (i.e. those not in signal PDF!)
    data = RooDataSet('data','data',tree,angles)
    ab = abasis(w,angles)
    moments = []
    for (i,l,m) in product(irange,lrange,mrange) :
          if abs(m)>l : continue
          #  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
          moments.append( Moment( ab.build(name+'_mom',i,0,l,m,1.), float(2*i+1)/2 ) )
    pdf = buildMomentPDF( w, name, data, moments )

    c = TCanvas()
    c.Divide(3,2)
    for (i,v) in enumerate( angles ) :
        c.cd(1+i)
        frame = v.frame()
        data.plotOn(frame)
        pdf.plotOn(frame)
        frame.Draw()

        c.cd(4+i)
        others = RooArgList( angles )
        others.remove( v )
        hist = pdf.createHistogram( others.names() )
        pdf.fillHistogram( hist,others,1., RooArgSet(v))
        hist.Draw('COLZ')
        # hist.Draw('sinusoidal')
    return c

fname = 'Bs2JpsiPhiTuple.root'
dataName = 'dataset'

f = TFile(fname)
tree = f.Get(dataName)

w = RooWorkspace("w")
declareObservables(w)

## TODO: invoke mass fit, and compute SWeights
#buildMassPDFs(w)

# transversity needs many more trphi moments
c2 = doit("bkg_trangles_pdf", w.set('transversityangles'), tree, range(3), range(21), range(-8,9) )
# helicity has phi almost flat...
c1 = doit("bkg_helangles_pdf",w.set('helicityangles'), tree, range(3), range(21), range(-2,3) )

