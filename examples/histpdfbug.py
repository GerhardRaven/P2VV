from math import  pi
from RooFitWrappers import *

obj  = RooObject( workspace = 'workspace')
cpsi   = RealVar( Name = 'trcospsi',   Title = 'cos(#psi)',        MinMax=(-1,1), nBins = 24 )
ctheta = RealVar( Name = 'trcostheta', Title = 'cos(#theta_{tr})', MinMax=(-1,1), nBins = 24 )
m      = RealVar( Name = 'mass',  Title = 'M(J/#psi#phi)', Unit = 'MeV',  MinMax = (5259, 5451), nBins =  48
                       ,  Ranges =  { 'leftsideband'  : ( None, 5330 )
                                    , 'signal'        : ( 5330, 5410 )
                                    , 'rightsideband' : ( 5410, None ) 
                                    } )
from P2VVGeneralUtils import readData
data = readData( '/data/bfys/dveijk/DataJpsiPhi/2012/Bs2JpsiPhi_ntupleB_for_fitting_20111220.root'
               , dataSetName = 'DecayTree'
               , NTuple = True
               , observables = [ m,cpsi,ctheta ]
               )

side       = data.reduce( CutRange = 'leftsideband'  )
side.append( data.reduce( CutRange = 'rightsideband' ) )

# data = side
# why doesn't:
# side = data.reduce( CutRange = 'leftsideband,rightsideband' ) 
# work? Note that specifying two CutRanges is the 'and' and we want the 'or'....
# All that is needed is that AbsArg.inRange('left,right') returns the 'or' of the two range...

pdf2 = HistPdf( Name = 'bkg_angles'
              , Observables = [ cpsi,ctheta ]
              , Data = data
              )

from ROOT import TCanvas, RooFit
_c = TCanvas('2D')
_c2 = TCanvas('2_x_1D')

for (_cc,a) in zip(_c.pads(2), [cpsi,ctheta] ) :
    f = a.frame()
    data.plotOn(f)
    pdf2.dataHist().plotOn(f, RooFit.MarkerStyle(4), RooFit.MarkerSize(2) )
    pdf2.plotOn(f)
    f.Draw(pad = _cc)

pdf1 = dict()
for a in [cpsi,ctheta] :
    pdf1[a.GetName()] = HistPdf( Name = 'bkg_angles_' + a.GetName()
                               , Observables = (a,)
                               , Data = data
                               )
for (_cc,a) in zip(_c2.pads(2), [cpsi,ctheta] ) :
    f = a.frame()
    data.plotOn(f)
    pdf1[a.GetName()].dataHist().plotOn(f, RooFit.MarkerStyle(4), RooFit.MarkerSize(2) )
    pdf1[a.GetName()].plotOn(f)
    f.Draw(pad = _cc)

