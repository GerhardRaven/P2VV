from ROOT import *
from array import array

#filename = 'TaggedProfileLikelihood'
filename = 'UntaggedProfileLikelihood'
name = 'ProfileLikelihood'

tfile =TFile('%s.root'%(filename))
ProfileLikelihood = tfile.Get(name)

x = ProfileLikelihood.GetXaxis()
y = ProfileLikelihood.GetYaxis()

gStyle.SetPalette(1)
gStyle.SetOptStat(0)

#contours = [2.30/2.,4.61/2.,5.99/2.,9.21/2.]
contours = [2.30/2.,4.61/2.,5.99/2.]
contourarray = array('d',contours)

Canvas = TCanvas('Canvas','Canvas')
#Canvas.Divide(2)

#Canvas.cd(1)
#ProfileLikelihood.Draw('surf2')
#ProfileLikelihood.Draw('CONT1Z SAME')
#ProfileLikelihood.Draw('COLZ')

#Canvas.cd(2)
ProfileLikelihood.SetContour(len(contours),contourarray)
ProfileLikelihood.Draw('CONT1 LIST')

ProfileLikelihood.SetLineWidth(2)
ProfileLikelihood.GetXaxis().SetTitle('#Phi_{s}')
ProfileLikelihood.GetYaxis().SetTitle('#Delta#Gamma (ps^{-1})')

## for i in range(1,x.GetNbins()+1):
##     xlinex = x.GetBinCenter(i) 
##     vars()['xline%i'%(i)] = TLine(xlinex,y.GetXmin(),xlinex,y.GetXmax())
##     vars()['xline%i'%(i)].Draw()
## for j in range(1,y.GetNbins()+1):
##     yliney = y.GetBinCenter(j)
##     vars()['yline%i'%(j)] = TLine(x.GetXmin(),yliney,x.GetXmax(),yliney)
##     vars()['yline%i'%(j)].Draw()

#This doesn't work yet 
## for i in range(1,x.GetNbins()+1):
##     for j in range(1,y.GetNbins()+1):
##         if (-0.1 <= ProfileLikelihood.GetBinContent(i,j)<=0.075):
##             pm= TPolyMarker(4)

##             pm.SetPoint(0, x.GetBinCenter(i),y.GetBinCenter(j))
##             pm.SetPoint(1, x.GetBinCenter(i),y.GetBinCenter(j))
##             pm.SetPoint(2, x.GetBinCenter(i),y.GetBinCenter(j))
##             pm.SetPoint(3, x.GetBinCenter(i),y.GetBinCenter(j))            

##             pm.SetMarkerSize(1)
##             pm.SetMarkerColor(1)
##             pm.SetMarkerStyle(2)

##             pm.Draw()
