"""sets ROOT style

Sets the ROOT style for plots and other graphics
"""

print "P2VV - INFO: setting ROOT style"
from ROOT import gStyle, gROOT

# set ROOT plot style
plotStyle = gROOT.GetStyle('Plain')

plotStyle.SetTitleBorderSize(0)

plotStyle.SetPadTopMargin(0.10)
plotStyle.SetPadBottomMargin(0.15)
plotStyle.SetPadLeftMargin(0.15)
plotStyle.SetPadRightMargin(0.10)

plotStyle.SetTitleSize(0.048, 'XY')
plotStyle.SetLabelSize(0.045, 'XY')

plotStyle.SetTitleOffset(1.2, 'X')
plotStyle.SetTitleOffset(1.6, 'Y')
plotStyle.SetLabelOffset(0.011, 'XY')

gROOT.SetStyle('Plain')
gROOT.ForceStyle()

gStyle.SetPalette(1)
gStyle.UseCurrentStyle()

