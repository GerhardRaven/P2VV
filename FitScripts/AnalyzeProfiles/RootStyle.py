from ROOT import (gStyle, TStyle, gROOT)
from array import array

#colors = [0,11,5,7,3,6,2,4,12,1]
#colors = [8,11,5,7,3,6,2,4,12,9]
#colours = array('i',colors)

def MyStyle():
    MyStyle = TStyle("MyStyle","MyStyle")

    #MyStyle.SetPalette(10,colours)
    #MyStyle.SetPalette(1)
    #MyStyle.SetOptFitW(0.3)

    MyStyle.SetFrameBorderMode(0)

    #MyStyle.SetOptStat(0)
    #MyStyle.SetOptFit(1)
    #MyStyle.SetOptFit(1111)
    #MyStyle.SetStatX(0.9)
    #MyStyle.SetStatY(0.9)
    #MyStyle.SetStatW(0.2)
    #MyStyle.SetStatH(0.1)

    fontnr = 42
    fontsize = 0.05

    MyStyle.SetStatFont(fontnr)

    MyStyle.SetTitleFont(fontnr,"XYZ") 
    MyStyle.SetLabelFont(fontnr,"XYZ")
    MyStyle.SetTitleSize(fontsize,"XYZ")
    MyStyle.SetLabelSize(fontsize,"XYZ")
    #MyStyle.SetTitleOffset(1.5,"XYZ")
    MyStyle.SetLabelOffset(0.015,"X")
    #MyStyle.SetCanvasBorderSize(2)
    MyStyle.SetCanvasColor(0)
    MyStyle.SetStatColor(0)
    MyStyle.SetTitleFillColor(0)
    MyStyle.SetFrameFillColor(0)
    MyStyle.SetPadColor(0)
    MyStyle.SetPadLeftMargin(0.13)
    MyStyle.SetPadBottomMargin(0.13)
    #MyStyle.SetLegendBorderSize(1)
    #MyStyle.SetTitleBorderSize(1)
    #MyStyle.SetStatBorderSize(1)
    
    #MyStyle.SetLineWidth(1)

    return MyStyle
