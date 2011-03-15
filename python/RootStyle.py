from ROOT import gStyle, TStyle, gROOT
from array import array

def getStyle():
    P2VVStyle = TStyle("P2VVStyle","P2VV's default style")
#    P2VVStyle.SetFrameBorderMode(1)
#    P2VVStyle.SetCanvasBorderMode(1)
#    P2VVStyle.SetPadBorderMode(1)
    P2VVStyle.SetPadColor(0)
    P2VVStyle.SetCanvasColor(0)
    #P2VVStyle.SetStatColor(0)
    P2VVStyle.SetFillColor(0)
    ### set the paper & margin sizes
    #P2VVStyle.SetPaperSize(20,26)
    #P2VVStyle.SetPadTopMargin(0.05)
    #P2VVStyle.SetPadRightMargin(0.05)
    P2VVStyle.SetPadBottomMargin(0.16)
    P2VVStyle.SetPadLeftMargin(0.12)

    ### use large Times-Roman fonts
    P2VVStyle.SetTextFont(132)
    P2VVStyle.SetTextSize(0.04)

    P2VVStyle.SetLabelFont(132,"x")
    P2VVStyle.SetLabelFont(132,"y")
    P2VVStyle.SetLabelFont(132,"z")

    P2VVStyle.SetLabelSize(0.04,"x")
    P2VVStyle.SetLabelSize(0.04,"y")
    P2VVStyle.SetLabelSize(0.04,"z")

    P2VVStyle.SetTitleSize(0.04,"x")
    P2VVStyle.SetTitleSize(0.03,"y")
    P2VVStyle.SetTitleSize(0.06,"z")

    P2VVStyle.SetTitleXOffset(1.5)
    P2VVStyle.SetTitleYOffset(1.7)

    P2VVStyle.SetTitleFont( 132, "xyz" )
    P2VVStyle.SetTitleSize(0.04 )

    P2VVStyle.SetLabelOffset(0.025,"x")
    P2VVStyle.SetLabelOffset(0.005,"y")
    ### use bold lines and markers
    P2VVStyle.SetMarkerStyle(1)
    #P2VVStyle.SetHistLineWidth(1.85)
    P2VVStyle.SetLineStyleString(2,"[12 12]") 
    P2VVStyle.SetOptTitle(1)
    P2VVStyle.SetOptStat(1)

    #P2VVStyle.SetOptFit(1)

    colors = [0,11,5,7,3,6,2,4,12,1]
    colours = array('i',colors)
    P2VVStyle.SetPalette(10,colours)
    P2VVStyle.SetOptFit(0111)
    P2VVStyle.SetOptStat(0)
    P2VVStyle.SetStatX(0.9)
    P2VVStyle.SetStatY(0.9)
    P2VVStyle.SetStatW(0.2)
    P2VVStyle.SetStatH(0.1)
    P2VVStyle.SetStatFont(12)
    
    return P2VVStyle

