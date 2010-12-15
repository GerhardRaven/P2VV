from ROOT import (gStyle, TStyle, gROOT)
from array import array

def plainstyle():
    daanStyle = TStyle("Daan","Daan's Style")
#    daanStyle.SetFrameBorderMode(1)
#    daanStyle.SetCanvasBorderMode(1)
#    daanStyle.SetPadBorderMode(1)
    daanStyle.SetPadColor(0)
    daanStyle.SetCanvasColor(0)
    #daanStyle.SetStatColor(0)
    daanStyle.SetFillColor(0)
    ### set the paper & margin sizes
    #daanStyle.SetPaperSize(20,26)
    #daanStyle.SetPadTopMargin(0.05)
    #daanStyle.SetPadRightMargin(0.05)
    daanStyle.SetPadBottomMargin(0.16)
    daanStyle.SetPadLeftMargin(0.12)

    ### use large Times-Roman fonts
    daanStyle.SetTextFont(132)
    daanStyle.SetTextSize(0.04)

    daanStyle.SetLabelFont(132,"x")
    daanStyle.SetLabelFont(132,"y")
    daanStyle.SetLabelFont(132,"z")

    daanStyle.SetLabelSize(0.04,"x")
    daanStyle.SetLabelSize(0.04,"y")
    daanStyle.SetLabelSize(0.04,"z")

    daanStyle.SetTitleSize(0.04,"x")
    daanStyle.SetTitleSize(0.03,"y")
    daanStyle.SetTitleSize(0.06,"z")

    daanStyle.SetTitleXOffset(1.5)
    daanStyle.SetTitleYOffset(1.7)

    daanStyle.SetTitleFont( 132, "xyz" )
    daanStyle.SetTitleSize(0.04 )

    daanStyle.SetLabelOffset(0.025,"x")
    daanStyle.SetLabelOffset(0.005,"y")
    ### use bold lines and markers
    daanStyle.SetMarkerStyle(1)
    #daanStyle.SetHistLineWidth(1.85)
    daanStyle.SetLineStyleString(2,"[12 12]") 
    daanStyle.SetOptTitle(1)
    daanStyle.SetOptStat(1)

    #daanStyle.SetOptFit(1)

    colors = [0,11,5,7,3,6,2,4,12,1]
    colours = array('i',colors)
    daanStyle.SetPalette(10,colours)
    daanStyle.SetOptFit(0111)
    daanStyle.SetOptStat(0)
    daanStyle.SetStatX(0.9)
    daanStyle.SetStatY(0.9)
    daanStyle.SetStatW(0.2)
    daanStyle.SetStatH(0.1)
    daanStyle.SetStatFont(12)
    
    return daanStyle

#In your file: 
#import RootStyle
#from ROOT import (gROOT,gStyle,TStyle)
#myStyle = RootStyle.plainstyle()
#gROOT.SetStyle(myStyle.GetName())
#gROOT.ForceStyle()
#gStyle.UseCurrentStyle()


#How to get all the commands in plainstyle()?
#save a canvas as a .root file, go to root, 'new TBrowser', find the file, open the editor and edit your plot. Then save as .c file (rootmacro). Open this .c file in emacs. Then you can see all the commands needed to create that particular Canvas.

