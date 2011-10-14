from ROOT import *
from ROOT import (gROOT,gStyle,TStyle)

black=1
red=2
green=3
blue=4
yellow=5 
magenta=6
cyan=7
purple=9

def lhcbStyle():
    
    lhcbStyle= TStyle("lhcbStyle","Standard LHCb plots style")

    #use helvetica-bold-r-normal, precision 2 (rotatable)
    lhcbFont = 62
    #line thickness
    lhcbWidth = 3

    #use plain black on white colors
    lhcbStyle.SetFrameBorderMode(0) 
    lhcbStyle.SetCanvasBorderMode(0) 
    lhcbStyle.SetPadBorderMode(0) 
    lhcbStyle.SetPadColor(0) 
    lhcbStyle.SetCanvasColor(0) 
    lhcbStyle.SetStatColor(0) 
    lhcbStyle.SetPalette(1) 

    #set the paper & margin sizes
    lhcbStyle.SetPaperSize(20,26) 
    lhcbStyle.SetPadTopMargin(0.05) 
    lhcbStyle.SetPadRightMargin(0.05)  #increase for colz plots
    lhcbStyle.SetPadBottomMargin(0.16) 
    lhcbStyle.SetPadLeftMargin(0.14) 

    #use large fonts
    lhcbStyle.SetTextFont(lhcbFont) 
    lhcbStyle.SetTextSize(0.08) 
    lhcbStyle.SetLabelFont(lhcbFont,"x") 
    lhcbStyle.SetLabelFont(lhcbFont,"y") 
    lhcbStyle.SetLabelFont(lhcbFont,"z") 
    lhcbStyle.SetLabelSize(0.05,"x") 
    lhcbStyle.SetLabelSize(0.05,"y") 
    lhcbStyle.SetLabelSize(0.05,"z") 
    lhcbStyle.SetTitleFont(lhcbFont) 
    lhcbStyle.SetTitleSize(0.06,"x") 
    lhcbStyle.SetTitleSize(0.06,"y") 
    lhcbStyle.SetTitleSize(0.06,"z") 

    #use bold lines and markers
    lhcbStyle.SetLineWidth(lhcbWidth) 
    lhcbStyle.SetFrameLineWidth(lhcbWidth) 
    lhcbStyle.SetHistLineWidth(lhcbWidth) 
    lhcbStyle.SetFuncWidth(lhcbWidth) 
    lhcbStyle.SetGridWidth(lhcbWidth) 
    lhcbStyle.SetLineStyleString(2,"[12 12]")  #postscript dashes
    lhcbStyle.SetMarkerStyle(20) 
    lhcbStyle.SetMarkerSize(1.5) 

    #label offsets
    lhcbStyle.SetLabelOffset(0.015) 

    #by default, do not display histogram decorations:
    lhcbStyle.SetOptStat(0)   
    lhcbStyle.SetOptStat("emr")   #show only nent -e , mean - m , rms -r
    #full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
    lhcbStyle.SetStatFormat("6.3g")  #specified as c printf options
    lhcbStyle.SetOptTitle(0) 
    lhcbStyle.SetOptFit(0) 
    #lhcbStyle.SetOptFit(1011)  # order is probability, Chi2, errors, parameters

    #look of the statistics box:
    lhcbStyle.SetStatBorderSize(0) 
    lhcbStyle.SetStatFont(lhcbFont) 
    lhcbStyle.SetStatFontSize(0.05) 
    lhcbStyle.SetStatX(0.9) 
    lhcbStyle.SetStatY(0.9) 
    lhcbStyle.SetStatW(0.25) 
    lhcbStyle.SetStatH(0.15) 
    #put tick marks on top and RHS of plots
    lhcbStyle.SetPadTickX(1) 
    lhcbStyle.SetPadTickY(1) 

    #histogram divisions: only 5 in x to avoid label overlaps
    lhcbStyle.SetNdivisions(505,"x") 
    lhcbStyle.SetNdivisions(510,"y") 

    #define style for text
    lhcbLabel = TText() 
    lhcbLabel.SetTextFont(lhcbFont) 
    lhcbLabel.SetTextColor(1) 
    lhcbLabel.SetTextSize(0.04) 
    lhcbLabel.SetTextAlign(12) 

    #define style of latex text
    lhcbLatex = TLatex() 
    lhcbLatex.SetTextFont(lhcbFont) 
    lhcbLatex.SetTextColor(1) 
    lhcbLatex.SetTextSize(0.04) 
    lhcbLatex.SetTextAlign(12) 

    #set this style
    gROOT.SetStyle("lhcbStyle") 
    gROOT.ForceStyle() 

    return lhcbStyle


#routine to print 'LHCb', 'LHCb Preliminary' on plots 
#options: optLR=L (top left) / R (top right) of plots
#         optPrelim= Final (LHCb), Prelim (LHCb Preliminary), Other
#         optText= text printed if 'Other' specified
def printLHCb(style,optLR="L",optPrelim="Final",optText=""):
    
    if (optLR=="R"):
        lhcbName = TPaveText(0.70 - style.GetPadRightMargin(),
                             0.75 - style.SetPadTopMargin(0.05),
                             0.95 - style.GetPadRightMargin(),
                             0.85 - style.SetPadTopMargin(0.05),
                             "BRNDC") 
  
    elif (optLR=="L"):
        lhcbName = TPaveText(style.GetPadLeftMargin() + 0.05,
                             0.87 - style.GetPadTopMargin(),
                             style.GetPadLeftMargin() + 0.30,
                             0.95 - style.GetPadTopMargin(),
                             "BRNDC") 
        
    else:
        print "printLHCb: option unknown" 
  
    if (optPrelim=="Final"):
        lhcbName.AddText("LHCb") 
  
    elif (optPrelim=="Prelim"):
        lhcbName.AddText("#splitline{LHCb}{#scale[1.0]{Preliminary}}")
        
    elif (optPrelim=="Other"):
        lhcbName.AddText(optText) 

    if optText:
        lhcbName.AddText(optText)

    else:
      print "printLHCb: option unknown " 

    lhcbName.SetFillColor(0) 
    lhcbName.SetTextAlign(12) 
    lhcbName.SetBorderSize(0) 
    #lhcbName.Draw()

    return lhcbName
    


