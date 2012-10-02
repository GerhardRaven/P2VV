import RootStyle
from ROOT import gROOT, gStyle, TStyle, TGraph, TCanvas, TPaveText
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

#AparPhase 60 bins
#jobnr = 86
#varname = "AparPhase"
#xtitle = "#delta_{||}"
#nsteps = 60

#deltams wide
#jobnr = 337
#varname = "dM"
#xtitle = "#Delta m_{s} [1/ps]"
#nsteps = 40

#deltams narrow
#jobnr = 338
#varname = "dM"
#xtitle = "#Delta m_{s} [1/ps]"
#nsteps = 40

#f_S crosscheck
jobnr = 395
varname = "f_S"
xtitle = "f_{S}"
nsteps = 40

path = '/project/bfys/dveijk/gangadir/workspace/dveijk/LocalXML'
filename = '%s_LL.txt'%(varname)

llref = []
varval = []
ll = []

for i in range(0,nsteps):
    filepath = path+'/%s/%s/output/'%(str(jobnr),str(i))+filename
    llFile = open(filepath)
    while True :
        # read next line
        line = llFile.readline()
        if not line : break

        # check for empty or comment lines
        line = line.strip()
        if not line or line[0] == '#' :
            continue
        else:

            # check moment format
            line = line.split()

            varval.append(float(line[1]))
            llref.append(float(line[2]))
            ll.append(float(line[3]))

llFile.close()

dll = []
for i in range(0,nsteps):
    dll.append(abs(ll[i]-llref[i]))

from array import array
dll_arr = array('f',dll)
varval_arr = array('f',varval)

canvas = TCanvas("%s"%(xtitle),"%s"%(xtitle),972,600)
graph = TGraph(nsteps,varval_arr,dll_arr)

graph.SetTitle("")
graph.GetXaxis().SetTitle(xtitle)
graph.GetXaxis().SetTitleOffset(1.2)
graph.GetYaxis().SetTitle('DLL')
graph.SetMarkerStyle(20)
#graph.SetMarkerSize(1)

graph.Draw("apc")

#Draw Legend
PT = TPaveText(0.4,0.75,0.9,0.85,"NDC")
#PT.AddText(0.,0.,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}")
PT.SetShadowColor(0)
PT.SetBorderSize(0)
PT.SetFillStyle(0)
PT.SetTextSize(0.045)
PT.SetTextAlign(12)
PT.SetTextFont(42)
PT.Draw('same')

canvas.SaveAs("%s.eps"%(varname))
