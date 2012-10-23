import RootStyle
from ROOT import gROOT,gStyle,TStyle,TCanvas,TGraph, TPaveText
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

jobbegin = 75
jobend = 84

joblist = range(jobbegin,jobend+1)

varnamelist = ["A0Mag2","ASPhase","AparPhase","AperpMag2","AperpPhase","Gamma","dGamma","dM","f_S","phi_s"]
xtitlelist = ["|A_{0}|^{2}","#delta_{S}","#delta_{||}","|A_{perp}|^{2}","#delta_{perp}","#Gamma_{s} [1/ps]","#Delta#Gamma_{s} [1/ps]","#Delta m_{s}[1/ps]","f_{S}","#phi_{s}"]

nsteps = 40
path = '/project/bfys/dveijk/gangadir/workspace/dveijk/LocalXML'

for n,j in enumerate(joblist):
    llref = []
    varval = []
    ll = []
    filename = '%s_LL.txt'%(varnamelist[n])
    for i in range(0,nsteps):
        filepath = path+'/%s/%s/output/'%(str(j),str(i))+filename
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

    vars()["canvas_%s"%(varnamelist[n])] = TCanvas("%s"%(varnamelist[n]),"%s"%(varnamelist[n]),972,600)
    vars()["graph_%s"%(varnamelist[n])] = TGraph(nsteps,varval_arr,dll_arr)

    vars()["graph_%s"%(varnamelist[n])].GetXaxis().SetTitle(xtitlelist[n])
    vars()["graph_%s"%(varnamelist[n])].GetXaxis().SetTitleOffset(1.2)
    vars()["graph_%s"%(varnamelist[n])].GetYaxis().SetTitle('DLL')
    vars()["graph_%s"%(varnamelist[n])].SetTitle("")
    vars()["graph_%s"%(varnamelist[n])].SetMarkerStyle(20)
    #vars()["graph_%s"%(varnamelist[n])].SetMarkerSize(1)

    vars()["graph_%s"%(varnamelist[n])].Draw("apl")

    #Draw Legend
    vars()["PT_%s"%(varnamelist[n])] = TPaveText(0.4,0.75,0.6,0.85,"NDC")
    #vars()["PT_%s"%(varnamelist[n])].AddText(0.,0.,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}")
    vars()["PT_%s"%(varnamelist[n])].SetShadowColor(0)
    vars()["PT_%s"%(varnamelist[n])].SetBorderSize(0)
    vars()["PT_%s"%(varnamelist[n])].SetFillStyle(0)
    vars()["PT_%s"%(varnamelist[n])].SetTextSize(0.05)
    vars()["PT_%s"%(varnamelist[n])].SetTextAlign(12)
    vars()["PT_%s"%(varnamelist[n])].SetTextFont(42)
    vars()["PT_%s"%(varnamelist[n])].Draw('same')

    vars()["canvas_%s"%(varnamelist[n])].SaveAs("%s.eps"%(varnamelist[n]))
