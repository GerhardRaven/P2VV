import RootStyle

from ROOT import gROOT, gStyle, TStyle, TCanvas, TH2D, TGraphErrors, TPaveText, TMarker, TLine, TPad, TFile
MyStyle = RootStyle.MyStyle()
gROOT.SetStyle(MyStyle.GetName())
gROOT.ForceStyle()
gStyle.UseCurrentStyle()

smallregion = False
inverted = True
if smallregion:
    jobbegin = 256
    jobend = 295

else:
    if inverted:
        jobbegin = 296
        jobend = 335
    else:
        jobbegin = 216
        jobend = 255


joblist = range(jobbegin,jobend+1)

npoints = 40
path = '/project/bfys/dveijk/gangadir/workspace/dveijk/LocalXML'

valuelist = []

filename = 'dGamma_phi_s_LL.txt'

for job in joblist:
    for subjob in range(0,npoints):
        filepath = path+'/%s/%s/output/'%(str(job),str(subjob))+filename
        #print filepath
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
                valuelist.append(line)
        llFile.close()

if smallregion:
    param1_min = -0.4
    param1_max = 0.4
    param2_min = 0
    param2_max = 0.2
else:
    param1_min = -1
    param1_max = 4
    param2_min = -0.25
    param2_max = 0.25
    
#Fill the TH2D
ProfileLikelihood = TH2D('ProfileLikelihood','ProfileLikelihood',npoints,param1_min,param1_max,npoints,param2_min,param2_max)
for v in valuelist:
    xstep = int(v[0])
    ystep = int(v[1])
    xvarval = float(v[2])
    yvarval = float(v[3])
    refLL = float(v[4])
    currentLL = float(v[5])

    dLL = currentLL - refLL
    ProfileLikelihood.SetBinContent(xstep+1,ystep+1,dLL)

if smallregion:
    from array import array
    colors = [2,3,4]
    colours = array('i',colors)
    gStyle.SetPalette(3,colours)

    gStyle.SetOptStat(0)

    #contours = [2.30/2.,4.61/2.,5.99/2.,9.21/2.]
    contours = [2.30/2.,4.61/2.,5.99/2.]
    contourarray = array('d',contours)

    Canvas = TCanvas('Canvas','Canvas',972,600)
    ProfileLikelihood.SetContour(len(contours),contourarray)
    ProfileLikelihood.Draw('CONT LIST')
    #The 'CONT LIST' is a trick to store the contours as TGraphs, such that the style can be modified later, if you don't want, you can use
    #ProfileLikelihood.Draw('CONT1')
    f = TFile('Profile_smallregion.root','recreate')
    ProfileLikelihood.Write()
    f.Close()
    
    Canvas.Update()
    mycontours = gROOT.GetListOfSpecials().FindObject("contours")
    ncontours     = mycontours.GetSize()

    contour1 = mycontours.At(0).First()
    contour2 = mycontours.At(1).First()
    contour3 = mycontours.At(2).First()

    contour1.SetLineColor(2)
    contour1.SetLineStyle(1)
    contour1.SetLineWidth(2)

    contour2.SetLineColor(3)
    contour2.SetLineStyle(2)
    contour2.SetLineWidth(2)

    contour3.SetLineColor(4)
    contour3.SetLineStyle(3)
    contour3.SetLineWidth(2)

    profilecanvas = TCanvas('profilecanvas','profilecanvas',972,600)
    pad = TPad('pad','pad',0,0,1,1)
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.SetFrameFillStyle(4000)
    pad.Draw()
    pad.cd()
    frame = pad.DrawFrame(param1_min,param2_min,param1_max,param2_max)
    frame.GetXaxis().SetTitle('#phi_{s}')
    frame.GetYaxis().SetTitle('#Delta#Gamma_{s} (ps^{-1})')
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetYaxis().SetTitleOffset(1.2)
    frame.SetTitle("")
    contour3.Draw('l same')
    contour1.Draw('l same')
    contour2.Draw('l same')

    #Draw SM point
    smphi = array("f", [-0.036])
    smphierr = array("f", [0.002])
    smdg = array("f", [0.087])
    smdgerr = array("f", [0.021])
    smpoint = TGraphErrors(1,smphi,smdg,smphierr,smdgerr)
    smpoint.Draw('same')

    #Draw best fit point
    marker1 = TMarker()
    marker1.SetMarkerStyle(8)
    marker1.SetMarkerColor(1)
    marker1.SetMarkerSize(0.9)
    marker1.DrawMarker(0.0025,0.1172)

    #Draw Legend
    myPT = TPaveText(0.1,0.16,0.3,0.20)
    #myPT.AddText(0.,0.,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}")
    myPT.SetShadowColor(0)
    myPT.SetBorderSize(0)
    myPT.SetFillStyle(0)
    myPT.SetTextSize(0.035)
    myPT.SetTextAlign(12)
    myPT.SetTextFont(42)
    myPT.Draw('same')

    #Draw Legend2
    myPT2 = TPaveText(0.05,0.01,0.35,0.07)
    myPT2.SetTextAlign(12)
    myPT2.AddText(0.35,0.9,"68% C.L.")
    myPT2.AddText(0.35,0.7,"90% C.L.")
    myPT2.AddText(0.35,0.5,"95% C.L.")
    myPT2.AddText(0.35,0.3,"SM expectation")
    myPT2.AddText(0.35,0.1,"Best fit point")
    myPT2.SetShadowColor(0)
    myPT2.SetFillStyle(0)
    myPT2.SetTextSize(0.035)
    myPT2.SetTextFont(42)

    myPT2.Draw('same')

    marker0 = TMarker()
    marker0.SetMarkerStyle(2)
    marker0.SetMarkerColor(1)
    marker0.SetMarkerSize(1.4)
    marker0.DrawMarker(0.1,0.029)

    marker2 = TMarker()
    marker2.SetMarkerStyle(8)
    marker2.SetMarkerColor(1)
    marker2.SetMarkerSize(0.9)
    marker2.DrawMarker(0.1,0.018)

    line1 = TLine(0.09,0.040,0.11,0.040)
    line1.SetLineColor(4)
    line1.SetLineStyle(3)
    line1.SetLineWidth(2)
    line1.Draw('same')

    line2 = TLine(0.09,0.0525,0.11,0.0525)
    line2.SetLineColor(3)
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)
    line2.Draw('same')

    line3 = TLine(0.09,0.065,0.11,0.065)
    line3.SetLineColor(2)
    line3.SetLineStyle(1)
    line3.SetLineWidth(2)
    line3.Draw('same')

    line4 = TLine(0.,0.,0,0.2)
    line4.SetLineColor(17)
    line4.SetLineStyle(7)
    #line4.SetLineWidth(1)
    line4.Draw('same')

    profilecanvas.SaveAs('2Dprofile_smallregion.eps')

else:
    from array import array
    colors = [2,3,4]
    colours = array('i',colors)
    gStyle.SetPalette(3,colours)

    gStyle.SetOptStat(0)

    #contours = [2.30/2.,4.61/2.,5.99/2.,9.21/2.]
    contours = [2.30/2.,4.61/2.,5.99/2.]
    contourarray = array('d',contours)

    Canvas = TCanvas('Canvas','Canvas',972,600)
    ProfileLikelihood.SetContour(len(contours),contourarray)
    ProfileLikelihood.Draw('CONT1')
    ProfileLikelihood.GetXaxis().SetTitle('#phi_{s}')
    ProfileLikelihood.GetYaxis().SetTitle('#Delta#Gamma_{s} (ps^{-1})')
    ProfileLikelihood.GetXaxis().SetTitleOffset(1.2)
    ProfileLikelihood.GetYaxis().SetTitleOffset(1.2)
    ProfileLikelihood.SetTitle("")

    if inverted:
        f = TFile('Profile_inverted.root','recreate')        
    else:
        f = TFile('Profile.root','recreate')        

    ProfileLikelihood.Write()
    f.Close()
    
    #Draw SM point
    smphi = array("f", [-0.036])
    smphierr = array("f", [0.002])
    smdg = array("f", [0.087])
    smdgerr = array("f", [0.021])
    smpoint = TGraphErrors(1,smphi,smdg,smphierr,smdgerr)
    smpoint.Draw('same')

    #Draw Legend
    myPT = TPaveText(2,0.16,3.5,0.20)
    #myPT.AddText(0.,0.,"#splitline{#font[42]{LHCb preliminary}}{#font[42]{#sqrt{s} = 7 TeV, L = 1.0 fb^{-1}}}")
    myPT.SetShadowColor(0)
    myPT.SetBorderSize(0)
    myPT.SetFillStyle(0)
    myPT.SetTextSize(0.035)
    myPT.SetTextAlign(12)
    myPT.SetTextFont(42)
    myPT.Draw('same')

    #Draw Legend2
    myPT2 = TPaveText(0.0,-0.15,1.5,-0.05)
    myPT2.SetTextAlign(12)
    myPT2.AddText(0.0,0.9,"Confidence levels")
    myPT2.AddText(0.35,0.7,"68% C.L.")
    myPT2.AddText(0.35,0.5,"90% C.L.")
    myPT2.AddText(0.35,0.3,"95% C.L.")
    myPT2.AddText(0.35,0.1,"SM expectation")
    myPT2.SetShadowColor(0)
    myPT2.SetFillStyle(0)
    myPT2.SetTextSize(0.035)
    myPT2.SetTextFont(42)

    myPT2.Draw('same')

    marker0 = TMarker()
    marker0.SetMarkerStyle(2)
    marker0.SetMarkerColor(1)
    marker0.SetMarkerSize(1.4)
    marker0.DrawMarker(0.25,-0.14)

    line1 = TLine(0.15,-0.12,0.35,-0.12)
    line1.SetLineColor(4)
    line1.SetLineWidth(2)
    line1.Draw('same')

    line2 = TLine(0.15,-0.10,0.35,-0.10)
    line2.SetLineColor(3)
    line2.SetLineWidth(2)
    line2.Draw('same')

    line3 = TLine(0.15,-0.08,0.35,-0.08)
    line3.SetLineColor(2)
    line3.SetLineWidth(2)
    line3.Draw('same')

    if inverted:
        Canvas.SaveAs('2Dprofile_inverted.eps')
    else:
        Canvas.SaveAs('2Dprofile.eps')
    
