model = 'lamb_phi'
scanParFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/NLLPlots/jobOutput/%s_%s/NLLVals_%s_1000???.par'
plotFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/NLLPlots/%s_plots.pdf' % model
nPointsPara = 1000
nllRange = ( 0., 6. )

from math import pi
if model == 'phi' :
    scanPars = [  'phiCP'
                , 'Gamma', 'dGamma', 'dM'
                , 'A0Mag2', 'AperpMag2'
                , 'f_S_bin0', 'f_S_bin1', 'f_S_bin2', 'f_S_bin3', 'f_S_bin4', 'f_S_bin5'
                , 'AparPhase', 'AperpPhase'
                , 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
               ]

    parSettings = dict(  phiCP              = dict( name = '#phi_{s}', min = -0.31, max = 0.19
                                                   , mean = -0.056308557, error = 0.049043162 )
                       , Gamma              = dict( name = '#Gamma_{s}', min = 0.643, max = 0.676
                                                   , mean = 0.65909601, error = 0.0031105883 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s}', min = 0.033, max = 0.124
                                                   , mean = 0.078519243, error = 0.0091211774 )
                       , dM                 = dict( name = '#Deltam_{s}', min = 17.44, max = 18.01
                                                   , mean = 17.696881, error = 0.060104838 )
                       , A0Mag2             = dict( name = '|A_{#kern[0.4]{0}}|^{2}', min = 0.506, max = 0.541
                                                   , mean = 0.52362287, error = 0.0034390042 )
                       , AperpMag2          = dict( name = '|A_{#perp#kern[1.8]{ }}|^{2}', min = 0.227, max = 0.276
                                                   , mean = 0.2512126, error = 0.0049162961 )
                       , f_S_bin0           = dict( name = 'F_{S1}', min = 0.16, max = 0.70
                                                   , mean = 0.42608908, error = 0.054057373 )
                       , f_S_bin1           = dict( name = 'F_{S2}', min = 1.e-5, max = 0.156
                                                   , mean = 0.058932781, error = 0.017982136 )
                       , f_S_bin2           = dict( name = 'F_{S3}', min = 1.e-5, max = 0.042
                                                   , mean = 0.0095124401, error = 0.0068361602 )
                       , f_S_bin3           = dict( name = 'F_{S4}', min = 1.e-5, max = 0.039
                                                   , mean = 0.007901675, error = 0.005714802 )
                       , f_S_bin4           = dict( name = 'F_{S5}', min = 1.e-5, max = 0.126
                                                   , mean = 0.045102864, error = 0.016058515 )
                       , f_S_bin5           = dict( name = 'F_{S6}', min = 0.064, max = 0.319
                                                   , mean = 0.19240477, error = 0.025483256 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0}', min = 2.64, max = 3.88
                                                   , mean = 3.2638391, error = 0.12336781 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0}', min = 2.38, max = 3.82
                                                   , mean = 3.0426437, error = 0.16087983 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.83992791, error = 0.19957076 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 2.1505793, error = 0.28832208 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.4827505, error = 0.23097498 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.40271698, error = 0.21455899 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.62348152, error = 0.17595703 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.90087493, error = 0.13933952 )
                      )

elif model == 'lamb_phi' :
    scanPars = [  'phiCP', 'lambdaCP'
                , 'Gamma', 'dGamma', 'dM'
                , 'A0Mag2', 'AperpMag2'
                , 'f_S_bin0', 'f_S_bin1', 'f_S_bin2', 'f_S_bin3', 'f_S_bin4', 'f_S_bin5'
                , 'AparPhase', 'AperpPhase'
                , 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
               ]

    parSettings = dict(  phiCP              = dict( name = '#phi_{s}', min = -0.31, max = 0.19
                                                   , mean = -0.057095686, error = 0.049722504 )
                       , lambdaCP           = dict( name = '|#lambda_{s}|', min = 0.869, max = 1.057
                                                   , mean = 0.96269479, error = 0.018785372 )
                       , Gamma              = dict( name = '#Gamma_{s}', min = 0.643, max = 0.676
                                                   , mean = 0.65916402, error = 0.0031062888 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s}', min = 0.033, max = 0.124
                                                   , mean = 0.078511129, error = 0.0091307288 )
                       , dM                 = dict( name = '#Deltam_{s}', min = 17.44, max = 18.01
                                                   , mean = 17.723016, error = 0.056704725 )
                       , A0Mag2             = dict( name = '|A_{#kern[0.4]{0}}|^{2}', min = 0.506, max = 0.541
                                                   , mean = 0.52367494, error = 0.0034391553 )
                       , AperpMag2          = dict( name = '|A_{#perp#kern[1.8]{ }}|^{2}', min = 0.227, max = 0.276
                                                   , mean = 0.25121522, error = 0.0049002128 )
                       , f_S_bin0           = dict( name = 'F_{S1}', min = 0.16, max = 0.70
                                                   , mean = 0.42628995, error = 0.054010232 )
                       , f_S_bin1           = dict( name = 'F_{S2}', min = 1.e-5, max = 0.156
                                                   , mean = 0.05875113, error = 0.017556867 )
                       , f_S_bin2           = dict( name = 'F_{S3}', min = 1.e-5, max = 0.042
                                                   , mean = 0.0095565319, error = 0.00656844 )
                       , f_S_bin3           = dict( name = 'F_{S4}', min = 1.e-5, max = 0.039
                                                   , mean = 0.0094082209, error = 0.0058415528 )
                       , f_S_bin4           = dict( name = 'F_{S5}', min = 1.e-5, max = 0.126
                                                   , mean = 0.04816402, error = 0.015469857 )
                       , f_S_bin5           = dict( name = 'F_{S6}', min = 0.064, max = 0.319
                                                   , mean = 0.19205916, error = 0.02548107 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0}', min = 2.64, max = 3.88
                                                   , mean = 3.2567147, error = 0.12416783 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0}', min = 2.38, max = 3.82
                                                   , mean = 3.0985034, error = 0.14397861 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.84325523, error = 0.19890243 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 2.1474486, error = 0.28112294 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.48300441, error = 0.22314765 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.36420991, error = 0.18210145 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.59230381, error = 0.15584043 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.9081375, error = 0.14030522 )
                      )
else :
    scanPars=[  'phiCPAv', 'phiCPRel_Apar', 'phiCPRel_AperpApar', 'phiCPRel_AS', 'CCPAv', 'CCPRel_Apar', 'CCPRel_Aperp', 'CCPAv_AS'
              , 'Gamma', 'dGamma', 'dM'
              , 'avA02', 'avAperp2'
              , 'avf_S_bin0', 'avf_S_bin1', 'avf_S_bin2', 'avf_S_bin3', 'avf_S_bin4', 'avf_S_bin5'
              , 'AparPhase', 'AperpPhase'
              , 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
             ]

    parSettings = dict(  phiCPAv            = dict( name = '#phi_{s}^{av}', min = -0.30, max = 0.21
                                                   , mean = -0.046577055, error = 0.050931782 )
                       , phiCPRel_Apar      = dict( name = '#Delta#phi_{s}^{#parallel}', min = -0.23, max = 0.19
                                                   , mean = -0.018693067, error = 0.042595888 )
                       , phiCPRel_AperpApar = dict( name = '#Delta#phi_{s}^{#perp\'}', min = -0.15, max = 0.14
                                                   , mean = -0.0026496495, error = 0.028615121 )
                       , phiCPRel_AS        = dict( name = '#Delta#phi_{s}^{S}', min = -0.30, max = 0.33
                                                   , mean = 0.014648547, error = 0.062238538 )
                       , CCPAv              = dict( name = 'C_{s}^{av}', min = -0.20, max = 0.19
                                                   , mean = -0.0063173106, error = 0.03865467 )
                       , CCPRel_Apar        = dict( name = '#DeltaC_{s}^{#parallel}', min = -0.63, max = 0.58
                                                   , mean = -0.024660592, error = 0.1217666 )
                       , CCPRel_Aperp       = dict( name = '#DeltaC_{s}^{#perp}', min = -0.77, max = 0.86
                                                   , mean = 0.043741257, error = 0.16232504 )
                       , CCPAv_AS           = dict( name = 'C_{s}^{avS}', min = -0.10, max = 0.22
                                                   , mean = 0.059909153, error = 0.032151267 )
                       , Gamma              = dict( name = '#Gamma_{s}', min = 0.643, max = 0.676
                                                   , mean = 0.65910705, error = 0.0031164336 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s}', min = 0.033, max = 0.124
                                                   , mean = 0.078376082, error = 0.009155929 )
                       , dM                 = dict( name = '#Deltam_{s}', min = 17.44, max = 18.01
                                                   , mean = 17.696386, error = 0.062072481 )
                       , avA02              = dict( name = 'A_{#kern[0.4]{0}}^{av2}', min = 0.506, max = 0.541
                                                   , mean = 0.52363993, error = 0.0034434589 )
                       , avAperp2           = dict( name = 'A_{#perp#kern[1.8]{ }}^{av2}', min = 0.227, max = 0.276
                                                   , mean = 0.25125825, error = 0.0049294865 )
                       , avf_S_bin0         = dict( name = 'F_{S1}^{av}', min = 0.16, max = 0.70
                                                   , mean = 0.42437739, error = 0.054153196 )
                       , avf_S_bin1         = dict( name = 'F_{S2}^{av}', min = 1.e-5, max = 0.156
                                                   , mean = 0.057217699, error = 0.01766972 )
                       , avf_S_bin2         = dict( name = 'F_{S3}^{av}', min = 1.e-5, max = 0.042
                                                   , mean = 0.0085864882, error = 0.0065480008 )
                       , avf_S_bin3         = dict( name = 'F_{S4}^{av}', min = 1.e-5, max = 0.039
                                                   , mean = 0.0092803772, error = 0.0056304117 )
                       , avf_S_bin4         = dict( name = 'F_{S5}^{av}', min = 1.e-5, max = 0.126
                                                   , mean = 0.047881399, error = 0.015387567 )
                       , avf_S_bin5         = dict( name = 'F_{S6}^{av}', min = 0.064, max = 0.319
                                                   , mean = 0.19101304, error = 0.025535031 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0}', min = 2.64, max = 3.88
                                                   , mean = 3.2461245, error = 0.13231374 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0}', min = 2.38, max = 3.82
                                                   , mean = 3.0369133, error = 0.16474461 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.86518475, error = 0.20273115 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 2.1238832, error = 0.30849018 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = 0.52856132, error = 0.25208206 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.35481346, error = 0.18334874 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.58572974, error = 0.15818632 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }}', min = -pi, max = pi
                                                   , mean = -0.90582695, error = 0.1473001 )
                      )

from P2VV.Load import LHCbStyle
from ROOT import gStyle
gStyle.SetColorModelPS(1)

from glob import glob
from array import array
from ROOT import TCanvas, TGraph, kBlack, kRed, kBlue, kFullDotLarge
graphs = { }
canvs = dict( dummy = TCanvas('dummy_canv') )
canvs['dummy'].Print( plotFilePath + '[' )
for parIt, par in enumerate(scanPars) :
    parFiles = glob( scanParFilePath % ( model, par, par ) )
    if not parFiles : continue

    parVals  = array( 'd', [ ] )
    NLLVals  = array( 'd', [ ] )
    profVals = array( 'd', [ ] )
    for filePath in parFiles :
        parFile = open(filePath)
        while True :
            line = parFile.readline()
            if not line : break
            line = line.split()

            assert line[0] == par and line[4] == 'NLL' and line[8] == 'profiled' and line[9] == 'NLL'

            pos = 0
            parVal = float(line[2])
            if par.startswith('ASOddPhase') and parVal > pi : parVal -= 2. * pi
            for val in parVals :
                if val > parVal : break
                pos += 1
            parVals.insert( pos, parVal )
            NLLVals.insert( pos, float(line[6]) )
            profVals.insert( pos, float(line[11]) )

        parFile.close()

    NLLMin  = min(NLLVals)
    profMin = min(profVals)
    for it in range( len(NLLVals) )  : NLLVals[it]  -= NLLMin
    for it in range( len(profVals) ) : profVals[it] -= profMin

    parMin  = parSettings[par]['min']
    parMax  = parSettings[par]['max']
    parMean = parSettings[par]['mean']
    parErr  = parSettings[par]['error']
    parValsPara = array( 'd', [ parMin + float(it) / float(nPointsPara - 1) * ( parMax - parMin ) for it in range(nPointsPara) ] )
    NLLValsPara = array( 'd', [ 0.5 * ( ( val - parMean ) / parErr )**2 for val in parValsPara ] )

    graphs[par] = (  TGraph( len(parValsPara), parValsPara, NLLValsPara )
                   , TGraph( len(parVals), parVals, NLLVals )
                   , TGraph( len(parVals), parVals, profVals )
                  )
    canvs[par] = TCanvas( '%s_canv' % par )
    canvs[par].SetLeftMargin(0.18)
    canvs[par].SetRightMargin(0.05)
    canvs[par].SetBottomMargin(0.20)
    canvs[par].SetTopMargin(0.05)
    for it, graph in enumerate(graphs[par]) :
        graph.SetLineWidth( 3 if it > 0 else 2 )
        graph.SetMarkerStyle(kFullDotLarge)
        graph.SetMarkerSize(0.6)
        graph.SetLineColor( kBlue if it == 1 else kRed if it == 2 else kBlack )
        graph.SetMarkerColor( kBlue if it == 1 else kRed if it == 2 else kBlack )
        graph.SetMinimum( nllRange[0] )
        graph.SetMaximum( nllRange[1] )
        graph.GetXaxis().SetLimits( parMin, parMax )
        graph.GetXaxis().SetTitle(parSettings[par]['name'])
        graph.GetYaxis().SetTitle('#Delta log(L)')
        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset(0.9)
        graph.Draw( 'AL' if it == 0 else 'SAMES PC' )
    canvs[par].Print(plotFilePath)
canvs['dummy'].Print( plotFilePath + ']' )
