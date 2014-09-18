from ROOT import kBlack, kGray, kRed, kMagenta, kBlue, kCyan, kGreen, kYellow, kFullDotLarge, kSolid, kDashed
model = 'polarDep'
scanParFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/NLLPlots/jobOutput/%s_%s/NLLVals_%s_100[2-3]???.par'
#scanParFilePath = '/project/bfys/jleerdam/softDevel/P2VV2/test/NLLPlots/jobOutput/%s_%s/NLLVals_%s_1002???.par'
plotFilePath = 'NLL_plots_%s.pdf' % model #'/project/bfys/jleerdam/softDevel/P2VV2/test/NLLPlots/%s_plots.pdf' % model
nPointsPara = 1000
NLLRange = ( 0., 14.)
NLLLabel = '-#Deltalog(L)'

SWavePhases = [ 'ASOddPhase_bin%d' % bin for bin in range(6) ]
SWavePhasesLabel = '#delta_{Si} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
SWavePlotLabels = { 0 : ( -0.5, 4.5, '990 - 1008 MeV' ), 1 : ( -0.5, 4.5, '1008 - 1016 MeV' ), 2 : ( -0.5, 4.5, '1016 - 1020 MeV' )
                   , 3 : ( 2.5, 4.5, '1020 - 1024 MeV' ), 4 : ( 2.5, 4.5, '1024 - 1032 MeV' ), 5 : ( 2.5, 4.5, '1032 - 1050 MeV' ) }
SWaveColours = { 0 : kCyan + 3, 1 : kBlue + 2, 2 : kMagenta + 2, 3 : kRed + 1, 4 : kYellow + 3, 5 : kGreen + 3 }
SWaveStyles = { 0 : 9, 1 : 7, 2 : 0, 3 : 0, 4 : 7, 5 : 9 }
SWaveLabelPos = { 0 : ( 0.8, -0.55 ), 1 : ( 2.5, 0.4 ), 2 : ( 0.33, -0.55 )
                 , 3 : ( -0.45, -0.55 ), 4 : ( -0.85, -0.55 ), 5 : ( -1.5, 0.4 ) }
SWaveNLLRange = ( 0., 8. )

from math import pi
if model == 'phi' :
    scanPars = [  'phiCP'
                , 'Gamma', 'dGamma', 'dM'
                , 'A0Mag2', 'AperpMag2'
                , 'f_S_bin0', 'f_S_bin1', 'f_S_bin2', 'f_S_bin3', 'f_S_bin4', 'f_S_bin5'
                , 'AparPhase', 'AperpPhase'
                , 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
               ]
    parSettings = dict(  #phiCP              = dict( name = '#phi_{s} [rad]', min = -0.35, max = 0.25
                         #                          , mean = -0.056308557, error = 0.049043162 )
                         phiCP              = dict( name = '#phi_{s} [rad]', min = -0.20, max = 0.40, offset = 0.1301
                                                   , mean = -0.056308557, error = 0.049043162 )
                       , Gamma              = dict( name = '#Gamma_{s} [ps^{-1}]', min = 0.64, max = 0.68
                                                   , mean = 0.65909601, error = 0.0031105883 )
                       #, dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.02, max = 0.13
                       #                            , mean = 0.078519243, error = 0.0091211774 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.03, max = 0.14, offset = 0.0100
                                                   , mean = 0.078519243, error = 0.0091211774 )
                       , dM                 = dict( name = '#Deltam_{s} [ps^{-1}]', min = 17.3, max = 18.1
                                                   , mean = 17.696881, error = 0.060104838 )
                       , A0Mag2             = dict( name = '|A_{#kern[0.4]{0}}|^{2}', min = 0.50, max = 0.55
                                                   , mean = 0.52362287, error = 0.0034390042 )
                       , AperpMag2          = dict( name = '|A_{#perp#kern[1.8]{ }}|^{2}', min = 0.22, max = 0.28
                                                   , mean = 0.2512126, error = 0.0049162961 )
                       , f_S_bin0           = dict( name = 'F_{S1}', min = 0.15, max = 0.70
                                                   , mean = 0.42608908, error = 0.054057373 )
                       , f_S_bin1           = dict( name = 'F_{S2}', min = 0., max = 0.16
                                                   , mean = 0.058932781, error = 0.017982136 )
                       , f_S_bin2           = dict( name = 'F_{S3}', min = 0., max = 0.06
                                                   , mean = 0.0095124401, error = 0.0068361602 )
                       , f_S_bin3           = dict( name = 'F_{S4}', min = 0., max = 0.06
                                                   , mean = 0.007901675, error = 0.005714802 )
                       , f_S_bin4           = dict( name = 'F_{S5}', min = 0., max = 0.15
                                                   , mean = 0.045102864, error = 0.016058515 )
                       , f_S_bin5           = dict( name = 'F_{S6}', min = 0.06, max = 0.35
                                                   , mean = 0.19240477, error = 0.025483256 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0} [rad]', min = 2.6, max = 3.7
                                                   , mean = 3.2638391, error = 0.12336781 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0} [rad]'
                                                   , min = 2.0, max = 4.0
                                                   , mean = 3.0426437, error = 0.16087983 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0., max = 3.
                                                   , mean = 0.83992791, error = 0.19957076 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0., max = 3.
                                                   , mean = 2.1505793, error = 0.28832208 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0., max = 3.
                                                   , mean = 0.4827505, error = 0.23097498 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -3., max = 0.
                                                   , mean = -0.40271698, error = 0.21455899 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -3., max = 0.
                                                   , mean = -0.62348152, error = 0.17595703 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -3., max = 0.
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

    parSettings = dict(  #phiCP              = dict( name = '#phi_{s} [rad]', min = -0.34, max = 0.34
                         #                          , mean = -0.057095686, error = 0.049722504 )
                         phiCP              = dict( name = '#phi_{s} [rad]', min = -0.20, max = 0.40, offset = 0.1301
                                                   , mean = -0.057095686, error = 0.049722504 )
                       , lambdaCP           = dict( name = '|#lambda_{s}|', min = 0.80, max = 1.20
                                                   , mean = 0.96269479, error = 0.018785372 )
                       , Gamma              = dict( name = '#Gamma_{s} [ps^{-1}]', min = 0.64, max = 0.68
                                                   , mean = 0.65916402, error = 0.0031062888 )
                       #, dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.02, max = 0.13
                       #                            , mean = 0.078511129, error = 0.0091307288 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.03, max = 0.14, offset = 0.0100
                                                   , mean = 0.078511129, error = 0.0091307288 )
                       , dM                 = dict( name = '#Deltam_{s} [ps^{-1}]', min = 17.25, max = 18.10
                                                   , mean = 17.723016, error = 0.056704725 )
                       , A0Mag2             = dict( name = '|A_{#kern[0.4]{0}}|^{2}', min = 0.50, max = 0.55
                                                   , mean = 0.52367494, error = 0.0034391553 )
                       , AperpMag2          = dict( name = '|A_{#perp#kern[1.8]{ }}|^{2}', min = 0.22, max = 0.28
                                                   , mean = 0.25121522, error = 0.0049002128 )
                       , f_S_bin0           = dict( name = 'F_{S1}', min = 0.15, max = 0.70
                                                   , mean = 0.42628995, error = 0.054010232 )
                       , f_S_bin1           = dict( name = 'F_{S2}', min = 0., max = 0.16
                                                   , mean = 0.05875113, error = 0.017556867 )
                       , f_S_bin2           = dict( name = 'F_{S3}', min = 0., max = 0.06
                                                   , mean = 0.0095565319, error = 0.00656844 )
                       , f_S_bin3           = dict( name = 'F_{S4}', min = 0., max = 0.06
                                                   , mean = 0.0094082209, error = 0.0058415528 )
                       , f_S_bin4           = dict( name = 'F_{S5}', min = 0., max = 0.16
                                                   , mean = 0.04816402, error = 0.015469857 )
                       , f_S_bin5           = dict( name = 'F_{S6}', min = 0.06, max = 0.35
                                                   , mean = 0.19205916, error = 0.02548107 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0} [rad]', min = 2.6, max = 3.7
                                                   , mean = 3.2567147, error = 0.12416783 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0} [rad]'
                                                   , min = 2., max = 4.
                                                   , mean = 3.0985034, error = 0.14397861 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0.5, max = 1.5
                                                   , mean = 0.84325523, error = 0.19890243 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 1.5, max = 2.5
                                                   , mean = 2.1474486, error = 0.28112294 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0., max = 1.
                                                   , mean = 0.48300441, error = 0.22314765 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1., max = 0.
                                                   , mean = -0.36420991, error = 0.18210145 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1., max = 0.
                                                   , mean = -0.59230381, error = 0.15584043 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1.5, max = -0.5
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

    parSettings = dict(#  phiCPAv            = dict( name = '#phi_{s}^{av} [rad]', min = -0.34, max = 0.34
                       #                            , mean = -0.046577055, error = 0.050931782 )
                       #, phiCPRel_Apar      = dict( name = '#Delta#phi_{s}^{#parallel} [rad]', min = -0.40, max = 0.40
                       #                            , mean = -0.018693067, error = 0.042595888 )
                       #, phiCPRel_AperpApar = dict( name = '#Delta#phi_{s}^{#perp}#kern[1]{\'} [rad]', min = -0.40, max = 0.40
                       #                            , mean = -0.0026496495, error = 0.028615121 )
                       #, phiCPRel_AS        = dict( name = '#Delta#phi_{s}^{S} [rad]', min = -0.70, max = 0.40
                       #                            , mean = 0.014648547, error = 0.062238538 )
                         phiCPAv            = dict( name = '#phi_{s}^{av} [rad]', min = -0.15, max = 0.45, offset = 0.1801
                                                   , mean = -0.046577055, error = 0.050931782 )
                       , phiCPRel_Apar      = dict( name = '#Delta#phi_{s}^{#parallel} [rad]', min = -0.50, max = 0.30, offset = -0.1089
                                                   , mean = -0.018693067, error = 0.042595888 )
                       , phiCPRel_AperpApar = dict( name = '#Delta#phi_{s}^{#perp}#kern[1]{\'} [rad]', min = -0.30, max = 0.50
                                                   , offset = 0.0960, mean = -0.0026496495, error = 0.028615121 )
                       , phiCPRel_AS        = dict( name = '#Delta#phi_{s}^{S} [rad]', min = -0.90, max = 0.20, offset = -0.2350
                                                   , mean = 0.014648547, error = 0.062238538 )
                       , CCPAv              = dict( name = 'C_{s}^{av}', min = -0.24, max = 0.24
                                                   , mean = -0.0063173106, error = 0.03865467 )
                       , CCPRel_Apar        = dict( name = '#DeltaC_{s}^{#parallel}', min = -0.80, max = 0.80
                                                   , mean = -0.024660592, error = 0.1217666 )
                       , CCPRel_Aperp       = dict( name = '#DeltaC_{s}^{#perp}', min = -1.0, max = 1.0
                                                   , mean = 0.043741257, error = 0.16232504 )
                       , CCPAv_AS           = dict( name = 'C_{s}^{avS}', min = -0.45, max = 0.45
                                                   , mean = 0.059909153, error = 0.032151267 )
                       , Gamma              = dict( name = '#Gamma_{s} [ps^{-1}]', min = 0.64, max = 0.68
                                                   , mean = 0.65910705, error = 0.0031164336 )
                       #, dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.02, max = 0.13
                       #                            , mean = 0.078376082, error = 0.009155929 )
                       , dGamma             = dict( name = '#Delta#Gamma_{s} [ps^{-1}]', min = 0.03, max = 0.14, offset = 0.0100
                                                   , mean = 0.078376082, error = 0.009155929 )
                       , dM                 = dict( name = '#Deltam_{s} [ps^{-1}]', min = 17.25, max = 18.10
                                                   , mean = 17.696386, error = 0.062072481 )
                       , avA02              = dict( name = '|A_{#kern[0.4]{0}}^{CP}|^{2}', min = 0.50, max = 0.55
                                                   , mean = 0.52363993, error = 0.0034434589 )
                       , avAperp2           = dict( name = '|A_{#perp#kern[1.8]{ }}^{CP}|^{2}', min = 0.22, max = 0.28
                                                   , mean = 0.25125825, error = 0.0049294865 )
                       , avf_S_bin0         = dict( name = 'F_{S1}^{av}', min = 0.15, max = 0.70
                                                   , mean = 0.42437739, error = 0.054153196 )
                       , avf_S_bin1         = dict( name = 'F_{S2}^{av}', min = 0., max = 0.15
                                                   , mean = 0.057217699, error = 0.01766972 )
                       , avf_S_bin2         = dict( name = 'F_{S3}^{av}', min = 0., max = 0.04
                                                   , mean = 0.0085864882, error = 0.0065480008 )
                       , avf_S_bin3         = dict( name = 'F_{S4}^{av}', min = 0., max = 0.04
                                                   , mean = 0.0092803772, error = 0.0056304117 )
                       , avf_S_bin4         = dict( name = 'F_{S5}^{av}', min = 0., max = 0.12
                                                   , mean = 0.047881399, error = 0.015387567 )
                       , avf_S_bin5         = dict( name = 'F_{S6}^{av}', min = 0.06, max = 0.35
                                                   , mean = 0.19101304, error = 0.025535031 )
                       , AparPhase          = dict( name = '#delta_{#parallel} - #delta_{0} [rad]', min = 2.6, max = 3.7
                                                   , mean = 3.2461245, error = 0.13231374 )
                       , AperpPhase         = dict( name = '#delta_{#kern[1.8]{#perp}#kern[1.8]{ }} - #delta_{0} [rad]'
                                                   , min = 2.0, max = 4.0
                                                   , mean = 3.0369133, error = 0.16474461 )
                       , ASOddPhase_bin0    = dict( name = '#delta_{S1} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0.5, max = 1.5
                                                   , mean = 0.86518475, error = 0.20273115 )
                       , ASOddPhase_bin1    = dict( name = '#delta_{S2} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 1.4, max = 2.4
                                                   , mean = 2.1238832, error = 0.30849018 )
                       , ASOddPhase_bin2    = dict( name = '#delta_{S3} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = 0., max = 1.
                                                   , mean = 0.52856132, error = 0.25208206 )
                       , ASOddPhase_bin3    = dict( name = '#delta_{S4} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1., max = 0.
                                                   , mean = -0.35481346, error = 0.18334874 )
                       , ASOddPhase_bin4    = dict( name = '#delta_{S5} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1., max = 0.
                                                   , mean = -0.58572974, error = 0.15818632 )
                       , ASOddPhase_bin5    = dict( name = '#delta_{S6} - #delta_{#kern[1.8]{#perp}#kern[1.8]{ }} [rad]'
                                                   , min = -pi, max = pi
                                                   #, min = -1.5, max = -0.5
                                                   , mean = -0.90582695, error = 0.1473001 )
                      )

from P2VV.Load import LHCbStyle
from ROOT import gStyle
#gStyle.SetColorModelPS(1)
#gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
#gStyle.SetLineStyleString( 7, ' 40 20'        )
#gStyle.SetLineStyleString( 9, ' 100 20'       )

from glob import glob
from array import array
from ROOT import TPad, TCanvas, TGraph, TLatex
graphs = { }
graphs0 = { }
graphs1 = { }
pads = { }
canvs = dict( dummy = TCanvas('dummy_canv') )
canvs['dummy'].Print( plotFilePath + '[' )
if SWavePhases :
    firstSWave = True
    canvs['SWavePhases_0'] = TCanvas( 'SWavePhases_canv0', '', 2 )
    canvs['SWavePhases_1'] = TCanvas('SWavePhases_canv1')
    canvs['SWavePhases_1'].SetLeftMargin(0.18)
    canvs['SWavePhases_1'].SetRightMargin(0.05)
    canvs['SWavePhases_1'].SetBottomMargin(0.20)
    canvs['SWavePhases_1'].SetTopMargin(0.05)

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
            parVal = float(line[2]) + parSettings[par].get( 'offset', 0. )
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
    NLLValsPara = array( 'd', [ 0.5 * ( ( val - parMean - parSettings[par].get( 'offset', 0. ) ) / parErr )**2 for val in parValsPara ] )

    graphs[par] = (  TGraph( len(parValsPara), parValsPara, NLLValsPara )
                   #, TGraph( len(parVals), parVals, NLLVals )
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
        graph.SetLineStyle( kDashed if it == 0 else kSolid )
        graph.SetLineColor( kBlue if it == 1 else kRed if it == 2 else kBlack )
        graph.SetMarkerColor( kBlue if it == 1 else kRed if it == 2 else kBlack )
        graph.SetMinimum( NLLRange[0] )
        graph.SetMaximum( NLLRange[1] )
        graph.GetXaxis().SetLimits( parMin, parMax )
        graph.GetXaxis().SetTitle(parSettings[par]['name'])
        graph.GetYaxis().SetTitle(NLLLabel)
        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset(0.7)
        graph.GetXaxis().SetNdivisions( 7, 5, 0 )
        graph.Draw( 'AL' if it == 0 else 'SAMES L' )
        #graph.Draw( 'AL' if it == 0 else 'SAMES PC' )
    canvs[par].Print(plotFilePath)

    if par in SWavePhases :
        ind = SWavePhases.index(par)

        padMargins = ( 0.08, 0.02 )
        padOverlap = 0.02
        padHeight = ( 1. - sum(padMargins) + float( len(SWavePhases) - 1 ) * padOverlap ) / float( len(SWavePhases) )
        it = float( len(SWavePhases) - 1 - ind )
        pads[par] = TPad( par, par, 0.2, padMargins[0] + it * ( padHeight - padOverlap )
                                  , 0.8, padMargins[0] + ( it + 1. ) * ( padHeight - padOverlap ) + padOverlap )
        pads[par].SetNumber( ind + 1 )
        pads[par].SetLeftMargin(0.07)
        pads[par].SetRightMargin(0.)
        pads[par].SetBottomMargin(0.20)
        pads[par].SetTopMargin(0.)
        canvs['SWavePhases_0'].cd()
        pads[par].Draw()

        canvs['SWavePhases_0'].cd( ind + 1 )
        graphs0[par] = graphs[par][-1].Clone()
        graphs0[par].SetMinimum( SWaveNLLRange[0] )
        graphs0[par].SetMaximum( SWaveNLLRange[1] )
        graphs0[par].SetLineWidth(2)
        graphs0[par].GetYaxis().SetNdivisions( 2, 5, 0 )
        graphs0[par].GetXaxis().SetTitle('')
        graphs0[par].GetYaxis().SetTitle('')
        graphs0[par].GetXaxis().SetLabelSize(0.18)
        graphs0[par].GetYaxis().SetLabelSize(0.18)
        graphs0[par].GetXaxis().SetLabelOffset( 0.030 if ind == len(SWavePhases) - 1 else 0.1 )
        graphs0[par].GetYaxis().SetLabelOffset(0.010)

        label0 = TLatex()
        label0.SetTextAlign(31)
        label0.SetTextSize(0.18)

        graphs0[par].Draw('AC')
        if ind in SWavePlotLabels :
            label0.DrawLatex( SWavePlotLabels[ind][0], SWavePlotLabels[ind][1], SWavePlotLabels[ind][2] )
        else :
            label0.DrawLatex( 1.6, 4.5, 'bin %d' % ( bin + 1 ) )

        colour = SWaveColours[ind] if ind in SWaveColours else kBlue
        style = SWaveStyles[ind] if ind in SWaveStyles else 0
        graphs1[par] = graphs[par][-1].Clone()
        graphs1[par].SetMinimum( SWaveNLLRange[0] )
        graphs1[par].SetMaximum( SWaveNLLRange[1] )
        graphs1[par].GetXaxis().SetTitle(SWavePhasesLabel)
        graphs1[par].SetLineColor(colour)
        graphs1[par].SetMarkerColor(colour)
        graphs1[par].SetLineStyle(style)
        graphs1[par].GetXaxis().SetNdivisions( 5, 5, 0 )

        label1 = TLatex()
        label1.SetTextAlign(11)
        label1.SetTextSize(0.06)
        label1.SetTextColor(colour)

        canvs['SWavePhases_1'].cd()
        graphs1[par].Draw( 'AC' if firstSWave else 'SAMES C' )
        label1.DrawLatex( SWaveLabelPos[ind][0], SWaveLabelPos[ind][1], 'b%d' % ( ind + 1 ) )
        firstSWave = False

if SWavePhases :
    canvs['SWavePhases_0'].cd()
    xLabel = TLatex()
    xLabel.SetTextAlign(31)
    xLabel.SetTextSize(0.035)
    xLabel.DrawLatex( 0.80, 0.04, SWavePhasesLabel )
    yLabel = TLatex()
    yLabel.SetTextAlign(11)
    yLabel.SetTextAngle(90)
    yLabel.SetTextSize(0.035)
    yLabel.DrawLatex( 0.20, 0.84, NLLLabel )

    canvs['SWavePhases_0'].Print(plotFilePath)
    canvs['SWavePhases_1'].Print(plotFilePath)
canvs['dummy'].Print( plotFilePath + ']' )
