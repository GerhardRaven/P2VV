###########################################################################################################################################
## set script parameters ##
###########################

from math import pi
from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_2011Analysis as PdfConfig
pdfConfig = PdfConfig()

# job parameters
doFit                   = False
pdfConfig['selection']  = 'paper2012'

parFileIn  = ''
parFileOut = ''

pdfConfig['nTupleName']     = 'DecayTree'
pdfConfig['nTupleFile']     = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs2JpsiPhi_ntupleB_for_fitting_20121012_MagDownMagUp.root'
pdfConfig['nominalDataSet'] = True

pdfConfig['timeEffHistFile']      = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Bs_HltPropertimeAcceptance_Data-20120816.root'
pdfConfig['timeEffHistUBName']    = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1DiMuon_Hlt2DiMuonDetached_Reweighted'
pdfConfig['timeEffHistExclBName'] = 'Bs_HltPropertimeAcceptance_PhiMassWindow30MeV_NextBestPVCut_Data_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
pdfConfig['angEffMomentsFile']    = '/project/bfys/jleerdam/data/Bs2Jpsiphi/hel_UB_UT_trueTime_BkgCat050_KK30_Basis_weights'

# fit options
fitOpts = dict(  NumCPU    = 2
               , Optimize  = 2
               , Minimizer = 'Minuit2'
#               , Offset    = True
#               , Hesse     = False
               , Timer     = True
#               , Verbose   = True
              )
pdfConfig['fitOptions'] = fitOpts

fitRange      = ''
corrSFitErr   = 'sumWeight' # '' / 'sumWeight' / ( 0.887, [ 0.566, 0.863, 0.956, 0.948, 0.855, 0.662 ] ) / 'matrix'
randomParVals = ( ) # ( 1., 12345 )
MinosPars     = [#  'AparPhase'
                 #, 'f_S_bin0',        'f_S_bin1',        'f_S_bin2',        'f_S_bin3',        'f_S_bin4',        'f_S_bin5'
                 #, 'ASOddPhase_bin0', 'ASOddPhase_bin1', 'ASOddPhase_bin2', 'ASOddPhase_bin3', 'ASOddPhase_bin4', 'ASOddPhase_bin5'
                ]

# PDF options
pdfConfig['multiplyByTimeEff']    = 'signal'
pdfConfig['timeEffType']          = 'paper2012'
pdfConfig['multiplyByAngEff']     = 'weights'
pdfConfig['parameterizeKKMass']   = 'simultaneous'
pdfConfig['SWeightsType']         = 'simultaneousFreeBkg'
pdfConfig['KKMassBinBounds']      = [ 990., 1020. - 12., 1020. -  4., 1020., 1020. +  4., 1020. + 12., 1050. ]
pdfConfig['SWaveAmplitudeValues'] = (  [ (0.23, 0.08), (0.067, 0.029), (0.008, 0.011), (0.016, 0.011), (0.055, 0.026), (0.17,  0.04) ]
                                     , [ (1.3,  0.7 ), (0.77,  0.28 ), (0.50,  0.47 ), (-0.51, 0.25 ), (-0.46, 0.21 ), (-0.65, 0.20) ] )
pdfConfig['CSPValues']            = [ 0.966, 0.956, 0.926, 0.926, 0.956, 0.966 ]

pdfConfig['sameSideTagging']    = True
pdfConfig['conditionalTagging'] = True
pdfConfig['continuousEstWTag']  = True
pdfConfig['constrainTagging']   = 'constrain'

pdfConfig['timeResType']           = 'eventNoMean'
pdfConfig['numTimeResBins']        = 40
pdfConfig['constrainTimeResScale'] = 'fixed'

pdfConfig['constrainDeltaM'] = 'constrain'

pdfConfig['lambdaCPParam'] = 'lambPhi'

pdfConfig['makePlots'] = False


###########################################################################################################################################
## read data and build PDF ##
#############################

# workspace
from P2VV.RooFitWrappers import RooObject
worksp = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

from P2VV.Parameterizations.FullPDFs import Bs2Jpsiphi_PdfBuilder as PdfBuilder
pdfBuild = PdfBuilder( **pdfConfig )
pdf = pdfBuild.pdf()

if not 'Optimize' in fitOpts or fitOpts['Optimize'] < 2 :
    # unset cache-and-track
    for par in pdfBuild['taggingParams'].parameters() : par.setAttribute( 'CacheAndTrack', False )

if parFileIn :
    # read parameters from file
    pdfConfig.readParametersFromFile( filePath = parFileIn )
    pdfConfig.setParametersInPdf(pdf)

# signal and background data sets
sigData = pdfBuild['sigSWeightData']
bkgData = pdfBuild['bkgSWeightData']
data    = pdfBuild['data']

# data set with weights corrected for background dilution: for phi_s fit only!
if corrSFitErr == 'sumWeight'\
        or ( type(corrSFitErr) != str and hasattr( corrSFitErr, '__iter__' ) and hasattr( corrSFitErr, '__getitem__' ) ) :
    from P2VV.GeneralUtils import correctSWeights
    fitData = correctSWeights( pdfBuild['sigSWeightData'], 'N_bkgMass_sw'
                              , 'KKMassCat' if pdfConfig['parameterizeKKMass'] == 'simultaneous' else ''
                              , CorrectionFactors = None if corrSFitErr == 'sumWeight' else corrSFitErr )

else :
    fitData = pdfBuild['sigSWeightData']

# get observables and parameters in PDF
pdfObs  = pdf.getObservables(fitData)
pdfPars = pdf.getParameters(fitData)


KKMass     = pdfBuild['observables']['KKMass']
mumuMass   = pdfBuild['observables']['mumuMass']


####################################################################################################
## Make the mumu mass plots
##################################################

from P2VV.Load import LHCbStyle

from P2VV.Parameterizations.MassPDFs import CB_Signal_Mass, DoubleCB_Signal_Mass, Linear_Background_Mass
#mumuSig = CB_Signal_Mass(  Name        = 'sig_mumu'
#                         , Prefix      = 'mumu_'
#                         , mass        = mumuMass
#                         , mumu_m_sig_mean  = dict( Value = 3100., MinMax = ( 3050., 3150. ) )
#                         , mumu_m_sig_sigma = dict( Value = 13.4,  MinMax = ( 12.5,  14.5  ) )
#                         , mumu_m_sig_alpha = dict( Value = 1.7,   MinMax = ( 1.0,   2.0   ) )
#                         , mumu_m_sig_n     = dict( Value = 2. )
#                        )
mumuSig = DoubleCB_Signal_Mass(  Name                = 'sig_mumu'
                               , Prefix              = 'mumu_'
                               , mass                = mumuMass
                               , mumu_m_sig_frac     = dict( Value = 0.81,   Error = 0.03, MinMax = ( 0.,    1.    ) )
                               , mumu_m_sig_mean     = dict( Value = 3099.7, Error = 0.1,  MinMax = ( 3050., 3150. ) )
                               , mumu_m_sig_sigma_1  = dict( Value = 12.5,   Error = 0.2,  MinMax = ( 10.,   15.   ) )
                               , mumu_m_sig_sigma_sf = dict( Value = 1.7,    Error = 0.1,  MinMax = ( 0.1,   5.0   ) )
                               , mumu_m_sig_alpha_1  = dict( Value = 1.90,   Error = 0.04, MinMax = ( 1.0,   3.0   ) )
                               , mumu_m_sig_alpha_sf = dict( Value = 1.4,    Error = 0.3,  MinMax = ( 0.1,   5.0   )
                                                            , ObjectType = 'RealVar' )
                               , mumu_m_sig_n_1      = dict( Value = 2. )
                               , mumu_m_sig_n_2      = dict( Value = 2. )
                              )
#mumuBkg = Linear_Background_Mass( Name = 'bkg_mumu', mass = mumuMass, Constant = True )

#Initialize Components
from P2VV.RooFitWrappers import buildPdf, Component
mumuMassSigComp = Component('mumuSig', ( mumuSig.pdf(), ), Yield = ( 27000.,                        0., sigData.sumEntries() ) )
#mumuMassBkgComp = Component('mumuBkg', ( mumuBkg.pdf(), ), Yield = ( sigData.sumEntries() - 27000., 0., sigData.sumEntries() ) )

#Build mumu Mass Pdf
#mumuMassPdf = buildPdf(Components = ( mumuMassSigComp, mumuMassBkgComp ), Observables = (mumuMass, ), Name='mumuMassPdf'  )
mumuMassPdf = buildPdf(Components = ( mumuMassSigComp, ), Observables = (mumuMass, ), Name='mumuMassPdf'  )

#Fit
#ws.var('m_bkg_arg').setRange(-1e-4,.5) 
fitResult = mumuMassPdf.fitTo( sigData, SumW2Error = False, Save = True, **fitOpts )
fitResult.PrintSpecial( text = True )

#Plot
from ROOT import kFullCircle, kRed, kGreen
mumuMassPlot = mumuMass.frame(60)

#sigArgs = {'Components':'sig_mumu', 'LineColor':kRed,     'LineStyle':10, 'LineWidth':3}
#bkgArgs = {'Components':'bkg_mumu', 'LineColor':kGreen+3, 'LineStyle': 2, 'LineWidth':3}

sigData.plotOn( mumuMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.7, LineWidth = 3)
mumuMassPdf.plotOn(mumuMassPlot, LineWidth = 3)
#mumuMassPdf.plotOn(mumuMassPlot, **sigArgs)
#mumuMassPdf.plotOn(mumuMassPlot, **bkgArgs)

binWidth = ( mumuMassPlot.GetXaxis().GetXmax() - mumuMassPlot.GetXaxis().GetXmin() ) / float( mumuMassPlot.GetNbinsX() )
mumuMassPlot.SetXTitle('m(#mu^{+}#mu^{-}) [MeV/c^{2}]')
mumuMassPlot.SetYTitle('Candidates / (%d MeV/c^{2})' % binWidth )
mumuMassPlot.SetTitleOffset( 1.1,  'x' )
mumuMassPlot.SetTitleOffset( 1.15, 'y' )

#LHCbLabel
from ROOT import TPaveText
LHCbLabel = TPaveText( 0.24, 0.81, 0.37, 0.89, 'BRNDC')
LHCbLabel.AddText('LHCb')
LHCbLabel.SetFillColor(0)
LHCbLabel.SetTextAlign(12)
LHCbLabel.SetTextSize(0.072)
LHCbLabel.SetBorderSize(0)

from ROOT import TCanvas
mumuMassCanv = TCanvas('mumuMassCanv')
mumuMassCanv.SetLeftMargin(0.18)
mumuMassCanv.SetRightMargin(0.05)
mumuMassCanv.SetBottomMargin(0.18)
mumuMassCanv.SetTopMargin(0.05)
mumuMassPlot.Draw()
LHCbLabel.Draw()
mumuMassCanv.Print('mumuMass.ps(')

mumuMassPlot.SetMinimum(3.)
mumuMassPlot.SetMaximum(3.e3)
mumuMassPlot.SetTitleOffset( 1.0, 'y' )

mumuMassCanvLog = TCanvas('mumuMassCanvLog')
mumuMassCanvLog.SetLeftMargin(0.18)
mumuMassCanvLog.SetRightMargin(0.05)
mumuMassCanvLog.SetBottomMargin(0.18)
mumuMassCanvLog.SetTopMargin(0.05)
mumuMassCanvLog.SetLogy(True)
mumuMassPlot.Draw()
LHCbLabel.Draw()
mumuMassCanvLog.Print('mumuMass.ps)')


##################################################
## Make the KK mass plots
##################################################
#Build KK mas pdf
from ROOT import RooRealVar, RooRelBreitWigner, RooConstVar, RooFFTConvPdf, RooGaussModel

#Mass Resolution Model
KKMassVar = KKMass._target_()
resMean  = RooRealVar( 'resMean'   , 'resMean'   ,  0.62, -5.,  5.  ) 
resSigma = RooRealVar( 'resSigma'  , 'resSigma'  ,  1.21, 0.01, 5.  )
resMean.setError(0.02)
resSigma.setError(0.04)

GaussModel = RooGaussModel( 'GaussModel', 'GaussModel', KKMassVar, resMean, resSigma )

#Build Phi Mass Pdf
mean   = RooConstVar( 'mean'  , 'mean'  , 1019.455 )
width  = RooConstVar( 'width' , 'width' , 4.26    )
spin   = RooConstVar( 'spin'  , 'spin'  ,  1.      )  
radius = RooConstVar( 'radius', 'radius', -1.      )
K1mass = RooConstVar( 'K1mass', 'K1mass', 493.68   )
K2mass = RooConstVar( 'K2mass', 'K2mass', 493.68   )

phiMassBW  = RooRelBreitWigner('phiMassBW',  'phiMassBW',  KKMassVar, mean, width, spin, radius, K1mass, K2mass )
phiMassPdf = RooFFTConvPdf    ('phiMassPdf', 'phiMassPdf', KKMassVar, phiMassBW, GaussModel )

#s-Wave Pdf
from ROOT import RooDstD0BG as PhaseSpaceFunc
dm0   = RooConstVar( 'dm0'  , 'dm0', 987.36 )#, 980., 1010. )
parA = RooConstVar(  'parA', 'parA', 1.     )#,  0., 20. )
parB = RooRealVar(  'parB', 'parB', 11.,  -20., 40. )
parC = RooRealVar(   'parC', 'parC', 1.4, 0., 10. )
parB.setError(5.)
parC.setError(1.)

KKBkgPHSP = PhaseSpaceFunc( 'KKBkgPHSP', 'KKBkgPHSP', KKMassVar, dm0, parC, parA, parB )
KKBkgPdf  = RooFFTConvPdf ( 'KKBkgPdf' , 'KKBkgPdf',  KKMassVar, KKBkgPHSP, GaussModel )

#Total Pdf
from ROOT import RooAddPdf
phiFrac = RooRealVar( 'phiFrac', 'phiFrac', 0.96, 0., 1. )
phiFrac.setError(0.003)
KKMassPdf = RooAddPdf( 'KKMassPdf', 'KKMassPdf', phiMassPdf, KKBkgPdf, phiFrac )

#Fit
fitResult = KKMassPdf.fitTo( sigData, SumW2Error = False, Save = True, **fitOpts )
fitResult.PrintSpecial( text = True )

#Plot
KKMassPlot = KKMassVar.frame(60)

#comps = { 'PhiMassPdf':dict( LineColor=kRed,      LineStyle=10, LineWidth=lineWidth ), 
#          'KKbkgPdf'  :dict( LineColor=kGreen +3, LineStyle= 2, LineWidth=lineWidth )
#        }
sigData.plotOn( KKMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.7, LineWidth = 3 )
KKMassPdf.plotOn( KKMassPlot, LineWidth = 3 )

binWidth = ( KKMassPlot.GetXaxis().GetXmax() - KKMassPlot.GetXaxis().GetXmin() ) / float( KKMassPlot.GetNbinsX() )
KKMassPlot.SetXTitle('m(K^{+}K^{-}) [MeV/c^{2}]')
KKMassPlot.SetYTitle('Candidates / (%d MeV/c^{2})' % binWidth )
KKMassPlot.SetTitleOffset( 1.1,  'x' )
KKMassPlot.SetTitleOffset( 1.15, 'y' )

KKMassPlot.SetMinimum(0.)

KKMassCanv = TCanvas('KKMassCanv')
KKMassCanv.SetLeftMargin(0.18)
KKMassCanv.SetRightMargin(0.05)
KKMassCanv.SetBottomMargin(0.18)
KKMassCanv.SetTopMargin(0.05)
KKMassPlot.Draw()
LHCbLabel.Draw()
KKMassCanv.Print('KKMass.ps(')

KKMassPlot.SetMinimum(2.)
KKMassPlot.SetMaximum(1.e4)
KKMassPlot.SetTitleOffset( 1.0, 'y' )

KKMassCanvLog = TCanvas('KKMassCanvLog')
KKMassCanvLog.SetLeftMargin(0.18)
KKMassCanvLog.SetRightMargin(0.05)
KKMassCanvLog.SetBottomMargin(0.18)
KKMassCanvLog.SetTopMargin(0.05)
KKMassCanvLog.SetLogy(True)
KKMassPlot.Draw()
LHCbLabel.Draw()
KKMassCanvLog.Print('KKMass.ps)')
