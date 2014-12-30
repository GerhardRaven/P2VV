###########################################################################################################################################
## Script options ##
####################

dataSetFile = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/P2VVDataSets20112012Reco14_I2Mass_6KKMassBins_2TagCats_20140619.root'

mumuPlotsFilePath = 'mumuMass.pdf'
KKPlotsFilePath   = 'KKMass.pdf'

drawTotalDists = True
fitOpts = dict( NumCPU = 8, Optimize = 2, Timer = True, Minimizer = 'Minuit2', Strategy = 1, Offset = True )

from ROOT import gStyle, kFullCircle, kSolid, kBlue, kRed, kMagenta
gStyle.SetLineStyleString( 5, ' 40 20 10 20'  )
gStyle.SetLineStyleString( 7, ' 40 20'        )
gStyle.SetLineStyleString( 9, ' 100 20'       )


###########################################################################################################################################
## Read data and get observables ##
###################################

from P2VV.Load import P2VVLibrary, LHCbStyle
from ROOT import gStyle
gStyle.SetColorModelPS(1)

# workspace
from P2VV.RooFitWrappers import RooObject
ws = RooObject( workspace = 'JpsiphiWorkspace' ).ws()

# read data set from file
from P2VV.Utilities.DataHandling import readData
if drawTotalDists:
    data = readData( filePath = dataSetFile, dataSetName = 'JpsiKK', NTuple = False )
    data = data.reduce('hlt2_biased==1 && mass>5338. && mass<5398.')
    data.Print()
sigData = readData( filePath = dataSetFile, dataSetName = 'JpsiKK_sigSWeight', NTuple = False )
sigData = sigData.reduce('hlt2_biased==1')
sigData.Print()

# get observables
from P2VV.RooFitWrappers import RealVar
mumuMass = RealVar('mdau1')
KKMass   = RealVar('mdau2')

# LHCb label
labelText = '' # 'LHCb'
if labelText :
    from ROOT import TLatex
    label = TLatex()
    label.SetTextAlign(12)
    label.SetTextSize(0.072)


###########################################################################################################################################
## Make mumu mass plots ##
##########################

from P2VV.Parameterizations.MassPDFs import CB_Signal_Mass, DoubleCB_Signal_Mass, Linear_Background_Mass
#mumuSig = CB_Signal_Mass(  Name             = 'sig_mumu'
#                         , ParNamePrefix    = 'mumu'
#                         , mass             = mumuMass
#                         , mumu_m_sig_mean  = dict( Value = 3100., MinMax = ( 3050., 3150. ) )
#                         , mumu_m_sig_sigma = dict( Value = 13.4,  MinMax = ( 12.5,  14.5  ) )
#                         , mumu_m_sig_alpha = dict( Value = 1.7,   MinMax = ( 1.0,   2.0   ) )
#                         , mumu_m_sig_n     = dict( Value = 2. )
#                        )
mumuSig = DoubleCB_Signal_Mass(  Name           = 'sig_mumu'
                               , ParNamePrefix  = 'mumu'
                               , mass           = mumuMass
                               , m_sig_frac     = dict( Value = 0.81,   Error = 0.03, MinMax = ( 0.,    1.    ) )
                               , m_sig_mean     = dict( Value = 3099.7, Error = 0.1,  MinMax = ( 3050., 3150. ) )
                               , m_sig_sigma_1  = dict( Value = 12.5,   Error = 0.2,  MinMax = ( 8.,    15.   ) )
                               , m_sig_sigma_sf = dict( Value = 1.7,    Error = 0.1,  MinMax = ( 0.1,   5.0   ) )
                               , m_sig_alpha_1  = dict( Value = 1.90,   Error = 0.04, MinMax = ( 1.0,   3.0   ) )
                               , m_sig_alpha_sf = dict( Value = 1.4,    Error = 0.3,  MinMax = ( 0.1,   5.0   ), ObjectType = 'RealVar' )
                               , m_sig_n_1      = dict( Value = 2. )
                               , m_sig_n_2      = dict( Value = 2. )
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
mumuMassPlot = mumuMass.frame(60)

#sigArgs = {'Components':'sig_mumu', 'LineColor':kRed,     'LineStyle':10, 'LineWidth':3}
#bkgArgs = {'Components':'bkg_mumu', 'LineColor':kGreen+3, 'LineStyle': 2, 'LineWidth':3}

if drawTotalDists :
    data.plotOn( mumuMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.5, MarkerColor = kRed, LineWidth = 2, LineColor = kRed )
sigData.plotOn( mumuMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.5, LineWidth = 2)
mumuMassPdf.plotOn(mumuMassPlot, LineWidth = 3)
#mumuMassPdf.plotOn(mumuMassPlot, **sigArgs)
#mumuMassPdf.plotOn(mumuMassPlot, **bkgArgs)

binWidth = ( mumuMassPlot.GetXaxis().GetXmax() - mumuMassPlot.GetXaxis().GetXmin() ) / float( mumuMassPlot.GetNbinsX() )
mumuMassPlot.SetXTitle('m(#mu^{+}#mu^{-}) [MeV/c^{2}]')
mumuMassPlot.SetYTitle('Candidates / (%.2g MeV/c^{2})' % binWidth )
mumuMassPlot.SetMinimum(0.)
mumuMassPlot.SetMaximum( 6200. if not drawTotalDists else 8000. )
mumuMassPlot.SetTitleOffset( 1.10, 'x' )
mumuMassPlot.SetTitleOffset( 1.15, 'y' )

from ROOT import TCanvas
mumuMassCanv = TCanvas('mumuMassCanv')
mumuMassCanv.SetLeftMargin(0.18)
mumuMassCanv.SetRightMargin(0.05)
mumuMassCanv.SetBottomMargin(0.18)
mumuMassCanv.SetTopMargin(0.05)
mumuMassPlot.Draw()
if labelText : label.DrawLatexNDC( 0.25, 0.85, labelText )
mumuMassCanv.Print( mumuPlotsFilePath + '(' )

mumuMassPlot.SetMinimum(2.e1)
mumuMassPlot.SetMaximum(1.e4)
mumuMassPlot.SetTitleOffset( 1.0, 'y' )

mumuMassCanvLog = TCanvas('mumuMassCanvLog')
mumuMassCanvLog.SetLeftMargin(0.18)
mumuMassCanvLog.SetRightMargin(0.05)
mumuMassCanvLog.SetBottomMargin(0.18)
mumuMassCanvLog.SetTopMargin(0.05)
mumuMassCanvLog.SetLogy(True)
mumuMassPlot.Draw()
if labelText : label.DrawLatexNDC( 0.25, 0.85, labelText )
mumuMassCanvLog.Print( mumuPlotsFilePath + ')' )


###########################################################################################################################################
## Make KK mass plots ##
########################

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

KKSWavePHSP = PhaseSpaceFunc( 'KKSWavePHSP', 'KKSWavePHSP', KKMassVar, dm0, parC, parA, parB )
KKSWavePdf  = RooFFTConvPdf ( 'KKSWavePdf' , 'KKSWavePdf',  KKMassVar, KKSWavePHSP, GaussModel )

#Total Pdf
from ROOT import RooAddPdf
phiFrac = RooRealVar( 'phiFrac', 'phiFrac', 0.96, 0., 1. )
phiFrac.setError(0.003)
KKMassPdf = RooAddPdf( 'KKMassPdf', 'KKMassPdf', phiMassPdf, KKSWavePdf, phiFrac )

#Fit
fitResult = KKMassPdf.fitTo( sigData, SumW2Error = False, Save = True, Minos = False, **fitOpts )
fitResult.PrintSpecial( text = True )

#Plot
KKMassPlot = KKMassVar.frame(120)

if drawTotalDists :
    data.plotOn( KKMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.5, MarkerColor = kRed, LineWidth = 2, LineColor = kRed )
sigData.plotOn( KKMassPlot, MarkerStyle = kFullCircle, MarkerSize = 0.5, LineWidth = 2 )
KKMassPdf.plotOn( KKMassPlot, LineStyle = kSolid, LineWidth = 3, LineColor = kBlue                                   )
#KKMassPdf.plotOn( KKMassPlot, LineStyle = 7,      LineWidth = 3, LineColor = kRed,         Components = 'phiMassPdf' )
#KKMassPdf.plotOn( KKMassPlot, LineStyle = 5,      LineWidth = 3, LineColor = kMagenta + 3, Components = 'KKSWavePdf' )

binWidth = ( KKMassPlot.GetXaxis().GetXmax() - KKMassPlot.GetXaxis().GetXmin() ) / float( KKMassPlot.GetNbinsX() )
KKMassPlot.SetXTitle('m(K^{+}K^{-}) [MeV/c^{2}]')
KKMassPlot.SetYTitle('Candidates / (%.2g MeV/c^{2})' % binWidth )
KKMassPlot.SetMinimum(0.)
KKMassPlot.SetMaximum( 6200. if not drawTotalDists else 7000. )
KKMassPlot.SetTitleOffset( 1.10, 'x' )
KKMassPlot.SetTitleOffset( 1.15, 'y' )

KKMassPlot.SetMinimum(0.)

KKMassCanv = TCanvas('KKMassCanv')
KKMassCanv.SetLeftMargin(0.18)
KKMassCanv.SetRightMargin(0.05)
KKMassCanv.SetBottomMargin(0.18)
KKMassCanv.SetTopMargin(0.05)
KKMassPlot.Draw()
if labelText : label.DrawLatexNDC( 0.25, 0.85, labelText )
KKMassCanv.Print( KKPlotsFilePath + '(' )

KKMassPlot.SetMinimum(1.e1)
KKMassPlot.SetMaximum(1.e4)
KKMassPlot.SetTitleOffset( 1.0, 'y' )

KKMassCanvLog = TCanvas('KKMassCanvLog')
KKMassCanvLog.SetLeftMargin(0.18)
KKMassCanvLog.SetRightMargin(0.05)
KKMassCanvLog.SetBottomMargin(0.18)
KKMassCanvLog.SetTopMargin(0.05)
KKMassCanvLog.SetLogy(True)
KKMassPlot.Draw()
if labelText : label.DrawLatexNDC( 0.25, 0.85, labelText )
KKMassCanvLog.Print( KKPlotsFilePath + ')' )
