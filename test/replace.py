from ROOT import TCanvas, RooCmdArg, RooFit, kDashed
sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
lw = RooCmdArg(RooFit.LineWidth(2))
ms = RooCmdArg(RooFit.MarkerSize(0.4))
xe = RooCmdArg(RooFit.XErrorSize(0))
dashed = RooCmdArg(RooFit.LineStyle(kDashed))


from RooFitWrappers import *

# setup (singleton) workspace
from ROOT import RooWorkspace
ws = RooObject()
ws.setWorkspace( RooWorkspace("myworkspace") )

# now (consistently!) create/declare observables
m = RealVar('m',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
t = RealVar('t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
#c = Category('tagdecision',{+1: 'B', -1 : 'Bbar'}, Value='B',  Observable=True )
#tau = RealVar('sig_tau',Observable=False,Blinded=('UnblindUniform','blindingString', 1.),MinMax=(1,2) )


# create signal and background
signal = Component('signal')
signal.setYield(100,50,150)
signal[m] = 'Gaussian(m,5300,15)'
signal[t] = 'Decay(t,sig_tau[1.5,1.0,2.0],TruthModel(t),SingleSided)'

background = Component('background')
background.setYield(1000,900,1100)
#background[m] = 'Exponential(m,-0.0004)' 
#background[t] = 'Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided)'
background[m,t] = 'PROD(Exponential(m,-0.004),Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided))'


pdf = buildPdf( (background,signal) , observables = (m,t), name='pdf' )

##########################################

data = pdf.generate((m,t))
pdf.fitTo(data)

def plot_pdf( pdf, data ):
    name = pdf.GetName() + '_canvas'
    c = TCanvas( name, name, 500, 1000 )
    c.Divide(1,2)
    for (cc,obs) in enumerate((m,t)) :
        plot( c.cd(1+cc),obs,data,pdf,{ 'signal*': (sigcolor,dashed), 'background*': (bkgcolor,dashed)}, logy = (cc==1), dataOpts = (xe,ms) )

plot_pdf( pdf, data )

tmp = RooObject()._declare( 'GaussModel::res_gauss( t, mean[ 0., -5, 5.], sigma[ 0.05, -0.2, 0.2], mean_sf[ 1. ], sigma_sf[ 1. ] )' )
gauss = Pdf( 'res_gauss' )

from ROOT import RooCustomizer, RooTruthModel
customizer = RooCustomizer(pdf._target_(), 'new_pdf' )
for c in pdf.getComponents() :
    if type(c) is not RooTruthModel : continue
    customizer.replaceArg( c, gauss )
pdf_tmp = customizer.build(True)
pdf_new = 
data_new = pdf_new.generate((m,t))
pdf_new.fitTo( data_new )

plot_pdf( pdf_new, data_new )

# create a continuous tagging variable...

# split in categories

