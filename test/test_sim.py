
class TestModule( object ):

    def test_import( self ):
        """ Test if importing the module works and it has the classes we need. """
        import RooFitWrappers as wrappers
        assert hasattr( wrappers, 'RooObject' )
        assert hasattr( wrappers, 'Category' )
        assert hasattr( wrappers, 'FormulaVar' )
        assert hasattr( wrappers, 'RealVar' )
        assert hasattr( wrappers, 'Pdf' )
        assert hasattr( wrappers, 'ProdPdf' )
        assert hasattr( wrappers, 'SumPdf' )
        assert hasattr( wrappers, 'ResolutionModel' )
        assert hasattr( wrappers, 'Component' )
        assert hasattr( wrappers, 'buildPdf' )

class TestPDFs( object ):

    def setup_class( cls ):
        # setup (singleton) workspace
        from ROOT import RooWorkspace
        from RooFitWrappers import RooObject
        ws = RooObject()
        ws.setWorkspace(RooWorkspace("test_workspace"))

    def test_variables( self ):
        # now (consistently!) create/declare observables
        from RooFitWrappers import RealVar
        m = RealVar('m',Observable=True,Unit='MeV/c^2',MinMax=(5000,6000))
        assert m[ 'Observable' ] == True
        assert m[ 'Unit' ] == 'MeV/c^2'
        assert m[ 'MinMax' ] == (5000,6000)

        res_mean = RealVar( 'res_mean', Observable = False, Unit = 'ps', Value = 0,
                            MinMax = ( -10, 10 ) )
        assert res_mean[ 'Observable' ] == False
        assert res_mean[ 'Value' ] == 0

        c = RealVar('c',Observable=False,Unit='ps',Value=0)
        assert c[ 'MinMax' ] == (-1e30, 1e30)

    def test_resolution( self ):
        from RooFitWrappers import RealVar
        from RooFitWrappers import ResolutionModel

        t = RealVar('t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)

        mean = RealVar( 'mean', Observable = False, Unit = 'ps', Value = 0,
                        MinMax = ( -10, 10 ) )
        sigma = RealVar( 'sigma', Observable = False, Unit = '1/ps', Value = 50,
                         MinMax = (  20, 60 ) )

        res = ResolutionModel( 'res', Type = 'RooGaussModel', Observables = [ t ],
                               Parameters = [ mean, sigma ] )
        assert res[ 'Type' ] == 'RooGaussModel'
        assert res[ 'Observables' ] == frozenset( [ t ] )
        assert res[ 'Parameters' ] == frozenset( [ mean, sigma ] )

    def test_PDF( self ):
        from RooFitWrappers import RealVar
        from RooFitWrappers import ResolutionModel
        from RooFitWrappers import Pdf

        t = RealVar('t',Observable=True,MinMax=(-1,14),Unit='ps',Value=0)
        assert t[ 'Observable' ] == True
        assert t[ 'Unit' ] == 'ps'
        assert t[ 'MinMax' ] == (-1,14)

        mean = RealVar( 'mean', Observable = False, Unit = 'ps', Value = 0,
                        MinMax = ( -10, 10 ) )
        sigma = RealVar( 'sigma', Observable = False, Unit = '1/ps', Value = 50,
                         MinMax = (  20, 60 ) )

        res = ResolutionModel( 'res', Type = 'RooGaussModel', Observables = [ t ],
                               Parameters = [ mean, sigma ] )
        assert res[ 'Type' ] == 'RooGaussModel'
        assert res[ 'Observables' ] == frozenset( [ t ] )
        assert res[ 'Parameters' ] == frozenset( [ mean, sigma ] )

        tau = RealVar( 'tau', Observable = False, Unit = 'ps', Value = 1.5, MinMax = ( 1., 2. ) )
        sig = Pdf( 'time', Type = 'Decay', Observables = ( t, ), Parameters = ( tau, ),
                   Options = ( 'SingleSided', ), ResolutionModel = res )

        assert sig[ 'Type' ] == 'Decay'
        assert sig[ 'Observables' ] == frozenset([t])
        assert sig[ 'Parameters' ] == frozenset([tau])

## mass_mean = RealVar( 'mass_mean', Observable = False, Unit = 'MeV', Value = 5300, MinMax = ( 5200, 5800 ) )
## mass_sigma = RealVar( 'mass_sigma', Observable = False, Unit = 'MeV', Value = 15, MinMax = ( 10, 20 ) )

## sig_m = Pdf( 'mass', Type = 'Gaussian', Observables = ( m, ), Parameters = ( mass_mean, mass_sigma ) )

## # create signal and background
## signal = Component('signal')
## signal.setYield(100,50,150)
## ## signal[m] = 'Gaussian(m,5300,15)'
## signal[m] = sig_m
## ## signal[t] = 'Decay(t,sig_tau[1.5,1.0,2.0],TruthModel(t),SingleSided)'
## signal[t] = sig_t

## background = Component('background')
## background.setYield(1000,900,1100)

## background_c = RealVar( 'background_c', Observable = False, Unit = '1/MeV', Value = -0.0004)
## background[ m ] = Pdf( 'background', Observables = ( m, ), Type = 'Exponential',
##                        Parameters = ( background_c, ) )

## background_tau = RealVar( 'background_tau', Observable = False, Unit = 'ps', Value = 0.4, MinMax = ( 0.1, 0.9 ) )
## background_res = ResolutionModel( 'background_res', Type = 'RooTruthModel', Observables = [ t ] )

## #background[t] = 'Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided)'
## background[ t ] = Pdf( 'background', Type = 'Decay', Observables = ( t, ),
##                        Parameters = ( background_tau, ), Options = ( 'SingleSided', ),
##                        ResolutionModel = background_res )
## # Exponential(m,-0.004)
## #background[m,t] = 'PROD(Exponential(m,-0.004),Decay(t,bkg_tau[0.4,0.1,0.9],TruthModel(t),SingleSided))'

## pdf = buildPdf( (background,signal) , observables = (m,t), name='pdf' )

## ##########################################

## data = pdf.generate((m,t))
## pdf.fitTo(data)

## from ROOT import TCanvas, RooCmdArg, RooFit, kDashed
## sigcolor = RooCmdArg( RooFit.LineColor(RooFit.kGreen ) )
## bkgcolor = RooCmdArg( RooFit.LineColor(RooFit.kRed))
## lw = RooCmdArg(RooFit.LineWidth(2))
## ms = RooCmdArg(RooFit.MarkerSize(0.4))
## xe = RooCmdArg(RooFit.XErrorSize(0))
## dashed = RooCmdArg(RooFit.LineStyle(kDashed))
## c = TCanvas()
## c.Divide(1,2)
## for (cc,obs) in enumerate((m,t)) :
##     plot( c.cd(1+cc),obs,data,pdf,{ 'signal*': (sigcolor,dashed), 'background*': (bkgcolor,dashed)}, logy = (cc==1), dataOpts = (xe,ms) )

## # create a continuous tagging variable...

## # split in categories

