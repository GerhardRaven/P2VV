from RooFitWrappers import Pdf

class TestPdf( Pdf ):
    def __init__( self, name, Observables, ResolutionModel = None ):
        self._init( name, 'RooAbsPdf' )
        self.ws
