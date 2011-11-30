###########################################################################################################################################
## P2VVParameterizations.TimeResolution: Time resolution models                                                                          ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##                                                                                                                                       ##
###########################################################################################################################################

class ResolutionModelLP2011 :
    def __init__( self, t ) :
        from RooFitWrappers import RealVar, ConstVar, ResolutionModel, AddModel
        mu = ConstVar('tres_mu', Value = -0.0027 )
        SF = RealVar('tres_SF', Value = 1.0, MinMax = (0.5,5) )
        from ROOT import RooGaussModel as GaussModel
        sigmas = [ (3,0.513 ), (2,0.0853), (1,0.0434) ]
        frac   = [ (3,0.0017), (2,0.165) ]
        self.Model = AddModel('tres', [ ResolutionModel('tres_%s'%n, Type = GaussModel, Observables = [ t ], Parameters = [ mu, ConstVar('tres_s%s'%n, Value = v  ), SF ] ) for n,v in sigmas ]
                                    , [ ConstVar('tres_f%s'%n,Value = v) for n,v in frac ]
                             )

