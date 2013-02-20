def legendre(x) :
        return (1,x,0.5*(3*x*x-1),0.5*(5*x*x*x-3*x) )
def spharmonic(ct,phi) :
        from math import sin,cos,pi,sqrt
        st = sqrt(1-ct*ct)
        x=st*cos(phi)
        y=st*sin(phi)
        z=ct
        n0 = sqrt(1./(4*pi))
        n1 = sqrt(3./(4*pi))
        return ( (n0,), (-n1*y,n1*z,-n1*x),  )  # 0,1,2  -> -1,0,+1

class efficiency :
    def __init__(self,*args) :
        (self.cpsi,self.ctheta,self.phi) = args[0] if len(args)==1 else args
    def accept(self) :
        p = legendre( self.cpsi.getVal() )
        y = spharmonic( self.ctheta.getVal(), self.phi.getVal() )
        from random import random
        assert ( p[0]+0.4*p[2] )/1.4  < 1.0001
        return random() < ( p[0]+0.4*p[2] )/1.4 
        #return random() < ( p[0]+0.2*p[1]+0.4*p[2] )/3 * ( y[0][0] + 0.1*y[1][-1+1] + 0.2*y[1][0+1] + 0.3*y[1][1+1] )


from P2VV.RooFitWrappers import *

ws = RooObject( workspace = 'myws' )

from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
t         = RealVar(  't', Title = 'decay time', Unit='ps',                  Observable = True,  MinMax=(0,14)  )
iTag      = Category( 'tagdecision' , Title = 'initial state flavour tag',   Observable = True,  States = { 'B': +1, 'Bbar': -1 } ) # , 'untagged' : 0 } )
eta       = RealVar(   'eta', Title = 'estimated mis tag', Observable = False, Constant = True, Value = 0.3, MinMax=(0,0.5) )
observables = [ i for i in angles.angles.itervalues() ] + [ t,iTag ]

for i in angles.angles.itervalues() : i.setBins(16)
t.setBins(48)

from P2VV.Parameterizations.CPVParams import LambdaSqArg_CPParam
CP = LambdaSqArg_CPParam( phiCP = { 'Name': 'HelloWorld', 'Value': -0.04, 'MinMax': (-3.2,3.2) }, lambdaCPSq = ConstVar(Name = 'one',Value=1) )

# polar^2,phase transversity amplitudes, with Apar^2 = 1 - Aperp^2 - A0^2, and delta0 = 0
from P2VV.Parameterizations.DecayAmplitudes import JpsiVPolar_AmplitudeSet
amplitudes = JpsiVPolar_AmplitudeSet( A0Mag2 = 0.60, A0Phase = 0
                                    , AperpMag2 = 0.160, AperpPhase = -0.17
                                    , AparPhase = 2.5
                                    , ASMag2 = { 'Value' : 0, 'Constant': True} , ASPhase = { 'Value': 0, 'Constant':True } )
#amplitudes.setConstant('.*AS.*',True)

from P2VV.Parameterizations.TimePDFs import JpsiphiBTagDecayBasisCoefficients
# need to specify what, and in which order, to traverse...
basisCoefficients = JpsiphiBTagDecayBasisCoefficients( angles.functions, amplitudes,CP, ['A0','Apar','Aperp','AS'] ) 

#from P2VV.Parameterizations.TimePDFs import JpsiphiBDecayBasisCoefficients
#basisCoefficients = JpsiphiBDecayBasisCoefficients( angles.functions, amplitudes,CP, iTag,  ['A0','Apar','Aperp','AS'] ) 

from P2VV.Parameterizations.FlavourTagging import WTag_TaggingParams
taggingParams = WTag_TaggingParams( wTag = eta ) # FormulaVar('wTag','@2 + @3*(@0-@1)',[eta,etaAverage,p0,p1] ) )

# now build the actual signal PDF...
from P2VV.Parameterizations.LifetimeParams import Gamma_LifetimeParams
lifetimeParams = Gamma_LifetimeParams( Gamma = 0.68, dGamma = 0.05, dM = dict( Value = 17.8, MinMax = (16,19), Constant = True) )

from P2VV.Parameterizations.TimeResolution import Truth_TimeResolution
args = { 'time'      : t
       , 'iTag'      : iTag
       , 'dm'        : lifetimeParams['dM']
       , 'tau'       : lifetimeParams['MeanLifetime']
       , 'dGamma'    : lifetimeParams['dGamma']
       , 'resolutionModel' : Truth_TimeResolution(time = t)['model']
       , 'coshCoef'  : basisCoefficients['cosh']
       , 'cosCoef'   : basisCoefficients['cos']
       , 'sinhCoef'  : basisCoefficients['sinh']
       , 'sinCoef'   : basisCoefficients['sin']
       , 'avgCEven'  : taggingParams['avgCEven'] 
       , 'avgCOdd'   : taggingParams['avgCOdd']
       , 'dilution'  : taggingParams['dilution']
       , 'ADilWTag'  : taggingParams['ADilWTag']
       }

# TODO: should be able to write BTagDecay('mypdf', **lifetimeParams.BTagDecay() + **basisCoefficients.BTagDecay() + **taggingParams.BTagDecay() )
mcpdf = BTagDecay( 'mc_pdf', **args  )


data = mcpdf.generate(observables,10000)
print 'got dataset with %s entries' % data.numEntries()
allObs = mcpdf.getObservables( data.get() )
allObs.Print()
print mcpdf.Observables()

#assert False
eps = efficiency( angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'] )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
origdata = data
data = inEffData;


print 'computing efficiency moments'
##################################
from P2VV.GeneralUtils import RealMomentsBuilder
# eff = RealMomentsBuilder( Moments = ( RealEffMoment( i, 1, mcpdf, angles.angles.itervalues() ) for v in angles.functions.itervalues() for i in v if i ) )
eff = RealMomentsBuilder()
indices  = [ ( i, l, m ) for i in range(3)
                         for l in range(1)
                         for m in range(-l,l+1) ]
indices += [ ( i, 2, m ) for i in range(3,3)  # 3,20)
                         for m in (-2,1) ] # these are for the 'infinite' series in the signal PDF
eff.appendPYList( angles.angles, indices, PDF = mcpdf, NormSet = allObs)
#eff.compute(data)

#from math import sqrt,pi
#eff.Print( MinSignificance = 0., Names = '.*_ang_.*',   Scales = ( 1. / (16*sqrt(pi)), 1. / (16*sqrt(pi)), 1. ) )
#eff.Print( MinSignificance = 0., Names = 'p2vvab.*',    Scales = ( 1. / ( 2*sqrt(pi)), 1. / ( 2*sqrt(pi)), 1. ) )
#eff.Print( MinSignificance = 0., Names = 'p2vvab.*' )
##################################
class abasis :
    def __init__(self,w,*args) :
       self.w = w
       (cpsi,cheta,phi) = args if len(args)==3 else args[0]

       def _f(x) :
            if type(x) is str : return w[x]
            if not self.w.function(x.GetName()) : x = w.put(x)
            return x
       self.cpsi   = _f(cpsi)
       self.ctheta = _f(cheta)
       self.phi    = _f(phi)
       print 'using %s,%s,%s' % (self.cpsi,self.ctheta,self.phi)

    def angles(self) : 
        return RooArgList(self.cpsi,self.ctheta,self.phi)
    def build(self,label,i,j,k,l,c) :
        name = "%s_%d_%d_%d_%d" % (label,i,j,k,l)
        name.replace("-","m")
        b = self.w.function(name) # workaround a bug in ROOT 5.26 -- if name not present, w.obj(name) will SEGV...
        from ROOT import RooP2VVAngleBasis
        if not b : 
            b = self.w.put( RooP2VVAngleBasis(name,name,self.cpsi,self.ctheta,self.phi,i,j,k,l,c) )
        return b

ab = abasis(mcpdf.ws(),angles.angles['cpsi']._var,angles.angles['ctheta']._var,angles.angles['phi']._var)
# if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
# Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
from ROOT import EffMoment
moments = [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, mcpdf._var, allObs ) for (i,l,m)  in indices ]

# loop over all data, determine moments
def computeMoments( data, moments ) :
    if not moments : return None
    from ROOT import std
    vecmom = std.vector('IMoment*')()
    for m in moments : vecmom.push_back(m)
    from ROOT import _computeMoments
    return _computeMoments( data, vecmom )

computeMoments(data,moments)
for i in moments : i.Print()

##################################
eff.compute(data)
eff.Print( MinSignificance = 0., Names = 'p2vvab.*' )


#pdf.Print("T")
pdf = eff * mcpdf
print pdf.GetName(), type(pdf)
#pdf2.Print("T")

# create generic PDF to describe the angular distribution
#moms = RealMomentsBuilder()
#moms.appendPYList( angles.angles, indices )
#moms.compute(data)
#mom_pdf = moms.createPDF( Name = 'mom_pdf' )


from ROOT import TCanvas
c = TCanvas()
#for (cc,(lab,ind)) in zip(c.pads(1,len(iTag.states().items())),iTag.states().items()) :
for (cc) in c.pads(1,1):
    dataCuts = {} # dict( Cut = '%s == %s' % ( iTag.GetName(), ind ) )
    pdfCuts  = {} # dict( Slice = ( iTag, lab ) )
    obs = [ o for o in observables if hasattr(o,'frame') and o.GetName() != 't' ]
    for (ccc,o) in zip(cc.pads(len(obs)),obs) :
        from P2VV.GeneralUtils import plot
        plot( ccc, o, data, pdf, addPDFs = [ mcpdf ]
                               , dataOpts = dict( MarkerSize = 0.8, MarkerColor = RooFit.kBlack, **dataCuts )
                               , pdfOpts  = dict( LineWidth = 2, **pdfCuts)
                               , addPDFsOpts = [ dict( LineWidth = 2, LineColor = RooFit.kRed, **pdfCuts) ]
                               )
        #mom_pdf.plotOn( f, LineColor = RooFit.kGreen)
