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
        e =  ( p[0]+p[2] )/2
        #e =  ( p[0]+p[1] )/2
        #e =  ( p[0]+0.5*p[1]+p[2])/3
        # e =  (y[0][0]+y[1][-1+1])/2
        #assert e < 1.0001 and e>= 0
        #return random() < e
        return random() < ( p[0]+0.2*p[1]+0.4*p[2] )/3 * ( y[0][0] + 0.1*y[1][-1+1] + 0.2*y[1][0+1] + 0.3*y[1][1+1] )


from P2VV.RooFitWrappers import *

ws = RooObject( workspace = 'myws' )

from P2VV.Parameterizations.AngularFunctions import JpsiphiHelicityAngles as HelAngles, JpsiphiTransversityAngles as TrAngles
#angles    = HelAngles( cpsi = 'helcthetaK', ctheta = 'helcthetaL', phi = 'helphi' )
angles    = TrAngles( cpsi = 'trcospsi', ctheta = 'trcostheta', phi = 'trphi' )
observables = [ i for i in angles.angles.itervalues() ] 
for i in angles.angles.itervalues() : i.setBins(16)



# build new PDFs with angular coefficients
from P2VV.Parameterizations.AngularPDFs import AngleBasis_AngularPdfTerms
indices  = [ ( PIndex, YIndex0, YIndex1 ) for PIndex in range(3) for YIndex0 in range(3) for YIndex1 in range( -YIndex0, YIndex0 + 1 ) ]
indices = [ (0,0,0, 1), (1,0,0,0.5), (2,1,0,0.1) ]
from math import sqrt
indices = [ ( 0,  0,  0 , 4. ), ( 0,  2,  0 , -sqrt( 16. / 5.)), ( 2,  0,  0 ,   8.) , ( 2, 2,  0 ,  -sqrt( 64. / 5.)) ]
cnvrtInd = lambda i : ('%s'%i).replace('-','m') 
coefPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angles.angles
                                          , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                     , {  'Name'    : 'Cab%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                        , 'Value'   : inds[3]
                                                        , 'MinMax'  : ( -10., +10. )
                                                        , 'Indices' : inds
                                                       }
                                                    ) for inds in indices
                                                  )
                                         )
mcpdf = coefPDFTerms.buildSumPdf('angCoefsPDF')
#mcpdf.Print("T")

eindices = [ (0,0,0,1),(1,0,0,0.5),(2,0,0,1) ]
effPDFTerms = AngleBasis_AngularPdfTerms(  Angles = angles.angles
                                          , **dict( (  'C%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                     , {  'Name'    : 'C_ab%d%d%s' % ( inds[0], inds[1], cnvrtInd(inds[2]) )
                                                        , 'Value'   : inds[3]
                                                        , 'MinMax'  : ( -10., +10. )
                                                        , 'Indices' : inds[0:3]
                                                       }
                                                    ) for inds in eindices
                                                  )
                                         )
epdf = effPDFTerms.buildSumPdf('angeCoefsPDF')
epdf.Print("T")


data = mcpdf.generate(observables,10000)
print 'got dataset with %s entries' % data.numEntries()
allObs = mcpdf.getObservables( data.get() )
allObs.Print()
#print mcpdf.Observables()

#assert False
eps = efficiency( angles.angles['cpsi'],angles.angles['ctheta'],angles.angles['phi'] )
inEffData = RooDataSet( "inEffData","inEffData", allObs )
for event in  data :
    allObs.assignValueOnly( event )
    if eps.accept() : inEffData.add( allObs )
origdata = data
data = inEffData;


#print 'computing efficiency moments'
###################################
from P2VV.GeneralUtils import RealMomentsBuilder
eff = RealMomentsBuilder()
indices  = [ ( i, l, m ) for i in range(4)
                         for l in range(2)
                         for m in range(-l,l+1) ]
eff.appendPYList( angles.angles, indices, PDF = mcpdf, NormSet = allObs)
###################################
#class abasis :
#    def __init__(self,w,*args) :
#       self.w = w
#       (cpsi,cheta,phi) = args if len(args)==3 else args[0]

#       def _f(x) :
#            if type(x) is str : return w[x]
#            if not self.w.function(x.GetName()) : x = w.put(x)
#            return x
#       self.cpsi   = _f(cpsi)
#       self.ctheta = _f(cheta)
#       self.phi    = _f(phi)
#       print 'using %s,%s,%s' % (self.cpsi,self.ctheta,self.phi)

#    def angles(self) : 
#        return RooArgList(self.cpsi,self.ctheta,self.phi)
#    def build(self,label,i,j,k,l,c) :
#        name = "%s_%d_%d_%d_%d" % (label,i,j,k,l)
#        name.replace("-","m")
#        b = self.w.function(name) # workaround a bug in ROOT 5.26 -- if name not present, w.obj(name) will SEGV...
#        from ROOT import RooP2VVAngleBasis
#        if not b : 
#            b = self.w.put( RooP2VVAngleBasis(name,name,self.cpsi,self.ctheta,self.phi,i,j,k,l,c) )
#        return b

#ab = abasis(mcpdf.ws(),angles.angles['cpsi']._var,angles.angles['ctheta']._var,angles.angles['phi']._var)
## if we want to write it as efficiency, i.e. eps_ijk * P_i * Y_jk * PDF then we need the marginal..
## Warning: the Y_lm are orthonormal, but the P_i are orthogonal, but the dot product is (2*i+1)/2
#from ROOT import EffMoment
#moments = [ EffMoment( ab.build("mom",i,0,l,m,1. ),float(2*i+1)/2, mcpdf._var, allObs ) for (i,l,m)  in indices ]

# loop over all data, determine moments
#def computeMoments( data, moments ) :
#    if not moments : return None
#    from ROOT import std
#    vecmom = std.vector('IMoment*')()
#    for m in moments : vecmom.push_back(m)
#    from ROOT import _computeMoments
#    return _computeMoments( data, vecmom )

#computeMoments(data,moments)
#for i in moments : i.Print()

##################################
eff.compute(data)
from math import sqrt,pi
eff.Print( MinSignificance = 0., Names = 'p2vvab.*', Scales = ( 1. / (2*sqrt(pi)), 1. / (2*sqrt(pi)), 1. ) )


pdf = eff * mcpdf
print pdf.GetName(), type(pdf)
pdf.Print("T")
#pdf2.Print("T")

# create generic PDF to describe the angular distribution
#moms = RealMomentsBuilder()
#moms.appendPYList( angles.angles, indices )
#moms.compute(data)
#mom_pdf = moms.createPDF( Name = 'mom_pdf' )


from ROOT import TCanvas
c = TCanvas()
obs = [ o for o in observables if hasattr(o,'frame') and o.GetName() != 't' ]
for (ccc,o) in zip(c.pads(len(obs)),obs) :
    from P2VV.GeneralUtils import plot
    plot( ccc, o, data, pdf, addPDFs = [ mcpdf ]
                           , dataOpts = dict( MarkerSize = 0.8, MarkerColor = RooFit.kBlack )
                           , pdfOpts  = dict( LineWidth = 2 , LineColor = RooFit.kGreen)
                           , addPDFsOpts = [ dict( LineWidth = 2, LineColor = RooFit.kRed ) ]
                           )
    #mom_pdf.plotOn( f, LineColor = RooFit.kGreen)
