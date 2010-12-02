from ROOT import *
gSystem.Load('libp2vv')

def createMomentPDF(w,name,data,moments) :
    if not moments : return None
    allObs = moments[0].basis().getObservables(data)
    for i in range( data.numEntries() ) :
        allObs.assignValueOnly( data.get(i) )
        for m in moments : m.inc()
    coef = RooArgList()
    fact = RooArgList()
    for m in moments :
        print m.coefficient()
        print m.basis()
        C = 'C_%f' % m.coefficient()
        w.factory( '%s[%f]'%(C,m.coefficient() ) )
        coef.add( w.var( C ) )
        fact.add( m.basis() )
    return RooRealSumPdf(name,name,fact,coef)

#def determineMoments(fname='p2vv_3.root',dataName = 'pdfData', wsname='w') :
fname = 'p2vv_3.root'
dataName = 'pdfData'
wsname = 'w'

f = TFile(fname)
w = f.Get(wsname)
data = w.data(dataName)

ab = abasis(w,'cpsi','ctheta','phi')
moments = []
for i in range(5) :
   for l in range(5) :
      for m in range(-l,l+1) :
           #  Warning: the Y_lm are orthonormal, but the P_i are orthogonal, with dot product 2/(2*i+1)
           moments.append( Moment( ab('mom',i,0,l,m,1.), float(2*i+1)/2 ) )

pdf = createMomentPDF( w, "bkg_angles_pdf", data, moments )
c = TCanvas()
c.Divide(3,1)
for (v,i) in zip( ['cpsi','ctheta','phi'], range(1,100) )  :
    c.cd(i)
    frame = w.var(v).frame()
    data.plotOn(frame)
    pdf.plotOn(frame)
    frame.Draw()

