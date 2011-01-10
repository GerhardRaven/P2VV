from ROOT import RooArgSet, RooArgList, RooDataSet
from ROOT import gStyle,gROOT
gStyle.SetPalette(1)
gROOT.SetStyle("Plain")

# needed to get RooFit.Name, RooFit.Components....
# how to get just the RooFit namespace ?
from ROOT import * 

def _RooDataSetIter(self) :
    for i in range( self.numEntries() ) :
        yield self.get(i)

RooDataSet.__iter__ = _RooDataSetIter

def _RooArgSetIter(self) :
    z = self.createIterator()
    while True :
        a = z.Next()
        if not a : return
        yield a



RooArgList.__iter__ = _RooArgSetIter
RooArgList.__len__  = lambda s   : s.getSize()
RooArgList.__contains__  = lambda s,i : s.contains(i)
#RooArgSet.__iadd__ = lambda s,y : s.add(y)
#RooArgSet.__radd__ = lambda y,s : s.add(y)
#RooArgSet.__add__  = lambda x,y : RooArgSet(x).__iadd__(y)
RooArgList.nameList = lambda s : [ j.GetName() for j in s ] 
RooArgList.names    = lambda s : ','.join( s.nameList() )

RooArgSet.__iter__  = _RooArgSetIter
RooArgSet.__len__   = lambda s   : s.getSize()
RooArgSet.__contains__  = lambda s,i : s.contains(i)
RooArgSet.nameList  = lambda s : [ j.GetName() for j in s ]
RooArgSet.names     = lambda s : ','.join( s.nameList() )


from ROOT import RooWorkspace
RooWorkspace.__getitem__ = lambda s,i : s.obj(i)
RooWorkspace.__contains__ = lambda s,i : s.obj(i) is not None
#RooWorkspace.__setitem__ = lambda s,k,v : s.put('%s[%s]'%(k,v))
def _RooWorkspacePut( self ,x ) :
    _import = getattr(RooWorkspace,'import')
    if _import(self,x) : return None
    return self[x.GetName()]
RooWorkspace.put = _RooWorkspacePut

def setConstant(ws, pattern, constant = True, value = None):
    import re
    nrexp = re.compile(pattern)
    rc = 0
    for arg in ws.allVars() :
        if not nrexp.match(arg.GetName()) : continue
        arg.setConstant( constant )
        if constant and value :
            if value < arg.getMin() : arg.setMin(value) 
            if value > arg.getMax() : arg.setMax(value) 
            arg.setVal(value) 
        rc += 1
    return rc


RooWorkspace.setConstant = setConstant

# plot -- example usage:
# _c1 = plot( c.cd(1),mpsi,data,pdf
#           , { 'psi'    : ( RooFit.LineColor(RooFit.kGreen),RooFit.LineStyle(kDashed) ) 
#             , 'nonpsi' : ( RooFit.LineColor(RooFit.kBlue),RooFit.LineStyle(kDashed) )
#             }
#           , frameOpts = ( RooFit.Bins(30), )
#           , dataOpts = ( RooFit.MarkerSize(0.4), RooFit.XErrorSize(0) )
#           , pdfOpts = ( RooFit.LineWidth(2), ) 
#           )
global  _stash
_stash = [] #keep the relevant objects alive by keeping a reference to them
def plot( c, obs, data, pdf, components, frameOpts = None, dataOpts = None, pdfOpts = None, logy = False, normalize = True, symmetrize = True ) :
    from ROOT import TLine, TPad
    #
    _obs = obs.frame( *frameOpts )  if frameOpts else obs.frame()
    _stash.append(_obs)
    data.plotOn(_obs,RooFit.Name('data'),*dataOpts)
    for comp,opt in components.iteritems() :
        z = opt + pdfOpts
        pdf.plotOn(_obs,RooFit.Components(comp), *z )
    pdf.plotOn(_obs,RooFit.Name('pdf'),*pdfOpts)
    _obs.drawAfter('pdf','data')
    #TODO: add chisq/nbins
    #chisq = _obs.chiSquare('pdf','data')
    #nbins = _obs.GetNbinsX()
    rh = _obs.residHist('data','pdf',normalize)

    _stash.append(rh)
    xa = _obs.GetXaxis()
    rh.GetXaxis().SetLimits(xa.GetXmin(),xa.GetXmax())
    rp = _obs.emptyClone( _obs.GetName() + '_resid' )
    _stash.append(rp)
    if logy : _obs.SetMinimum(0.1)
    for i in dataOpts : 
        if i.opcode() == 'MarkerSize'  : rh.SetMarkerSize(  i.getDouble(0) )
        if i.opcode() == 'MarkerStyle' : rh.SetMarkerStyle( i.getInt(0) )
        if i.opcode() == 'MarkerColor' : rh.SetMarkerColor( i.getInt(0) )
        if i.opcode() == 'Title' : rp.SetTitle( i.getString(0) )
    rp.addPlotable( rh, 'p' )
    if symmetrize :
        m  = max( abs( rh.getYAxisMin() ),abs( rh.getYAxisMax() ) )
        rp.SetMaximum(  m )
        rp.SetMinimum( -m )
    xa = rp.GetXaxis()
    l = TLine(xa.GetXmin() ,0,xa.GetXmax() ,0)
    rp.addObject(l)
    #TODO: improve (remove?) axis labels from rp, move up against the initial plot

    # only now start drawing...
    c.cd()
    hname = obs.GetName() + '_plot1'
    h = TPad(hname,hname,0,0.2,1,1)
    _stash.append(h)
    if logy: h.SetLogy(1)
    h.SetNumber(1)
    h.Draw()
    c.cd(1)
    _obs.Draw()

    c.cd()
    rname = obs.GetName() + '_resid1'
    r = TPad(rname,rname,0,0,1,0.2)
    _stash.append(r)
    r.SetNumber(2)
    r.Draw()
    c.cd(2)
    rp.Draw()

    c.Update()
    return c
