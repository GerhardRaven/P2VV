from ROOT import RooArgSet, RooArgList, RooDataSet, RooWorkspace, RooFitResult, RooFit
from ROOT import gStyle,gROOT
gStyle.SetPalette(1)
gROOT.SetStyle("Plain")

# needed to get RooFit.Name, RooFit.Components.... kRed
# how to get just the RooFit namespace ?
#from ROOT import * 

###### decorate TPad with pads...
from ROOT import TPad
def _pads(p,n=None,m=None) :
    if n : 
        if m : p.Divide(n,m)
        else : p.Divide(n)
    i=1
    while p.GetPad(i) :
        yield p.cd(i) 
        i=i+1

TPad.pads = _pads
##########################################

from ROOT import RooPlot
_RooPlot_Draw = RooPlot.Draw
def _RooPlotDraw(self,*args,**kw) :
   pad = kw.pop('pad') if 'pad' in kw else None
   if kw : raise RuntimeError('unknown keyword arguments: %s' % kw )
   if pad : pad.cd()
   _RooPlot_Draw(self,*args)
   if pad : pad.Update()
RooPlot.Draw = _RooPlotDraw

# RooDataSet functions
def _RooDataSetIter(self) :
    for i in range( self.numEntries() ) : yield self.get(i)

from ROOT import RooDataSet
RooDataSet.__iter__ = _RooDataSetIter

# RooAbsCollection/RooArgSet/RooArgList functions
def _RooAbsCollectionIter(self) :
    z = self.createIterator()
    while True :
        a = z.Next()
        if not a : return
        yield a

from ROOT import RooAbsCollection
RooAbsCollection.__iter__ = _RooAbsCollectionIter
RooAbsCollection.__len__  = lambda s   : s.getSize()
RooAbsCollection.__contains__  = lambda s,i : s.contains(i)
RooAbsCollection.__iadd__ = lambda s,x : s if s.add(x)    else s  # else None??
RooAbsCollection.__isub__ = lambda s,x : s if s.remove(x) else s  # else None??
RooAbsCollection.nameList = lambda s : [ j.GetName() for j in s ] 
RooAbsCollection.names    = lambda s : ','.join( s.nameList() )
RooAbsCollection.__eq__   = lambda s,x : s.equals(x)
RooAbsCollection.__ne__   = lambda s,x : not s.equals(x)

def _RooTypedUnary2Binary( t,op ) :
    return lambda x,y : getattr(t,op)(t(x),y)

from ROOT import RooArgSet, RooArgList
for t in [ RooArgSet,RooArgList ] :
    t.__sub__  = _RooTypedUnary2Binary( t, '__isub__' )
    t.__add__  = _RooTypedUnary2Binary( t, '__iadd__' )



# RooWorkspace functions

from ROOT import RooWorkspace
RooWorkspace.__getitem__ = lambda s,i : s.obj(i)
RooWorkspace.__contains__ = lambda s,i : bool( s.obj(i) )
#RooWorkspace.__setitem__ = lambda s,k,v : s.put('%s[%s]'%(k,v))
def _RooWorkspacePut( self ,x ) :
    _import = getattr(RooWorkspace,'import')
    if _import(self,x,RooFit.Silence()) : return None
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


# RooFitResult functions

def _RooFitResultGet(self, parList) :
  # get parameter indices in fit result
  indices = {}
  floatPars = self.floatParsFinal()
  for par in parList :
    index = floatPars.index(par)
    if index >= 0 :
      indices[par] = index
    else :
      print 'P2VV - ERROR: RooFitResult::result(): fit result does not contain parameter', par
      return None

  covMatrix = self.covarianceMatrix()

  values = tuple([floatPars[indices[par]].getVal() for par in parList])
  covariances = tuple([tuple([covMatrix[indices[row]][indices[col]]\
      for col in parList]) for row in parList])

  return (tuple(parList), values, covariances)

RooFitResult.result = _RooFitResultGet





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
def plot( c, obs, data, pdf, components, frameOpts = (), dataOpts = (), pdfOpts = (), logy = False, normalize = True, symmetrize = True, usebar =True ) :
    from ROOT import TLine, TPad, RooFit
    #
    _obs = obs.frame( *frameOpts )  if frameOpts else obs.frame()
    _stash.append(_obs)
    data.plotOn(_obs,RooFit.Name('data'),*dataOpts)
    if components :
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
    #TODO: if normalize : plot rh as a filled histogram with fillcolor blue...
    #      or, maybe, with the 'bar chart' options: 'bar' or 'b'
    for i in dataOpts : 
        if i.opcode() == 'MarkerSize'  : rh.SetMarkerSize(  i.getDouble(0) )
        if i.opcode() == 'MarkerStyle' : rh.SetMarkerStyle( i.getInt(0) )
        if i.opcode() == 'MarkerColor' : rh.SetMarkerColor( i.getInt(0) )
        if i.opcode() == 'Title' : rp.SetTitle( i.getString(0) )
    # rp.addPlotable( rh, 'p' if not usebar else 'b' )
    # zz.plotOn(f,RooFit.DrawOption('B0'), RooFit.DataError( RooAbsData.None ) )
    #rp.SetBarWidth(1.0)
    #rh.SetDrawOption("B HIST")

    rp.addPlotable( rh), 'p'  # , 'B HIST' )
    #rp.setDrawOptions(rh.GetName(),'B')
    if symmetrize :
        m  = max( abs( rh.getYAxisMin() ),abs( rh.getYAxisMax() ) )
        rp.SetMaximum(  m )
        rp.SetMinimum( -m )
    if normalize :
        if rh.getYAxisMin() > -5 : rp.SetMinimum(-5)
        if rh.getYAxisMax() <  5 : rp.SetMaximum(5)
    xa = rp.GetXaxis()
    l = TLine(xa.GetXmin() ,0,xa.GetXmax() ,0)
    from ROOT import kRed
    l.SetLineColor(kRed)
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

def writeCorrMatrixLatex(roofitresult,fname = 'corrtable.tex'):
    parlist = roofitresult.floatParsFinal()
    npar = parlist.getSize()
    corr = []
    for i in range(npar):
        corr.append(roofitresult.correlation(parlist[i]))

    layoutstring = '|c'*(npar+1)
    layoutstring += '|}'

    string = '\\documentclass{article}\n'
    string += '\\begin{document}\n'
    string += '\\begin{table}[h]\n'
    string += '\\begin{center}\n'
    string += '\\begin{tabular}{'
    string += layoutstring+'\n' 
    string += '\\hline\n'
    string += 'Parameter '
    for i in range(npar): string += ' & %s'%(parlist[i].GetName())
    string += ' \\\\\n'
    string += '\\hline\n\\hline\n'


    for j in range(npar):
        string += parlist[j].GetName() 
        for i in range(npar):
            if i>=j:
                if abs(corr[j][i].getVal()) < 0.005:
                    string += ' & -- '
                else:
                    string += ' & ' + str(round(corr[j][i].getVal(),3)) 
            else : 
                string += ' & '
        string += ' \\\\\n'
    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{center}\n'
    string += '\\end{table}\n'
    string += '\\end{document}\n'

    if fname :
        f = open(fname)
        f.write(string)
        f.close

    return {'string':string}


def writeFitParamsLatex(result,name,toys):
    f = open('params_%s.tex'%name,'w')

    string = '\\documentclass{article}\n'
    string += '\\begin{document}\n'

    if toys:
        string += '\\begin{tabular}{|c|c|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result & original value & $\sigma$ from original \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(result.floatParsFinal(),result.floatParsInit()):
            string += '%s & '%i.GetName()
            string += '%s $\pm$ %s & '%(round(i.getVal(),3),round(i.getError(),3))
            string += '%s & '%round(j.getVal(),3)
            string +=  '%s \\\\ \n'%(round((j.getVal()-i.getVal())/i.getError(),3))
    else:
        string += '\\begin{tabular}{|c|c|}\n'
        string += '\\hline\n'
        string += 'parameter & result \\\\ \n'
        string += '\\hline\n'
        string += '\\hline\n'
        for i,j in zip(result.floatParsFinal(),result.floatParsInit()):
            string += '%s & '%i.GetName()
            string += '%s $\pm$ %s \\\\ \n'%(round(i.getVal(),3),round(i.getError(),3))

    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{document}\n'

    f.write(string)
    f.close()
    

def MakeProfile(name,data,pdf,npoints,param1,param1min,param1max,param2,param2min,param2max,NumCPU=8,Extend=True):
    print '**************************************************'
    print 'making profile for %s and %s'%(param1.GetName(),param2.GetName())
    print '**************************************************'

    nll = pdf.createNLL(data,RooFit.NumCPU(NumCPU),RooFit.Extended(Extend))
    profile = nll.createProfile(RooArgSet( param1,param2))
    return profile.createHistogram(name,         param1, RooFit.Binning(npoints,param1_min,param1_max)
                                  , RooFit.YVar( param2, RooFit.Binning(npoints,param2_min,param2_max))
                                  , RooFit.Scaling(False)
                                  )
