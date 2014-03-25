
###########################################################################################################################################
## Utilities.Plotting: P2VV utilities for making plots                                                                                   ##
##                                                                                                                                       ##
## authors:                                                                                                                              ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl                                                                       ##
##   VS,  Vasilis Syropoulos, Nikhef,      v.syropoulos@nikhef.nl                                                                        ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl                                                                       ##
##   RA,  Roel Aaij,          Nikhef,      r.aaij@nikhef.nl                                                                              ##
##                                                                                                                                       ##
###########################################################################################################################################

# plot stash: keep the relevant objects alive by keeping a reference to them
global _P2VVPlotStash
_P2VVPlotStash = []


# plotting function
def plot(  canv, obs, data = None, pdf = None, addPDFs = [ ], components = None, xTitle = '', yTitle = '', xTitleOffset = None
           , yTitleOffset = None, yScale = ( None, None ), yScaleRel = ( None, None ), frameOpts = { }, dataOpts = { }, pdfOpts = { }
           , addPDFsOpts = [ { } ], plotResidHist = False, logy = False, logx = False, normalize = True, symmetrize = True, usebar = True
        ) :
    """makes a P2VV plot

    example usage:

    canvas = plot( canvas.cd(1), observable, data, pdf
                  , {  'psi'    : { 'LineColor' : RooFit.kGreen, 'LineStyle' : RooFit.kDashed }
                     , 'nonpsi' : { 'LineColor' : RooFit.kBlue,  'LineStyle' : RooFit.kDashed }
                    }
                  , xTitle = 'M (MeV/c)'
                  , frameOpts = { 'Title'      : 'B mass', 'Bins        : 30 }
                  , dataOpts  = { 'MarkerSize' : 0.4,      'XErrorSize' : 0  }
                  , pdfOpts   = { 'LineWidth'  : 2                           }
                 )
    """
    from ROOT import TLine, TPad

    # create frame for observable
    obsFrame = obs.frame(**frameOpts)  if frameOpts else obs.frame()
    xAxis = obsFrame.GetXaxis()
    yAxis = obsFrame.GetYaxis()
    _P2VVPlotStash.append(obsFrame)

    frames = []
    frames.append(obsFrame)

    # plot data
    if data :
        rooPlot = data.plotOn( obsFrame, Name = 'data', **dataOpts )
        # Set negative bins to 0 if logy is requested
        if logy:
            minimum = 0.
            hist = rooPlot.getHist()
            from ROOT import Double
            x = Double(0.)
            y = Double(0.)
            for i in range(hist.GetN()):
                r = hist.GetPoint(i, x, y)
                if y < minimum:
                    minimum = y
            #hist.SetMinimum(minimum + 0.1)
            obsFrame.SetMinimum( max( ( minimum, 0.1 ) ) )

    # plot PDF
    if pdf :
        # define function that parces the 'Slice(s)' argument and plots the pdf
        def plotPDFWithSlices( pdf, frame, name, **pdfOpts ) :
            if 'Slice' in pdfOpts or 'Slices' in pdfOpts :
                # get 'Slice(s)' argument from plot options
                origSlices = pdfOpts.pop( 'Slices', [ ] )
                if 'Slice' in pdfOpts : origSlices += [ pdfOpts.pop('Slice') ]

                # parse 'Slices' argument
                slicesList = [ [ ] ]
                for slice in origSlices :
                    tempList = [ ]
                    for slices in slicesList : tempList += [ slices + [( slice[0], catType.strip() )] for catType in slice[1].split(',') ]
                    slicesList = tempList

                for num, slices in enumerate(slicesList) :
                    # plot pdf for all slices
                    if num == 0 and len(slicesList) == 1 :
                        opts = dict( Name = name, Slices = slices, **pdfOpts )
                    elif num == 0 :
                        opts = dict( Name = name + '0', Invisible = None, Slices = slices, **pdfOpts )
                    elif num == len(slicesList) - 1 :
                        opts = dict( Name = name, AddTo = (name + '%d' % (num - 1), 1., 1.), Slices = slices, **pdfOpts )
                    else :
                        opts = dict(  Name = name + '%d' % num, AddTo = (name + '%d' % (num - 1), 1., 1.), Invisible = None
                                    , Slices = slices, **pdfOpts
                                   )

                    pdf.plotOn( obsFrame, **opts )

            else :
                pdf.plotOn( obsFrame, Name = name, **pdfOpts )

        if components :
            # plot separate components of the pdf
            for num, comp in enumerate( components.keys() ) :
                drawOpts = components[comp].copy()
                for opt, optVal in pdfOpts.iteritems() :
                    if opt not in drawOpts : drawOpts[opt] = optVal
                plotPDFWithSlices( pdf, obsFrame, 'comp%d' % num, Components = comp, **drawOpts )

        # plot total pdf
        drawOpts = pdfOpts.copy()
        plotPDFWithSlices( pdf, obsFrame, 'pdf', **drawOpts )

        # draw data after drawing the PDF
        if data and 'Asymmetry' not in pdfOpts : obsFrame.drawAfter( 'data', 'pdf' )

    # plot additional PDFs
    if addPDFs :
        for num, addPDF in enumerate(addPDFs) :
            addPDF.plotOn( obsFrame, Name = 'addPDF' + str(num), **(addPDFsOpts[num]) )
            if data and 'Asymmetry' not in addPDFsOpts[num] : obsFrame.drawAfter( 'data', 'addPDF' + str(num) )

    #TODO: add chisq/nbins
    #chisq = obsFrame.chiSquare( 'pdf', 'data' )
    #nbins = obsFrame.GetNbinsX()

    # set y scale
    if yScale[0]    != None : obsFrame.SetMinimum(yScale[0])
    if yScale[1]    != None : obsFrame.SetMaximum(yScale[1])
    if yScaleRel[0] != None : obsFrame.SetMinimum( yScaleRel[0] * obsFrame.GetMinimum() )
    if yScaleRel[1] != None : obsFrame.SetMaximum( yScaleRel[1] * obsFrame.GetMaximum() )
    if logy and obsFrame.GetMinimum() <= 0 : obsFrame.SetMinimum(0.1)

    # set axis titles
    if xTitle : xAxis.SetTitle(xTitle)
    if yTitle : yAxis.SetTitle(yTitle)

    # set axis title offsets
    if yTitleOffset: yAxis.SetTitleOffset(yTitleOffset)
    if xTitleOffset: xAxis.SetTitleOffset(xTitleOffset)

    # get residuals histogram
    if plotResidHist and data and pdf :
        residHist = obsFrame.residHist( 'data', 'pdf', normalize )
        residHist.GetXaxis().SetLimits( xAxis.GetXmin(), xAxis.GetXmax() )

        _P2VVPlotStash.append(residHist)

        xAxis.SetLabelOffset(0.1)
        #yAxis.SetTitleSize(0.10)
        #yAxis.SetLabelSize(0.08)
        yAxis.SetTitleOffset( 0.7 * yAxis.GetTitleOffset() )

        # create residuals frame
        residFrame = obsFrame.emptyClone( obsFrame.GetName() + '_resid' )
        #if 'Title' in frameOpts: residFrame.SetTitle(frameOpts['Title'])
        #residFrame.SetTitle('')
        # set axis titles
        xAxis = residFrame.GetXaxis()
        xAxis.SetLabelSize(0.15)
        xAxis.SetTitleSize(0.15)
        xAxis.SetLabelOffset(0.02)
        xAxis.SetTitleOffset(1.1)
        yAxis = residFrame.GetYaxis()
        yAxis.SetTitle('')
        yAxis.SetLabelSize(0.13)
        yAxis.SetLabelOffset(0.01)
        _P2VVPlotStash.append(residFrame)
        frames.append(residFrame)

        if xTitle:
            xAxis.SetTitle(obsFrame.GetXaxis().GetTitle())

        # set minimum for observable's frame if there is a log scale for y
        #if logy : obsFrame.SetMinimum(0.1)

        # set residual plot options
        #TODO: if normalize : plot residHist as a filled histogram with fillcolor blue...
        #      or, maybe, with the 'bar chart' options: 'bar' or 'b'
        if dataOpts :
            fun = {  'MarkerSize'  : lambda x : residHist.SetMarkerSize(x)
                   , 'MarkerStyle' : lambda x : residHist.SetMarkerStyle(x)
                   , 'MarkerColor' : lambda x : residHist.SetMarkerColor(x)
                   , 'LineWidth'   : lambda x : residHist.SetLineWidth(x)
                   , 'Title'       : lambda x : residFrame.SetTitle(x)
                  }
            for k, v in dataOpts.iteritems() :
                if k in fun : fun[k](v)

        # residFrame.addPlotable( residHist, 'p' if not usebar else 'b' )
        # zz.plotOn(f,RooFit.DrawOption('B0'), RooFit.DataError( RooAbsData.None ) )
        ## Draw residuals as a histogram, set the bar width a bit wider than one so they look nice.
        ## Now to fix the root bug to make it work properly with a variable binning.
        residFrame.SetBarWidth(1.02)
        residFrame.addPlotable( residHist, 'BX' if not type(plotResidHist) == str else plotResidHist )

        if symmetrize :
            # symmetrize y-axis residuals histogram
            maxY = max( abs(residHist.getYAxisMin()), abs(residHist.getYAxisMax()) )
            residFrame.SetMaximum(maxY)
            residFrame.SetMinimum(-maxY)

        if normalize :
            if residHist.getYAxisMin() > -5.5 : residFrame.SetMinimum(-5.5)
            if residHist.getYAxisMax() < +5.5 : residFrame.SetMaximum(+5.5)

        # add a line at y=0
        zeroLine = TLine( xAxis.GetXmin(), 0, xAxis.GetXmax(), 0 )
        from ROOT import kRed
        zeroLine.SetLineColor(kRed)
        residFrame.addObject(zeroLine)
        #TODO: improve (remove?) axis labels from residFrame, move up against the initial plot

        # draw observable frame
        canv.cd()
        obsName = obs.GetName() + '_plot1'
        obsPad = TPad( obsName, obsName, 0, 0.32, 1, 1 )
        _P2VVPlotStash.append(obsPad)
        if logy: obsPad.SetLogy(1)
        if logx: obsPad.SetLogx(1)
        obsPad.SetNumber(1)
        #obsPad.SetLeftMargin(0.12)
        obsPad.SetTopMargin(0.04)
        obsPad.SetBottomMargin(0.04)
        obsPad.Draw()
        canv.cd(1)
        if 'Title' in frameOpts and not frameOpts['Title'] : obsFrame.SetTitle('')
        obsFrame.Draw()

        # draw residuals frame
        canv.cd()
        residName = obs.GetName() + '_resid1'
        residPad = TPad( residName, residName, 0, 0, 1, 0.32 )
        if logx: residPad.SetLogx(1)
        _P2VVPlotStash.append(residPad)
        residPad.SetNumber(2)
        #residPad.SetLeftMargin(0.12)
        residPad.SetTopMargin(0.)
        residPad.SetBottomMargin(0.4)
        residPad.Draw()
        canv.cd(2)
        if 'Title' in frameOpts and not frameOpts['Title']:
            residFrame.SetTitle("")
        residFrame.Draw()

    else :
        # draw observable frame
        canv.cd()
        if logy: canv.SetLogy(1)
        if logx: canv.SetLogx(1)
        title = frameOpts.get("Title", "")
        if title:
            obsFrame.SetTitle(title)
        obsFrame.Draw()

    canv.Update()
    return frames


# function for plotting a PDF in KK mass bins
def getCondObsPlotsInKKbins(pdf, data, canv):
    from ROOT import TCanvas, gPad, TH1F
    from P2VV.RooFitWrappers import Category

    condObs = pdf.ConditionalObservables()
    nPads = len(condObs)
    nBins =  pdf.indexCat().numTypes()

    #canv = TCanvas('CondObsCanvInKKbins')
    if (nPads & 1 == 0): canv.Divide(nPads/2, 2)
    else: canv.Divide((nPads+1)/2, 2)
    pad = 1
    #Dictionary with all the histogrms, the form is:  { observable, {'bin_i',hist}  }
    Histos = dict( (obs.GetName(), dict( ('bin{0}'.format(bin), TH1F() ) for bin in xrange(nBins)  ) )  for obs in condObs )

    t = data.buildTree()
    entries = t.GetEntries()

    for obs in condObs:
        canv.cd(pad)
        obsName = obs.GetName()
        for KKbin in xrange(nBins):
            binName = 'bin{0}'.format(KKbin)
            sImpose = 'same' if not KKbin == 0 else ''
            if isinstance(obs, Category): #Is CondObs discreate?  (It makes difference in the way Categories are listed in the TTree)
                if (obsName=='iTagOS' or 'iTagSS'):Histos[obsName][binName].SetBins(3,-1,2)
                else: Histos[obsName][binName].SetBins(2,0,2)
                Histos[obsName][binName].SetDefaultSumw2(True)
                Histos[obsName][binName].SetName(obs.GetName() + 'bin{0}'.format(KKbin))
                Histos[obsName][binName].SetLineColor(KKbin + 1)
                Histos[obsName][binName].SetLineWidth(1)
                Histos[obsName][binName].SetAxisRange(0, 350,'Y')
                Histos[obsName][binName].SetXTitle(obsName)

                for event in xrange(entries):
                    t.GetEntry(event)
                    whichKKbin = t.__getattr__('KKMassCat_idx')
                    if (whichKKbin == KKbin ):
                        Histos[obsName][binName].Fill( t.__getattr__(obsName+'_idx'), t.__getattr__('weightVar') )
            else:
                if (obsName=='sigmat'): Histos[obs.GetName()][binName].SetBins(15,0,0.12)
                else:                   Histos[obsName][binName].SetBins(15,0.05,0.55)
                Histos[obsName][binName].SetDefaultSumw2(True)
                Histos[obsName][binName].SetName(obs.GetName() + 'bin{0}'.format(KKbin))
                Histos[obsName][binName].SetLineColor(KKbin + 1)
                Histos[obsName][binName].SetLineWidth(1)
                Histos[obsName][binName].SetAxisRange(.01, 1e3,'Y')
                Histos[obsName][binName].SetXTitle(obsName)

                for event in xrange(entries):
                    t.GetEntry(event)
                    whichKKbin = t.__getattr__('KKMassCat_idx')
                    if (whichKKbin == KKbin ):
                        Histos[obsName][binName].Fill( t.__getattr__(obsName), t.__getattr__('weightVar') )


            int1 = Histos[obsName][binName].GetSumOfWeights()
            if     KKbin==0 : int0 = Histos[obsName]['bin0'].GetSumOfWeights()
            if not KKbin==0 : Histos[obsName][binName].Scale(int0 / int1)

            Histos[obsName][binName].Draw(sImpose)

        if not isinstance(obs, Category): gPad.SetLogy()
        pad += 1
        #assert(False)
    #return  Histos
    if obsName=='iTagSS': assert(False)

    return canv

def compareDataSets( canv, obs, data={}, dataOpts={}, frameOpts={}, **kwargs):  
    """ Superimpose RooDataSets
    example usage:  compareDistributions(     data = dict( dataSet1   = mcBefore
                                                           dataSet2   = mcAfter                                  ),
                                          dataOpts = dict( dataSet1   = dict( MarkerColor = 2, MarkerSize = 1 ),
                                                           dataSet2   = dict( MarkerColor = 4, MarkerSize = 1 )  )
                                         frameOpts = dict( Bins = 100, Range=( lo, high ) )             
                                                           """
    logy          = kwargs.pop( 'logy',           False ) 
    RangeY        = kwargs.pop( 'RangeY',            [] )
    titleY        = kwargs.pop( 'titleY',            '' )
    statOn        = kwargs.pop( 'statOn',          False)
    legendEntries = kwargs.pop( 'includeInLegend',   [] )
    Xtitle        = kwargs.pop( 'xTitle',           ''  )
    save          = kwargs.pop('save',             False)

    # get dataset scales (scale down the largest samples)
    entries = [ data[k].sumEntries() for k in data.keys() ]
    for k in data.keys(): dataOpts[k]['Rescale'] = min(entries) / data[k].sumEntries()
        
    # create frame
    obsFrame = obs.frame(**frameOpts)

    # plot data on frame
    for key in data.keys():                                                                                                                  
        if key not in dataOpts.keys():                                                                                                           
            assert  False, 'P2VV - EROR: compareDataSets: data and dataOpts dictionaries must have the same keys.'                              

    if statOn: # make legend with mean and rms                                                                                                                     
        from ROOT import TPaveText                                                                                                               
        legend = TPaveText(.455, .818, .944, .951, 'NDC' )                                                                                       
        legend.SetFillColor(0)                                                                                                                   
        legend.SetShadowColor(0)           
                 
    YaxisMaxima, YaxisMinima = [], []                                                                                                                     
    for d in data.keys():                                                                                                                        
        frame = data[d].plotOn( obsFrame, **dataOpts[d] )                                                                                        
        if statOn and d in legendEntries:
            mean = data[d].meanVar(obs).getVal()                                                                                                 
            rms  = data[d].rmsVar(obs).getVal()                                                                                                  
            legend.AddText( '#color[%s]{#bullet} mean=%.3e, rms=%.3e'%(dataOpts[d]['MarkerColor'],mean,rms) )                                    
        YaxisMaxima += [frame.GetMaximum()]                                                                                                      
        YaxisMinima += [frame.GetMinimum()]                                                                                                      
        del frame

    # title and axis ranges
    if RangeY: obsFrame.SetAxisRange( RangeY[0], RangeY[1], 'Y' )
    else:
        scales = []
        for k in data.keys(): scales += [dataOpts[k]['Rescale']]
        obsFrame.SetAxisRange( min(YaxisMinima),  1.15 * min(scales) * max(YaxisMaxima) , 'Y' )
    obsFrame.SetYTitle(titleY)
    obsFrame.SetTitle('')
    xtitle = Xtitle if Xtitle else obsFrame.GetXaxis().GetTitle().replace('(M','[M').replace('c^2)','c^{2}]') 
    obsFrame.SetXTitle(xtitle)
    obsFrame.SetName( obs.GetName() )

    # set logy
    if logy:
        if obsFrame.GetMinimum() <= 0:
            obsFrame.SetMinimum(.1)
        canv.SetLogy(1)
    
    # draw / save
    canv.cd()
    obsFrame.Draw()
    if statOn: # save intemediate canvas.
        legend.Draw()
        from ROOT import TCanvas
        c = TCanvas( obs.GetName() +'_'+ canv.GetName(), obs.GetName() +'_'+ canv.GetName())
        obsFrame.Draw()
        legend.Draw()
        c.Print(obs.GetName() +'_'+ canv.GetName() + '.pdf')
        if save:
            names = ''
            for n in legendEntries: names += n + '_' 
            saveInRootFile(c, obsFrame.getPlotVar().GetName() +'_'+ names)
        return obsFrame, legend
    return obsFrame, None

def makeAssymetryPlot( canv, frame, refHist, numOfFrames, yRange=[], save=False ):
    # get list of RooHist objects from frame
    from ROOT import RooHist
    HistList = []
    for idx in range( int(frame.numItems()) ): 
        if type(frame.getObject(idx))==RooHist: HistList.append( frame.getObject(idx) )
  
    # grab reference histogram
    for h in HistList: 
        if refHist in h.GetName(): hRef = h
    if not hRef: assert False, 'P2VV - ERROR: makeAssymetryPlot(): Cannot find reference histogram.'
    HistList.remove(hRef)
    
    # get bin content and errors of reference histogram
    ylistRef = hRef.GetY()
    yErrHRef = hRef.GetEYhigh()
    yErrLRef = hRef.GetEYlow()

    # functions for error caclulation 
    from array import array
    from math import sqrt
    _l2a      = lambda List: array('d', List) # python list to array
    _assym    = lambda x0,x1: (x0-x1) / (x0+x1) # assymetry value
    _assymErr = lambda x0,e0,x1,e1: 2*sqrt((x1*e1)**2 + (x0*e0)**2) / (x1+x0)**2 # assymentry error
      
    # make assymetry plots
    from ROOT import TGraphAsymmErrors
    assymPlotsList = []
    for hist in HistList: # loop over histograms present in frame 
        ylist = hist.GetY()      # array of y-axis values
        yErrH = hist.GetEYhigh() # array of y-axis high errors 
        yErrL = hist.GetEYlow()  # array of y-axis low errors
        assymYlist, assymYerrH, assymYerrL = [],[],[] # prepare relevant assymetry lists 
        
        # loop over histogram and calculate per bin assymetry
        for point in xrange(frame.GetNbinsX()):
            if ylistRef[point]==0 and ylist[point]==0 :
                assymYlist.append(.0)
                assymYerrH.append(.0)
                assymYerrL.append(.0)
            else:
                assymYlist.append( _assym(ylistRef[point], ylist[point]) )
                assymYerrH.append( _assymErr(ylistRef[point], yErrHRef[point], ylist[point], yErrH[point]) )
                assymYerrL.append( _assymErr(ylistRef[point], yErrLRef[point], ylist[point], yErrL[point]) )

        graph = TGraphAsymmErrors( frame.GetNbinsX(),   # bins in x-axis
                                   hRef.GetX(),         # array of x-axis values
                                   _l2a(assymYlist),    # array of y-axis values
                                   _l2a( 100*[.0] ),    # array of x-axis high errors 
                                   _l2a( 100*[.0] ),    # array of x-axis low errors 
                                   _l2a( assymYerrH ),  # array of y-axis high errors 
                                   _l2a( assymYerrL )   # array of y-axis low errors 
                                   )
        
        if yRange: 
            graph.SetMinimum(yRange[0])
            graph.SetMaximum(yRange[1])
        graph.SetName(hist.GetName() + frame.getPlotVar().GetName())
        graph.SetTitle(frame.getPlotVar().GetName())
        graph.GetXaxis().SetTitle( frame.GetXaxis().GetTitle() )
        graph.GetYaxis().SetTitle('Normalised Assymetry')
        graph.SetMarkerColor(hist.GetMarkerColor())
        graph.SetMarkerStyle(20)
        markerSize = 1.5 if save else .5
        graph.SetMarkerSize(markerSize)
        assymPlotsList.append(graph)
        _P2VVPlotStash.append(graph)

    canv.cd()
    onePlot = assymPlotsList.pop()
    onePlot.Draw('APZ')

    if save:
        from ROOT import TCanvas
        cw = TCanvas(canv.GetName() + '.pdf')    
        cw.cd()
        onePlot.Draw('APZ')
        histNames = ''
        for h in HistList: histNames += h.GetName().replace('h_','') 
        saveInRootFile(cw,'assym_' + frame.getPlotVar().GetName() +'_'+ histNames)

    for g in assymPlotsList: 
        canv.cd()
        g.Draw('ZPsame')
        if save: 
            cw.cd()
            g.Draw('ZPsame')
            cw.Print('assym_' + frame.getPlotVar().GetName() +'_'+ histNames + '.pdf')
    
    del assymPlotsList
    

# class for plotting a PDF in KK mass bins
class CPcomponentsPlotingToolkit():
    def __init__(self, pdf, data):
        #Initializer builds the CP component pdfs
        #Create objects
        self._data = data
        self._tPdf = pdf
        self._CpCompPdfs = dict(total = self._tPdf)
        self._comps = ['even','odd','swave']
        self._condObservables = self._tPdf.ConditionalObservables()
        self._observables = self._tPdf.Observables() - self._condObservables
        self._CPnormFracs = {}
        self._lineColors = dict(total=1,even=4,odd=4,swave=2)
        self._lineStyles = dict(total=1,even=9,odd=3,swave=5)
        self._lineWidth  = 2

        # Check if KK mass binning feature is active
        from P2VV.RooFitWrappers import SimultaneousPdf
        self._flagKKbin = isinstance( self._tPdf, SimultaneousPdf )
        if (self._flagKKbin):
            self._nKKbins = self._tPdf.indexCat().numTypes()
            self._binNames = [('bin{0}'.format(bin))for bin in xrange(self._nKKbins)]
            self._pdfsSuperDict = {}
            self._sliceNormFracs = {}

        #Start CP components projection of pdf
        # Create pdf paramters dictionary and a set for the original paramters.
        parsDict    = dict((p.GetName(),p) for p in self._tPdf.Parameters()  )
        originalSet = set([ self._tPdf.ws().function("AparMag2"), parsDict['AperpMag2'], parsDict['A0Mag2']])
        if self._flagKKbin: # In the case of KK binning there are several f_S fractions
            f_sSet = set([ parsDict['f_S_bin{0}'.format(b)] for b in xrange(self._nKKbins) ])
        else: f_sSet = set([parsDict['f_S']])
        originalSet.update(f_sSet)

        #Construct CP pdf components
        from P2VV.RooFitWrappers import ConstVar, Customizer
        for Comp in self._comps:
            #Dictionary with the values  of the parameters to be replaced
            replacementDict = {}
            if (Comp == "even"):
                replacementDict.update(dict(A0Mag2    = parsDict['A0Mag2'].getVal(),
                                            AperpMag2 =            0               ,
                                            AparMag2  = 1-parsDict['A0Mag2'].getVal()-parsDict['AperpMag2'].getVal()
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), 0 ) for f_S_i in f_sSet)  )

            if (Comp == "odd"):
                replacementDict.update(dict(A0Mag2    =                  0             ,
                                            AperpMag2 = parsDict['AperpMag2'].getVal() ,
                                            AparMag2  =                  0
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), 0 ) for f_S_i in f_sSet)  )

            if (Comp == "swave"):
                replacementDict.update(dict(A0Mag2    = 0 ,
                                            AperpMag2 = 0 ,
                                            AparMag2  = 0
                                            ))
                replacementDict.update(dict( (f_S_i.GetName(), f_S_i.getVal() ) for f_S_i in f_sSet)  )

            replacementSet = set([ConstVar(Name=k + Comp, Value=v) for k,v in replacementDict.iteritems()])
            CPcompPDF = Customizer( Pdf = self._tPdf, OriginalArgs = originalSet, SubstituteArgs = replacementSet
                                   , ReplaceByName = True, ArgumentSuffix = Comp )

            if (Comp == "even") : self._CpCompPdfs.update(dict(even  = CPcompPDF ))
            if (Comp == "odd")  : self._CpCompPdfs.update(dict(odd   = CPcompPDF ))
            if (Comp == "swave"): self._CpCompPdfs.update(dict(swave = CPcompPDF ))

        if self._flagKKbin:
            #Split the Cp Components further into individual KK mass copmponents
            tPdfs = dict( (bin, self._CpCompPdfs['total'].getPdf(bin)) for bin in self._binNames)
            ePdfs = dict( (bin, self._CpCompPdfs['even' ].getPdf(bin)) for bin in self._binNames)
            oPdfs = dict( (bin, self._CpCompPdfs['odd'  ].getPdf(bin)) for bin in self._binNames)
            sPdfs = dict( (bin, self._CpCompPdfs['swave'].getPdf(bin)) for bin in self._binNames)
        self._pdfsSuperDict.update(dict( (bin, dict(total = tPdfs[bin],
                                                    even  = ePdfs[bin],
                                                    odd   = oPdfs[bin],
                                                    swave = sPdfs[bin] ) ) for bin in self._binNames ))
        #Sneek these lines that set a unit for the phi, it helps RooFit
        #  when creating the RooPlot. for some reason the unit of phi was not set up
        try: list(self._observables)[[a.GetName() for a in list(self._observables)].index('helphi')].setUnit('rad')
        except ValueError: print 'helphi observable not in angles lsit, Failed to set, rad, as a unit. '

    #End of the initialazation

    #Internal methods
    def calculateNormFracs(self, pdfDict):
        from ROOT import RooArgSet
        obs = RooArgSet(o._target_() for o in self._observables )
        totInt = pdfDict['total'].getNorm(obs)
        fEven  = pdfDict['even' ].getNorm(obs) / totInt
        fOdd   = pdfDict['odd'  ].getNorm(obs) / totInt
        fSwave = pdfDict['swave'].getNorm(obs) / totInt
        return dict(even=fEven, odd=fOdd, swave=fSwave)

    def calculateCPnormFracs(self):
        print 'P2VV - INFO: Calculating relative normaliation fractions of the CP components.'
        if not self._flagKKbin:
             print 'P2VV - INFO: Finished calculating relative normaliation fractions of the CP components.'
             return calculateNormFracs(self._CpCompPdfs)
        if self._flagKKbin:
            self._CPnormFracs = dict( (bin ,self.calculateNormFracs(dict(
                                 total = self._CpCompPdfs['total'].getPdf(bin),
                                 even  = self._CpCompPdfs['even'].getPdf(bin),
                                 odd   = self._CpCompPdfs['odd'].getPdf(bin),
                                 swave = self._CpCompPdfs['swave'].getPdf(bin)
                                 )))for bin in self._binNames )
            print 'P2VV - INFO: Finished calculating relative normaliation fractions of the CP components.'
            return self._CPnormFracs

    def calculateKKslicesNormFracs(self):
        table = self._data.table(self._tPdf.indexCat())
        total = float(self._data.sumEntries()) if   self._data.isWeighted() \
                                               else float(self_.data.numEntries())
        self._sliceNormFracs =  dict( (bin,table.get(bin)/total)for bin in self._binNames )
        return self._sliceNormFracs

    def getProJWdata(self,bin,Bins):
        #Helping internal function to aviodavoid dublicating code,
        #  Usefull in the case where you make 6x4 observable plots
        if bin:
            projData = self.binDataSet(Bins)
            from ROOT import RooRealVar
            projVars = []
            for pV in projData.get():
                if isinstance(pV,RooRealVar): projVars.append(pV)
        else :
            projVars = list(self._condObservables)
            if self._flagKKbin: projVars.append(self._tPdf.indexCat())
            projData = self._data.reduce(ArgSet=projVars)

        return dict(data=projData, vars=projVars)

    def binDataSet(self, nBins):
        if self._flagKKbin: projVars = list(self._condObservables) + [self._tPdf.indexCat()]
        else              : projVars = list(self._condObservables)

        from P2VV.RooFitWrappers import Category
        from ROOT import RooArgSet, RooDataHist
        binnedVarsList = []
        #Bin only the continous observables
        for pV in list(self._condObservables):
            if    isinstance(pV,Category):pass
            else: binnedVarsList.append(pV)
        for pV in binnedVarsList: pV.setBins(nBins)

        binnedVars =  RooArgSet(self._tPdf.indexCat(), *binnedVarsList)
        return RooDataHist('RDH', 'RDH', binnedVars, self._data.reduce(RooArgSet(*projVars)))

    #Interface
    def getCPcompPdf(self):        return self._CpCompPdfs
    def getNumKKbins(self):        return self._nKKbins
    def getCPcompPdfKKbins(self):  return self._pdfsSuperDict
    def getKKbinNames(self):       return self._binNames
    def getCpCompNames(self):      return self._comps
    def getCPnormFracs(self):
        if not self._CPnormFracs: self.calculateCPnormFracs()
        return self._CPnormFracs

    def getKKslicesNormFracs(self):
        if not self._sliceNormFracs: self.calculateKKslicesNormFracs()
        return self._sliceNormFracs

    def getPdfOpts(self, BinData=True,bins=20):
        if   BinData: projDataSet=self.binDataSet(bins)
        else:
            projVars = list(self._condObservables) + [self._tPdf.indexCat()]
            projDataSet = self._data.reduce(ArgSet=projVars)
        return dict( LineWidth = self._lineWidth
                   , LineColor = self._lineColors['total']
                   , ProjWData = (projDataSet, False)
                     )

    def getAddPdfs(self):
        return [self._CpCompPdfs['even' ].getPdf(b)for b in self._binNames] +\
               [self._CpCompPdfs['odd'  ].getPdf(b)for b in self._binNames] +\
               [self._CpCompPdfs['swave'].getPdf(b)for b in self._binNames]

    def getAddPdfsOpts(self, BinData=True,bins=20):
        if not self._CPnormFracs:    self.calculateCPnormFracs()
        if not self._sliceNormFracs: self.calculateKKslicesNormFracs()
        if BinData:
            data     = self.binDataSet(bins)
            projVars = self.getProJWdata(BinData,bins)['vars']
        else:
            data     = self._data
            projVars = list(self._condObservables)
        opts = []
        for comp in self._comps:
            for bin in self._binNames:
                binInd = self._binNames.index(bin)
                addPdfOpt_i = dict( ProjWData     = (data.reduce('KKMassCat==KKMassCat::' + bin),False),
                                    Normalization =  self._CPnormFracs[bin][comp] * self._sliceNormFracs[bin] )
                if not binInd==self._nKKbins-1:addPdfOpt_i.update(dict( Invisible = ()                       ))
                if     binInd==self._nKKbins-1:addPdfOpt_i.update(dict( LineColor = self._lineColors[comp],
                                                                        LineStyle = self._lineStyles[comp],
                                                                        LineWidth = self._lineWidth          ))
                if not binInd==0: #odd   Components First index = len(self._binNames)
                                  #swave Components First index = 2 *( len(self._binNames)
                                  if comp=='even' : argAddTo = ( 'addPDF{0}'.format(binInd-1),1.,1.)
                                  if comp=='odd'  : argAddTo = ( 'addPDF{0}'.format(  len(self._binNames) + binInd-1 ,1.,1.) )
                                  if comp=='swave': argAddTo = ( 'addPDF{0}'.format(2*len(self._binNames) + binInd-1 ,1.,1.) )
                                  addPdfOpt_i.update(dict(AddTo = argAddTo))
                opts.append(addPdfOpt_i)
        return opts

    def getPdfOptsSixKKbins(self, BinData=True, bins=20):
        projecting = self.getProJWdata(BinData,bins)
        projData   = projecting['data']
        projVars   = projecting['vars']
        if not BinData:  projVars.remove( self._tPdf.indexCat() )
        opts = {}
        KKCat = 'KKMassCat==KKMassCat::'
        for b in self._binNames:
            opts.update( {b : dict(  ProjWData = (projData.reduce(KKCat+b).reduce(ArgSet=projVars), False)
                                   , LineWidth = self._lineWidth
                                   , LineStyle = self._lineStyles['total']
                                   , LineColor = self._lineColors['total'])
                          }  )
        return opts

    def getAddPdfsOptsSixKKbins(self,BinData=True,bins=20):
        if not self._CPnormFracs:    self.calculateCPnormFracs()
        projecting = self.getProJWdata(BinData,bins)
        projData   = projecting['data']
        projVars   = projecting['vars']
        if not BinData:  projVars.remove( self._tPdf.indexCat() )
        opts = []
        for bin in self._binNames:
            ith_binOpts = { }
            for comp in self._comps:
                opt = dict(  ProjWData     = (projData.reduce('KKMassCat==KKMassCat::'+bin).reduce(ArgSet=projVars),False)
                           , LineColor     =  self._lineColors[comp]
                           , LineStyle     =  self._lineStyles[comp]
                           , LineWidth     =  self._lineWidth
                           , Normalization =  self._CPnormFracs[bin][comp]
                             )
                ith_binOpts.update( {comp:opt} )
            opts.append( ith_binOpts  )
        return opts

    def setLineColors(self,colors): self._lineColors = colors
    def setLineStyles(self,styles): self._lineStyles = styles
    def setLineWidth(self, width ): self._lineWidth  = width


# function for plotting the S-wave parameters versus the (binned) KK mass
def plotSWaveBins( **kwargs ) :
    mode = kwargs.pop( 'Mode', 'phases' )
    if mode not in [ 'phases', 'fractions', 'events' ] :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: possible plot modes: "phases", "fractions", "events"'
    if any( key not in kwargs for key in [ 'SValues', 'SLowErrors', 'SHighErrors' ] ) :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: "SValues", "SLowErrors" and "SHighErrors" arguments are required'

    defRange  = ( None, None ) if mode == 'phases' else ( 0., None )
    defSLabel = '#delta_{S} - #delta_{#perp}    [rad]' if mode == 'phases' else\
                'F_{S}' if mode == 'fractions' else\
                'N_{S} / (MeV/c^{2})'
    yAxisRange  = kwargs.pop( 'SAxisRange',    defRange                      )
    KKMassLabel = kwargs.pop( 'KKMassLabel',   'm(K^{+}K^{-}) [MeV/c^{2}]'   )
    SLabel      = kwargs.pop( 'SLabel',        defSLabel                     )
    plotTitle   = kwargs.pop( 'PlotTitle',     ''                            )
    LHCbText1   = kwargs.pop( 'LHCbTextLine1', 'LHCb'                        )
    LHCbText2   = kwargs.pop( 'LHCbTextLine2', ''                            )
    drawLegend  = kwargs.pop( 'DrawLegend',    False                         )
    massBins    = kwargs.pop( 'MassBins',      [ 988., 1008., 1032., 1050. ] )
    theoryVals  = kwargs.pop( 'TheoryValues',  None                          )
    SVals       = kwargs.pop( 'SValues'                                      )
    SLowErrs    = kwargs.pop( 'SLowErrors'                                   )
    SHighErrs   = kwargs.pop( 'SHighErrors'                                  )
    gray        = kwargs.pop( 'GrayScale',     False                         )

    if kwargs :
        raise KeyError, 'P2VV - ERROR: plotSWaveBins: unexpected keyword arguments: %s' % kwargs

    from array import array
    KKMass         = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMassErr      = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] )                for it in range( len(massBins) - 1 ) ] )

    offs = 0.35 if mode == 'phases' else 0.
    KKMass1        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass1LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs                for it in range( len(massBins) - 1 ) ] )
    KKMass1HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs                for it in range( len(massBins) - 1 ) ] )

    KKMass2        = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs + massBins[it] for it in range( len(massBins) - 1 ) ] )
    KKMass2LowErr  = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) - offs                for it in range( len(massBins) - 1 ) ] )
    KKMass2HighErr = array( 'd', [ 0.5 * ( massBins[it + 1] - massBins[it] ) + offs                for it in range( len(massBins) - 1 ) ] )

    from ROOT import TGraphAsymmErrors
    SGraphs = [ ]

    S1        = array( 'd', [ SVals[it]     / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SVals     )
    S1LowErr  = array( 'd', [ SLowErrs[it]  / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SLowErrs  )
    S1HighErr = array( 'd', [ SHighErrs[it] / (massBins[it+1] - massBins[it]) for it in range(len(massBins)-1) ] if mode == 'events' else SHighErrs )
    SGraphs.append( TGraphAsymmErrors( len(KKMass1), KKMass1, S1, KKMass1LowErr, KKMass1HighErr, S1LowErr, S1HighErr ) )

    if mode == 'phases' :
        from math import pi
        S2        = array( 'd', [ pi - val for val in S1 ] )
        S2LowErr  = array( 'd', S1HighErr )
        S2HighErr = array( 'd', S1LowErr  )
        SGraphs.append( TGraphAsymmErrors( len(KKMass2), KKMass2, S2, KKMass2LowErr, KKMass2HighErr, S2LowErr, S2HighErr ) )

    theory = None
    if theoryVals :
        from ROOT import TGraphErrors
        theory    = array( 'd', theoryVals           )
        theoryErr = array( 'd', [ 0. ] * len(KKMass) )
        SGraphs.append( TGraphErrors( len(KKMass), KKMass, theory, KKMassErr, theoryErr ) )

    SMin0 = min( val for val in theory ) if theory else +1.e32
    SMax0 = max( val for val in theory ) if theory else -1.e32
    SMin1 = min( val - err for val, err in zip( S1, S1LowErr  ) )
    SMax1 = max( val + err for val, err in zip( S1, S1HighErr ) )
    SMin  = min( [ SMin0, SMin1, ( pi - SMax1 ) if mode == 'phases' else +1.e32 ] )
    SMax  = max( [ SMax0, SMax1, ( pi - SMin1 ) if mode == 'phases' else -1.e32 ] )

    from ROOT import gStyle, kSolid
    gStyle.SetEndErrorSize(3)
    gStyle.SetLineStyleString( 11, ' 30 15' )
    SGraphs[0].SetLineStyle(kSolid)
    if len(SGraphs) > 1 : SGraphs[1].SetLineStyle(kSolid)
    if len(SGraphs) > 2 : SGraphs[2].SetLineStyle(11)

    from ROOT import kBlack, kBlue, kRed
    SGraphs[0].SetLineColor( kBlue + 2 )
    if len(SGraphs) > 1 : SGraphs[1].SetLineColor( kRed - 4 )
    if len(SGraphs) > 2 : SGraphs[2].SetLineColor(kBlack)

    SGraphs[0].SetMarkerColor( kBlue + 2 )
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerColor( kRed - 4 )
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerColor(kBlack)

    SGraphs[0].SetLineWidth(3)
    if len(SGraphs) > 1 : SGraphs[1].SetLineWidth(3)
    if len(SGraphs) > 2 : SGraphs[2].SetLineWidth(3)

    from ROOT import kFullCircle, kFullTriangleDown
    SGraphs[0].SetMarkerStyle(kFullCircle)
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerStyle(kFullTriangleDown)
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerStyle(kFullCircle)
    SGraphs[0].SetMarkerSize(1.2)
    if len(SGraphs) > 1 : SGraphs[1].SetMarkerSize(1.4)
    if len(SGraphs) > 2 : SGraphs[2].SetMarkerSize(0.7)

    for graph in SGraphs :
        graph.SetMinimum( yAxisRange[0] if yAxisRange[0] != None else SMin - 0.10 * ( SMax - SMin ) )
        graph.SetMaximum( yAxisRange[1] if yAxisRange[1] != None else SMax + 0.15 * ( SMax - SMin ) )

        graph.GetXaxis().SetTitle(KKMassLabel)
        graph.GetYaxis().SetTitle(SLabel)

        graph.GetXaxis().SetTitleOffset(1.1)
        graph.GetYaxis().SetTitleOffset( 0.8 if mode == 'phases' else 1.1 if mode == 'fractions' else 0.9 )

        graph.SetTitle(plotTitle)

        _P2VVPlotStash.append(graph)

    if mode == 'phases' and drawLegend :
        from ROOT import gStyle, TLegend
        leg = TLegend( 0.65, 0.46, 0.91, 0.66 )
        leg.SetTextFont( gStyle.GetTextFont() )
        leg.SetTextSize(0.072)
        leg.SetMargin(0.45)
        leg.AddEntry( SGraphs[0], '#Delta#Gamma_{s} > 0', 'LPE' )
        leg.AddEntry( SGraphs[1], '#Delta#Gamma_{s} < 0', 'LPE' )
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        _P2VVPlotStash.append(leg)
    else :
        leg = None

    if LHCbText1 or LHCbText2 :
        from ROOT import TPaveText
        LHCbText = TPaveText( 0.24, 0.73 if LHCbText1 and LHCbText2 else 0.81, 0.70 if LHCbText1 and LHCbText2 else 0.37, 0.89, 'BRNDC' )
        if LHCbText1 : LHCbText.AddText(LHCbText1)
        if LHCbText2 : LHCbText.AddText(LHCbText2)
        LHCbText.SetShadowColor(0)
        LHCbText.SetFillStyle(0)
        LHCbText.SetBorderSize(0)
        LHCbText.SetTextAlign(12)
        LHCbText.SetTextSize(0.072)
        _P2VVPlotStash.append(LHCbText)
    else :
        LHCbText = None

    from ROOT import TCanvas
    canv = TCanvas( 'SWavePhaseCanv' if mode == 'phases' else 'SWaveFracCanv' if mode == 'fractions' else 'SWaveEventsCanv', 'S-Wave' )
    canv.SetLeftMargin(0.18)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.18)
    canv.SetTopMargin(0.05)
    if gray :
        canv.SetGrayscale()
    if len(SGraphs) > 1 :
        SGraphs[1].Draw('AP')
        SGraphs[0].Draw('P sames')
    else :
        SGraphs[0].Draw('AP')
    if len(SGraphs) > 2 :
        SGraphs[2].Draw('P sames')
    if leg :
        leg.Draw()
    if LHCbText :
        LHCbText.Draw()

    return canv

def splot( pdf, sdata ) :
    # switch off all yields, except current one
    from contextlib import contextmanager
    @contextmanager
    def __select_component( i, yields ):
        orig = dict( (j,j.getVal()) for j in yields )
        [ j.setVal(0) for j in orig.iterkeys() if j!=i ]
        try     : yield
        finally : [ j.setVal(v) for (j,v) in orig.iteritems() ]
    from ROOT import TCanvas, kDashed, kRed, kGreen, kBlue, kBlack
    canvas = TCanvas(pdf.GetName() + '_splot')
    obs = [ o for o in pdf.Observables() if hasattr(o,'frame') and o not in sdata.usedObservables() ]
    for (p,o) in zip( canvas.pads(len(obs)), obs ) :
        # select yields
        _yields = [ y for y in pdf.Parameters() if y.getAttribute('Yield') ]
        # loop over components
        for (pp,i) in zip( p.pads(1,len(_yields)), _yields ) :
            # switch off all yields, except current one
            with __select_component( i, _yields ) :
                # plot both weighed data and PDF
                # TODO: add the same color coding as above...
                c_name = i.GetName()[2:]
                c_opts = { 'signal'             : dict( LineColor = kGreen )
                         , 'psi_background'     : dict( LineColor = kRed )
                         , 'cmb_background'     : dict( LineColor = kBlue )
                         }
                from P2VV.Utilities.Plotting import plot
                plot( pp, o, sdata.data( c_name ), pdf, pdfOpts = c_opts[c_name] if c_name in c_opts else {})
    return canvas

def Oldsplot( canv, sdata, pdf, frameOpts = dict(), dataOpts = dict() , pdfOpts = dict() ) :
    obs = [ o for o in pdf.Observables() if hasattr(o,'frame') and o not in sdata.usedObservables() ]
    for (p,o) in zip( canv.pads(len(obs)), obs ) :
        # snapshot yeilds
        _yields = dict( (p,p.getVal()) for p in pdf.Parameters() if p.getAttribute('Yield')  )
        # loop over components
        for (pp,i) in zip( p.pads(1,len(_yields)), _yields.iterkeys() ) :
            # switch off all yields, except current one
            for j in filter( lambda x: x!=i, _yields.iterkeys() ) : j.setVal(0)
            # and plot both weighed data and PDF
            from P2VV.Utilities.Plotting import plot
            c_name = i.GetName()[2:]
            if 'Title' not in frameOpts : frameOpts['Title'] =  '%s : %s' % ( c_name , o.GetTitle() )
            plot( pp, o, sdata.data( c_name ), pdf
                , pdfOpts = pdfOpts[c_name] if c_name in pdfOpts else {}
                , frameOpts = frameOpts
                , dataOpts = dataOpts
                )
            # and put back the original value!
            for (i,v) in _yields.iteritems() : i.setVal(v)

#To help sort the pads in the plot function
class Sorter(object):
    def __init__(self, d):
        self.__d = d
    def __call__(self, o):
        if o in self.__d:  return self.__d[o]
        else:              return len(self.__d) + 1

# Save a plot/Canvas in a root file
def saveInRootFile(plot, fileName=''):
    from ROOT import TFile
    if '.C' in fileName: name = fileName
    elif fileName: name = fileName + '.C'
    else:          name = plot.GetName() + '.C'
    
    f = TFile.Open(name,'recreate')
    plot.Write()
    f.Close()
    print 'P2VV -INFO:saveInRootFile: Wrote %s to file: %s.'%(plot.GetName(),name)
