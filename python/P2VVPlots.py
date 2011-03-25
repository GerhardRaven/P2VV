###############################################################################
## P2VVPlots: P2VV plotting tools                                            ##
##                                                                           ##
## authors:                                                                  ##
##   GR,  Gerhard Raven,      Nikhef & VU, Gerhard.Raven@nikhef.nl           ##
##   JvL, Jeroen van Leerdam, Nikhef,      j.van.leerdam@nikhef.nl           ##
##                                                                           ##
###############################################################################


from ROOT import gStyle, gROOT

# set ROOT plot style
plotStyle = gROOT.GetStyle('Plain')

plotStyle.SetTitleBorderSize(0)

plotStyle.SetPadTopMargin(0.10)
plotStyle.SetPadBottomMargin(0.15)
plotStyle.SetPadLeftMargin(0.15)
plotStyle.SetPadRightMargin(0.10)

plotStyle.SetTitleSize(0.048, 'XY')
plotStyle.SetLabelSize(0.045, 'XY')

plotStyle.SetTitleOffset(1.2, 'X')
plotStyle.SetTitleOffset(1.4, 'Y')
plotStyle.SetLabelOffset(0.011, 'XY')

gROOT.SetStyle('Plain')
gROOT.ForceStyle()

gStyle.SetPalette(1)
gStyle.UseCurrentStyle()

def plot(config, obs, canv, data, pdf, components = {}, xTitle = '',
    frameOpts = None, dataOpts = None, pdfOpts = None, logy = False,
    drawRes = False, normalize = True, symmetrize = True) :
  """makes a P2VV plot

  example usage:

  canv = plot(config, 'mpsi', canvas, data, pdf, xTitle = 'M (MeV/c)',
      components = {'psi'    : [RooFit.LineColor(RooFit.kGreen),
                                RooFit.LineStyle(RooFit.kDashed)],
                    'nonpsi' : [RooFit.LineColor(RooFit.kBlue),
                                RooFit.LineStyle(RooFit.kDashed)]},
      frameOpts = [RooFit.Bins(30)],
      dataOpts  = [RooFit.MarkerSize(0.4), RooFit.XErrorSize(0)],
      pdfOpts   = [RooFit.LineWidth(2)])
  """

  from ROOT import RooFit, TLine, TPad
  from P2VVConfiguration import P2VVSetting

  # get observable
  obsVar = config.workspace().function(config[obs].name())

  # get frame options
  frameDrawOpts = frameOpts[:] if frameOpts else None
  P2VVFrameDrawOpts = config.value('P2VVFrameDrawOpts')
  if P2VVFrameDrawOpts :
    if frameDrawOpts :
      for P2VVOpt in P2VVFrameDrawOpts :
        addOpt = True
        for opt in frameDrawOpts :
          if opt.opcode() == P2VVOpt.opcode() :
            addOpt = False
            break

        if addOpt : frameDrawOpts.append(P2VVOpt)

    else :
      frameDrawOpts = P2VVFrameDrawOpts[:]

  # get data options
  dataPlotOpts = dataOpts[:] if dataOpts else None
  P2VVDataPlotOpts = config.value('P2VVDataPlotOpts')
  if P2VVDataPlotOpts :
    if dataPlotOpts :
      for P2VVOpt in P2VVDataPlotOpts :
        addOpt = True
        for opt in dataPlotOpts :
          if opt.opcode() == P2VVOpt.opcode() :
            addOpt = False
            break

        if addOpt : dataPlotOpts.append(P2VVOpt)

    else :
      dataPlotOpts = P2VVDataPlotOpts[:]

  # get PDF options
  pdfPlotOpts = pdfOpts[:] if pdfOpts else None
  P2VVPDFPlotOpts = config.value('P2VVPDFPlotOpts')
  if P2VVPDFPlotOpts :
    if pdfPlotOpts :
      for P2VVOpt in P2VVPDFPlotOpts :
        addOpt = True
        for opt in pdfPlotOpts :
          if opt.opcode() == P2VVOpt.opcode() :
            addOpt = False
            break

        if addOpt : pdfPlotOpts.append(P2VVOpt)

    else :
      pdfPlotOpts = P2VVPDFPlotOpts[:]

  # create frame for observable
  obsFrame = obsVar.frame(*frameDrawOpts) if frameDrawOpts else obsVar.frame()
  config.addPlotObj(obsFrame)

  # set title of x-axis
  if type(xTitle) is str and len(xTitle) > 0 :
    obsFrame.GetXaxis().SetTitle(xTitle)

  # plot data
  if dataPlotOpts : data.plotOn(obsFrame, RooFit.Name('data'), *dataPlotOpts)
  else :            data.plotOn(obsFrame, RooFit.Name('data'))

  # plot PDF
  for comp, opt in components.iteritems() :
    compOpts = opt + pdfPlotOpts if pdfPlotOpts else opt
    pdf.plotOn(obsFrame, RooFit.Components(comp), *compOpts)

  if pdfPlotOpts: pdf.plotOn(obsFrame, RooFit.Name('pdf'), *pdfPlotOpts)
  else :          pdf.plotOn(obsFrame, RooFit.Name('pdf'))

  obsFrame.drawAfter('pdf', 'data')

  #TODO: add chisq/nbins
  #chisq = obsFrame.chiSquare('pdf','data')
  #nbins = obsFrame.GetNbinsX()

  if drawRes :
    # plot residuals
    resHist = obsFrame.residHist('data', 'pdf', normalize)
    config.addPlotObj(resHist)
    resHist.GetXaxis().SetLimits(obsFrame.GetXaxis().GetXmin(),
        obsFrame.GetXaxis().GetXmax())

    resFrame = obsFrame.emptyClone(obsFrame.GetName() + '_resid')
    config.addPlotObj(resFrame)

    if dataPlotOpts :
      for opt in dataPlotOpts : 
        if opt.opcode() == 'MarkerSize' :
          resHist.SetMarkerSize(opt.getDouble(0))
        if opt.opcode() == 'MarkerStyle' :
          resHist.SetMarkerStyle(opt.getInt(0))
        if opt.opcode() == 'MarkerColor' :
          resHist.SetMarkerColor(opt.getInt(0))
        if opt.opcode() == 'Title' :
          resFrame.SetTitle(opt.getString(0))

    resFrame.addPlotable(resHist, 'PZ')

    if symmetrize :
      resMax = max(abs(resHist.getYAxisMin()), abs(resHist.getYAxisMax()))
      resFrame.SetMaximum(resMax)
      resFrame.SetMinimum(-resMax)

    if normalize :
      if resHist.getYAxisMin() > -5 : resFrame.SetMinimum(-5)
      if resHist.getYAxisMax() <  5 : resFrame.SetMaximum(5)

    baseLine = TLine(resFrame.GetXaxis().GetXmin(), 0, resFrame.GetXaxis().GetXmax(), 0)
    baseLine.SetLineColor(RooFit.kRed)
    resFrame.addObject(baseLine)
    # TODO: improve (remove?) axis labels from resFrame,
    # move up against the initial plot

  # draw plots
  canv.cd()
  if drawRes :
    # draw plot (data + PDF)
    plotName = obs + '_plot1'
    plotPad = TPad(plotName, plotName, 0, 0.2, 1, 1)
    config.addPlotObj(plotPad)
    plotPad.SetNumber(1)
    plotPad.Draw()
    canv.cd(1)
    obsFrame.Draw()

    # draw residuals
    canv.cd()
    resName = obs + '_resid1'
    resPad = TPad(resName, resName, 0, 0, 1, 0.2)
    config.addPlotObj(resPad)
    resPad.SetNumber(2)
    resPad.Draw()
    canv.cd(2)
    resFrame.Draw()
  else :
    obsFrame.Draw()

  canv.Update()
  return canv

