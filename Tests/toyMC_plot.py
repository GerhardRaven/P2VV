from ROOT import (RooPullVar, RooConstVar,
                  RooFormulaVar, TFile,
                  TCanvas, RooArgList,
                  RooErrorVar, RooFit)

root_file = TFile.Open('data_100.root')
data = root_file.Get('data')
args = data.get()
tau = args.find('tau')

tau_true = RooConstVar('tau_true', 'tau_true', 2.)
tau_error = RooErrorVar('tau_error', 'tau_error', tau)

_pull = RooFormulaVar('pull', 'pull', '(@0 - @1) / @2', RooArgList(tau, tau_true, tau_error))
pull = data.addColumn(_pull)

from ROOT import TCanvas
canvas = TCanvas('canvas', 'canvas', 1000, 1000)
canvas.Divide(2, 2)
canvas.cd(1)

frame = pull.frame(RooFit.Range(-3, 3))
data.plotOn(frame)
frame.Draw()
