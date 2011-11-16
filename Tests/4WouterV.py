from ROOT import *
wsfile = TFile('4WouterWS.root')
ws = wsfile.Get('ws')
ws.Print()

