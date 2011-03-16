from ROOT import *
from math import pi,sqrt

RDSfile = TFile('MergedRDS.root')
RDS = RDSfile.Get('dlogLLDataSet')

numentries = RDS.numEntries()
NbinsX = int(sqrt(numentries))
NbinsY = int(sqrt(numentries))

param1_min = 0
param1_max = 2*pi

param2_min = -0.7
param2_max = 0.7

ProfileLikelihood = TH2D('ProfileLikelihood','ProfileLikelihood',NbinsX,param1_min,param1_max,NbinsY,param2_min,param2_max)
x = ProfileLikelihood.GetXaxis()
y = ProfileLikelihood.GetYaxis()

for k in range(numentries):
    event = RDS.get(k)
    i = int(event.getRealValue('stepphis'))
    print 'i = ', i
    j = int(event.getRealValue('stepdeltaGamma'))
    print 'j = ', j

    BaseNLL = event.getRealValue('BaseNLL')
    GridpointNLL = event.getRealValue('GridPointNLL')
    dLL = GridpointNLL - BaseNLL
    ProfileLikelihood.SetBinContent(i,j,dLL)

    print '######################'
    print 'Check'
    print '######################'
    print 'Bincenter from RDS = ',event.getRealValue('valuephis')
    print 'GetBinCenter i =', x.GetBinCenter(i)
    print 'Bincenter from RDS = ',event.getRealValue('valuedeltaGamma')
    print 'GetBinCenter j =', y.GetBinCenter(j)
    print 'dLL from RDS = ',dLL
    print 'Value at (i,j) =', ProfileLikelihood.GetBinContent(i,j)


tfile = TFile('UntaggedProfileLikelihood.root','RECREATE')
ProfileLikelihood.Write() 
tfile.Close()
