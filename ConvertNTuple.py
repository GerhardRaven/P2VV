#!/usr/bin/env python
from ROOT import TTree, TFile, AddressOf, gROOT
from math import cos

def main():
  # Make a tree

  oldf = TFile('Bs2JpsiPhi_stripping10_1cand.root')
  oldt = oldf.Get('DecayTree')
  numentries = oldt.GetEntries()
  print 'numentries =',numentries
  newf = TFile('reduced_ntuple_NoJpsiMassCut.root','RECREATE')
  newt = TTree('MyTree','MyTree')
  gROOT.ProcessLine(\
    "struct MyStruct{\
    Double_t Double_m;\
    Double_t Double_mdau1;\
    Double_t Double_mdau2;\
    Int_t Int_pid;\
    Int_t Int_piddau1;\
    Int_t Int_piddau2;\
    Double_t Double_t;\
    Double_t Double_sigmat;\
    Int_t Int_tagdecision;\
    Double_t Double_tagomega;\
    Double_t Double_trcospsi;\
    Double_t Double_trcostheta;\
    Double_t Double_trphi;\
    Double_t Double_helcosthetaK;\
    Double_t Double_helcosthetaL;\
    Double_t Double_helphi;\
    };")

  from ROOT import MyStruct

  # Create branches in the tree
  s = MyStruct()
  newt.Branch('m',AddressOf(s,'Double_m'),'Double_m/D')
  newt.Branch('mdau1',AddressOf(s,'Double_mdau1'),'Double_mdau1/D')
  newt.Branch('mdau2',AddressOf(s,'Double_mdau2'),'Double_mdau2/D')
  
  newt.Branch('pid',AddressOf(s,'Int_pid'),'Int_pid/I')
  newt.Branch('piddau1',AddressOf(s,'Int_piddau1'),'Int_piddau1/I')
  newt.Branch('piddau2',AddressOf(s,'Int_piddau2'),'Int_piddau2/I')

  newt.Branch('t',AddressOf(s,'Double_t'),'Double_t/D')
  newt.Branch('sigmat',AddressOf(s,'Double_sigmat'),'Double_sigmat/D')

  newt.Branch('tagdecision',AddressOf(s,'Int_tagdecision'),'Int_tagdecision/I')
  newt.Branch('tagomega',AddressOf(s,'Double_tagomega'),'Double_tagomega/D')

  newt.Branch('trcospsi',AddressOf(s,'Double_trcospsi'),'Double_trcospsi/D')
  newt.Branch('trcostheta',AddressOf(s,'Double_trcostheta'),'Double_trcostheta/D')
  newt.Branch('trphi',AddressOf(s,'Double_trphi'),'Double_trphi/D')

  newt.Branch('helcosthetaK',AddressOf(s,'Double_helcosthetaK'),'Double_helcosthetaK/D')
  newt.Branch('helcosthetaL',AddressOf(s,'Double_helcosthetaL'),'Double_helcosthetaL/D')
  newt.Branch('helphi',AddressOf(s,'Double_helphi'),'Double_helphi/D')


  c = 0.299792458
  # Fill tree
  for i in range(numentries):
      oldt.GetEntry(i)

      s.Double_m = oldt.B_s0_LOKI_MASS_JpsiConstr
      s.Double_mdau1 = oldt.J_psi_1S_MM
      s.Double_mdau2 = oldt.phi_1020_MM
      
      s.Int_pid = oldt.B_s0_ID
      s.Int_piddau1 = oldt.J_psi_1S_ID
      s.Int_piddau2 = oldt.phi_1020_ID
      
      s.Double_t = (oldt.B_s0_LOKI_DTF_CTAU)/c#ctau is in mm standard LHCb is mm so c = 2.99792458*10^8 1000mm/10^9ns = 299.792458 mm/ns. We want tau in ps-> 0.299792458 mm/ps
      s.Double_sigmat = (oldt.B_s0_LOKI_DTF_CTAUERR)/c
      
      s.Int_tagdecision = oldt.B_s0_TAGDECISION
      s.Double_tagomega = oldt.B_s0_TAGOMEGA
      
      s.Double_trcospsi = cos(oldt.B_s0_ThetaVtr)
      s.Double_trcostheta = cos(oldt.B_s0_ThetaTr)
      s.Double_trphi = oldt.B_s0_PhiTr
      
      s.Double_helcosthetaK = cos(oldt.B_s0_ThetaK)
      s.Double_helcosthetaL = cos(oldt.B_s0_ThetaL)
      s.Double_helphi = oldt.B_s0_Phi

#      if (oldt.triggeredByUnbiasedHlt1AndHlt2 == 1 and (oldt.B_s0_MINIPCHI2NEXTBEST > 50 or oldt.B_s0_MINIPCHI2NEXTBEST == -1) and oldt.hasBestVtxChi2 == 1 and (abs((oldt.J_psi_1S_MM - 3096.9)/(oldt.J_psi_1S_MMERR)) < 1.4*3)):      
      #if ( ((oldt.Hlt1SingleMuonNoIPL0Decision == 1 or oldt.Hlt1DiMuonNoIPL0DiDecision == 1) and oldt.Hlt2DiMuonUnbiasedJPsiDecision == 1) and (oldt.B_s0_MINIPCHI2NEXTBEST > 50 or oldt.B_s0_MINIPCHI2NEXTBEST == -1) and oldt.hasBestVtxChi2 == 1 and (abs((oldt.J_psi_1S_MM - 3096.9)/(oldt.J_psi_1S_MMERR)) < 1.4*3)):
      if ( ((oldt.Hlt1SingleMuonNoIPL0Decision == 1 or oldt.Hlt1DiMuonNoIPL0DiDecision == 1) and oldt.Hlt2DiMuonUnbiasedJPsiDecision == 1) and (oldt.B_s0_MINIPCHI2NEXTBEST > 50 or oldt.B_s0_MINIPCHI2NEXTBEST == -1) and oldt.hasBestVtxChi2 == 1):
        newt.Fill()

  newf.Write()
  newf.Close()

if __name__=='__main__':main()

