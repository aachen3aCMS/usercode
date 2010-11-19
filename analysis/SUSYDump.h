#include <iostream>
#include <iomanip>

void SUSYAna::BasicDump(int i) {

  cout << setprecision(3);

  cout << endl;
  cout << "  ======================================================================" << endl;
  cout << endl;
  cout << "         Tree/Chain Entry " << setw(13) << i << endl;
  cout << endl;
  cout << "   Run " << setw(11) << global_run << "     Event " << setw(10) << global_event 
       << "     Lumi Section " << setw(8) << lumi_section << endl;
  cout << "   Store " << setw(9) << global_store << "     Orbit " << setw(10) << global_orbit 
       << "     BX " << setw(18) << global_bx << endl;
  cout << endl;
  cout << "   B Field " << setw(7) << global_bfield << "     Trigger List         " 
       << unpack(trig_HLTName) << endl;
  cout << endl;
  cout << "   Hottest ECAL Cell R9 [ pT eta phi ] [ time chi flag ] : " << setw(7) 
       << noise_ecal_r9 << " [" << setw(7) << fixed << noise_ecal_pt << setw(7) << noise_ecal_eta 
       << setw(7) << noise_ecal_phi << " ] [" << setw(7) << noise_ecal_time << setw(7) 
       << noise_ecal_chi << setw(3) << noise_ecal_flag << " ] " << endl;
  cout << endl;
  if (!global_isdata) {
    cout << "   PDF Information" << endl;
    cout << "   No   pdgid    x         pthat : " << setw(7) << pdf_scale << endl;
    cout << "    1 " << setw(6) << fixed << pdf_id1 << setw(8) << pdf_x1 << endl;
    cout << "    2 " << setw(6) << pdf_id2 << setw(8) << pdf_x2 << endl;
    cout << endl;
  }
  cout << "  ======================================================================" << endl;
  cout << endl;
}

void SUSYAna::TriggerDump(TString sel) {

  cout << setprecision(3);
  
  cout << "  ===  Trigger Dump - #Objects : " << trig_n << " - Selection : " << sel << "  === " << endl;

  if (trig_n<1) {
    cout << endl;
    return;
  }

  cout << "   No  Trigger Name                   Filter Name                                  Prescale L1  HLT       pT      eta      phi " << endl;

  for (int i=0; i<trig_n; i++) {
    TString tname =  unpack(trig_name[i]);
    if (sel != "*" && !tname.Contains(sel))
      continue;
    cout.flags(ios::right);
    cout << setw(5) << i << "  ";
    cout.flags(ios::left);
    cout << setw(30) << tname << " " << setw(50) << unpack(trig_filter[i]);
    cout.flags(ios::right);
    cout << " " << setw(5) << fixed << trig_L1prescale[i] << setw(5) << trig_HLTprescale[i] 
	 << setw(9) << trig_pt[i] << setw(9) << trig_eta[i] << setw(9) << trig_phi[i] << endl;
  }
  cout << endl;

}

void SUSYAna::MuonDump(bool full) {

  cout << setprecision(3);

  cout << "  ===  Muon Dump - #Objects : " << muo_n << "  === " << endl;

  if (muo_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi      Isolation                       CM                      TK                   ID      ";
  if (full)    
    cout << "      Valid Hits                  Trigger Matches" << endl;
  else      cout << endl;
  cout << "                                    Trk   ECal   HCal    Rel       chi hits      d0          ~                GBL STA TK";
  if (full) 
    cout << "    PXL TK MUO CHAM MEAS NOTMEAS " << endl;
  else
    cout << endl;
  for (int i=0; i<muo_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << i << setw(9) << fixed << muo_charge[i]*muo_pt[i] << setw(8) 
	 << muo_eta[i] << setw(8) << muo_phi[i];
    cout << setprecision(3) << fixed;
    cout << " [" << setw(7) << muo_TrkIso[i] << setw(7) << muo_ECalIso[i] << setw(7) << muo_HCalIso[i] 
	 << setw(7) << (muo_TrkIso[i]+muo_ECalIso[i]+muo_HCalIso[i])/muo_pt[i] << " ]";
    cout << " [" << setw(7) << muo_TrkChiNormCm[i] << setw(4) << muo_hitsCm[i] << setw(8) << muo_d0Cm[i] << " ]";
    cout << " [" << setw(7) << muo_TrkChiNormTk[i] << setw(4) << muo_hitsTk[i] << setw(8) << muo_d0Tk[i] << " ]";
    cout << " [  " << muo_ID[i][1] << "  " << muo_ID[i][2] << "  " << muo_ID[i][3] << "  ]";
    if (full) {
      cout << " [" << setw(3) << muo_ValidPixelHitsCm[i] << setw(3) << muo_ValidTrackerHitsCm[i] << setw(3) 
	   << muo_ValidMuonHitsCm[i] << setw(6) << muo_ChambersMatched[i] << setw(5) 
	   << muo_TrackerLayersMeasCm[i] << setw(5) << muo_TrackerLayersNotMeasCm[i] << " ]";
      cout << " [";
      for (int k=0; k<muo_trign[i]; k++) 
	cout << " " << muo_trig[i][k];
      cout << " ]" <<endl;
    }
    else
      cout << endl;
  }
  if (full) {
    cout << endl;
    cout << "   Cocktail" << endl;
    for (int i=0; i<muo_n; i++) {
      if (muo_Cocktail_pt[i]>-1) {
	cout.flags(ios::right);
	cout << setw(5) << i << setw(9) << fixed << muo_Cocktail_pt[i] << setw(8) 
	     << muo_Cocktail_eta[i] << setw(8) << muo_Cocktail_phi[i] << endl;
      }
    }
  }
  cout << endl;

}

void SUSYAna::CaloJetDump() {

  cout << setprecision(3);

  cout << "  ===  CaloJet Dump - #Objects : " << calojet_n << "  === " << endl;

  if (calojet_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi     fem  fHPD  fRBX n90 const     btag ID truth" << endl;
  for (int i=0; i<calojet_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << i << setw(9) << fixed << calojet_pt[i] << setw(8) << calojet_eta[i] 
	 << setw(8) << calojet_phi[i] << fixed << setprecision(3) << " [" << setw(6) << calojet_fem[i] 
	 << setw(6) << calojet_fHPD[i] << setw(6) << calojet_fRBX[i] 
	 << setw(4) << calojet_n90[i] << setw(4) << calojet_const[i] << " ]";
    cout << setw(9) << fixed << calojet_btag[i] << setw(3) << calojet_ID[i] << setw(4) 
	 << calojet_truth[i] << endl;
  }
  cout << endl;

}

void SUSYAna::PFJetDump() {

  cout << setprecision(3);

  cout << "  ===  PFJet Dump - #Objects : " << pfjet_n << "  === " << endl;

  if (pfjet_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi   n90 const     Energy fractions                               Multiplicities           btag truth" << endl;
  cout << "                                             c_had n_had  n_em  c_em    mu    el gamma              ~" << endl;
  for (int i=0; i<pfjet_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << i << setw(9) << fixed << pfjet_pt[i] << setw(8) << pfjet_eta[i] 
	 << setw(8) << pfjet_phi[i] << fixed << setprecision(3) << " [" << setw(4) 
	 << pfjet_n90[i] << setw(4) << pfjet_const[i] << " ]";
    cout << " [";
    cout << setprecision(3);
    for (int k=0; k<7; k++) 
      cout << setw(6) << fixed << pfjet_PFF[i][k];
    cout << " ] [";
    for (int k=0; k<7; k++) 
      cout << setw(3) << fixed << pfjet_PFN[i][k];
    cout << " ] " << setw(9) << fixed << pfjet_btag[i] << setw(4) << pfjet_truth[i] << endl;

  }
  cout << endl;

}

void SUSYAna::TruthJetDump() {

  cout << setprecision(3);

  cout << "  ===  TruthJet Dump - #Objects : " << truthjet_n << "  === " << endl;

  if (truthjet_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi  " << endl;
  for (int i=0; i<truthjet_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << fixed << i << setw(9) << truthjet_pt[i] << setw(8) << truthjet_eta[i] 
	 << setw(8) << truthjet_phi[i] << endl;
  }
  cout << endl;

}

void SUSYAna::TruthDump() {
  
  cout << setprecision(3);

  cout << "  ===  Truth Dump - #Objects : " << truth_n << " with " << truthl_n 
       << " final truth leptons  === " << endl;
  
  if (truth_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi  pdgid    beg -> end     truthl" << endl;
  for (int i=0; i<truth_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << fixed << i << setw(9) << truth_pt[i] << setw(8) << truth_eta[i] 
	 << setw(8) << truth_phi[i] << setw(7) << truth_pdgid[i];
    cout << "  [" << setw(5) << truth_bvtxid[i] << setw(5) << truth_evtxid[i] << " ]";
    for (int k=0; k<truthl_n; k++) {
      if (truthl_ori[k]==i) {
	cout << "  -->  " << k;
      }
    }
    cout << endl;
  }
  cout << endl;
  for (int i=0; i<truthl_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << fixed << i << setw(9) << truthl_pt[i] << setw(8) << truthl_eta[i] 
	 << setw(8) << truthl_phi[i] << setw(7) << truthl_pdgid[i] << endl;;
    
  }
  cout << endl;

}

void SUSYAna::VertexDump() {

  cout << setprecision(3);

  cout << "  ===  Vertex Dump - #Objects : " << vtx_n << "  === " << endl;
  if (vtx_n>0) {
    cout << "   No        x        y        z fake      chi      ndf   ntr" << endl;
    for (int i=0; i<vtx_n; i++) {
      cout.flags(ios::right);
      cout << setw(5) << fixed << i << setw(9) << vtx_x[i] << setw(9) << vtx_y[i] 
	   << setw(9) << vtx_z[i] << setw(5) << vtx_fake[i] << setw(9) << vtx_chi[i] 
	   << setw(9) << vtx_ndof[i] << setw(6) << vtx_ntr[i] << endl;  
    }
  }
  cout << endl;
  cout << "   BS " << setw(8) << fixed << bs_x << setw(9) << bs_y << setw(9) << bs_z;  
  cout << "   Tracks " << setw(4) << tracks_n << "    hqf " << setw(7) << tracks_hqf << endl;
  cout << endl;

}

void SUSYAna::METDump() {

  cout << setprecision(3);

  cout << "  ===  MET Dump  === " << endl;
  cout << "  Type          MET     phi       MEx     MEy        SumEt   SumEtSig" << endl;
  cout << "  CaloMET " << setw(9) << fixed << met_et[0] << setw(8) << met_phi[0] 
       << "  [" << setw(8) << met_ex[0] << setw(8) << met_ey[0] 
       << " ] [" << setw(9) << met_sumet[0] << setw(9) << met_sumetsig[0] << " ]" << endl;
  cout << "  PFMET   " << setw(9) << fixed << met_et[3] << setw(8) << met_phi[3] 
       << "  [" << setw(8) << met_ex[3] << setw(8) << met_ey[3] 
       << " ] [" << setw(9) << met_sumet[3] << setw(9) << met_sumetsig[3] << " ]" << endl;
  cout << "  TCMET   " << setw(9) << fixed << met_et[4] << setw(8) << met_phi[4] 
       << "  [" << setw(8) << met_ex[4] << setw(8) << met_ey[4] 
       << " ] [" << setw(9) << met_sumet[4] << setw(9) << met_sumetsig[4] << " ]" << endl;
  cout << "  GENMET1 " << setw(9) << fixed << met_et[1] << setw(8) << met_phi[1] 
       << "  [" << setw(8) << met_ex[1] << setw(8) << met_ey[1] 
       << " ] [" << setw(9) << met_sumet[1] << setw(9) << met_sumetsig[1] << " ]" << endl;
  cout << "  GENMET2 " << setw(9) << fixed << met_et[2] << setw(8) << met_phi[2] 
       << "  [" << setw(8) << met_ex[2] << setw(8) << met_ey[2] 
       << " ] [" << setw(9) << met_sumet[2] << setw(9) << met_sumetsig[2] << " ]" << endl;
  cout << endl;

}

void SUSYAna::SCDump() {

  cout << setprecision(3);

  cout << "  ===  SuperCluster Dump - #Objects : " << SC_n << "  === " << endl;

  if (SC_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       E      eta     phi  truth " << endl;
  for (int i=0; i<SC_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << fixed << i << setw(9) << SC_E[i] << setw(8) << SC_eta[i] 
	 << setw(8) << SC_phi[i] << setw(5) << SC_truth[i] << endl;    
  }
  cout << endl;

}

void SUSYAna::EleDump(bool full) {

  cout << setprecision(3);

  cout << "  ===  Electron Dump - #Objects : " << ele_n << "  === " << endl;

  if (ele_n<1) {
    cout << endl;
    return;
  }

  cout << "   No       pT     eta     phi     Isolation (0.3)      ECAL TK       chi hits      d0      fbrem  convr      HoE   SC truth";
  if (full) {
    cout << "   PXL EXP            Electron ID                     Trigger Matches" << endl;
  }
  else
    cout << endl;
  cout << "                                    Trk   ECal   HCal";
  if (full)
    cout << setw(119) << "deta   dphi  HCalDepth1/2   E5x5" << endl;
  else
    cout << endl;
  for (int i=0; i<ele_n; i++) {
    cout.flags(ios::right);
    cout << setw(5) << i << setw(9) << fixed << ele_charge[i]*ele_pt[i] << setw(8) 
	 << ele_eta[i] << setw(8) << ele_phi[i];
    cout << setprecision(3) << fixed;
    cout << " [" << setw(7) << ele_Dr03TkSumPt[i] << setw(7) << ele_Dr03ECalSumEt[i] 
	 << setw(7) << ele_Dr03HCalSumEt[i] << " ] [";
    cout << setw(3) << ele_isECal[i] << setw(3) << ele_isTracker[i] << " ] [";
    cout << setw(7) << ele_TrkChiNorm[i] << setw(4) << ele_hits[i] << setw(8) 
	 << ele_d0vtx[i] << " ] [";
    cout << setw(7) << ele_fbrem[i] << setw(7) << ele_convr[i] << " ]";
    cout << setw(7) << ele_HCalOverEm[i] << setw(5) << ele_SC[i] << setw(5) << ele_truth[i];
    if (full) {
      cout << "   [" << setw(3) << ele_ValidHitFirstPxlB[i] << setw(3) 
	   << ele_TrkExpHitsInner[i] << " ] [";
      cout << setw(7) << ele_dEtaSCTrackAtVtx[i] << setw(7) << ele_dPhiSCTrackAtVtx[i] 
	   << setw(7) << ele_dr03HcalDepth1[i] << setw(7) << ele_dr03HcalDepth2[i]
	   << setw(7) << ele_e5x5[i] << " ]";
      cout << " [";
      for (int k=0; k<ele_trign[i]; k++) 
	cout << " " << ele_trig[i][k];
      cout << " ]" <<endl;
    }
    else
      cout << endl;
  }
  cout << endl;

}
