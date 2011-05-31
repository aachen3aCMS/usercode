//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct  6 09:55:10 2009 by ROOT version 5.22/00a
// from TTree allData/data after cuts
// found on file: out.root
//////////////////////////////////////////////////////////

#ifndef TreeContent_h
#define TreeContent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#define INTSIZE 4

class TreeContent {
 public :
  
  union u64 {
    char c[INTSIZE];
    int i;
  };
  
  int get_size(const int* a) {
    int size=0;
    u64 tmp;
    for(int i=0; ; ++i ) {
      tmp.i = a[i];
      size++;
      for(int j=0; j < INTSIZE ; ++j )
        if(tmp.c[j] == '\x0')
          return size;
    }
    return -1;
  }
  
  std::string unpack(const int* a) {
    int size = get_size(a);
    u64 tmp;
    char* c = new char[size*INTSIZE];
    for(int i=0; i < size ; ++i) {
      tmp.i = a[i];
      for(int j=0; j < INTSIZE; ++j )
        c[i*INTSIZE+j] = tmp.c[j];
    }
    std::string ret( c );
    delete[] c;
    return ret;
  }

  
  int* pack(const char* c) {
    u64 tmp;
    int j=0,count=0;
    int size = strlen(c)/INTSIZE+1;
    int* ii=new int[size];
    for(int i=0 ; ; i++) {
      tmp.c[j++]=c[i];
      ii[count]=tmp.i;
      if(j==INTSIZE) {
        j=0;
        count++;
      }
      if(c[i] == '\x0')
        break;
    }
    return ii;
  }

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   Double_t        global_weight;
   Int_t           global_procID;
   Double_t        global_pthat;
   Double_t        global_bfield;
   Int_t           global_store;
   Int_t           global_run;
   Int_t           global_event;
   Int_t           global_bx;
   Int_t           global_orbit;
   Int_t           global_exp;
   Int_t           global_isdata;
   Int_t           lumi_section;
   Double_t        lumi_del;
   Double_t        lumi_rec;
   Double_t        lumi_delerr;
   Double_t        lumi_recerr;
   Int_t           pu_bunchx;
   Int_t           pu_n;
   Int_t           pu_num_int[10];   //[pu_n]
   Double_t        pu_inst_Lumi[10][100];   //[pu_n]
   Double_t        pu_zPos[10][100];   //[pu_n]
   Double_t        pu_sumPthi[10][100];   //[pu_n]
   Double_t        pu_sumPtlo[10][100];   //[pu_n]
   Int_t           noise_pLoose;
   Int_t           noise_pTight;
   Int_t           noise_pHigh;
   Double_t        noise_ecal_r9;
   Double_t        noise_ecal_E;
   Double_t        noise_ecal_pt;
   Double_t        noise_ecal_px;
   Double_t        noise_ecal_py;
   Double_t        noise_ecal_pz;
   Double_t        noise_ecal_eta;
   Double_t        noise_ecal_phi;
   Double_t        noise_ecal_time;
   Double_t        noise_ecal_chi;
   Int_t           noise_ecal_flag;
   Int_t           noise_ecal_ieta;
   Int_t           noise_ecal_iphi;
   Double_t        noise_hcal_eventChargeFraction;
   Double_t        noise_hcal_eventEMEnergy;
   Double_t        noise_hcal_eventEMFraction;
   Double_t        noise_hcal_eventHadEnergy;
   Double_t        noise_hcal_eventTrackEnergy;
   Double_t        noise_hcal_flatNoiseSumE;
   Double_t        noise_hcal_flatNoiseSumEt;
   Int_t           noise_hcal_HasBadRBXTS4TS5;
   Double_t        noise_hcal_isolatedNoiseSumE;
   Double_t        noise_hcal_isolatedNoiseSumEt;
   Double_t        noise_hcal_max10GeVHitTime;
   Double_t        noise_hcal_max25GeVHitTime;
   Double_t        noise_hcal_maxE10TS;
   Double_t        noise_hcal_maxE2Over10TS;
   Double_t        noise_hcal_maxE2TS;
   Int_t           noise_hcal_maxHPDHits;
   Int_t           noise_hcal_maxHPDNoOtherHits;
   Int_t           noise_hcal_maxRBXHits;
   Int_t           noise_hcal_maxZeros;
   Double_t        noise_hcal_min10GeVHitTime;
   Double_t        noise_hcal_min25GeVHitTime;
   Double_t        noise_hcal_minE10TS;
   Double_t        noise_hcal_minE2Over10TS;
   Double_t        noise_hcal_minE2TS;
   Double_t        noise_hcal_minHPDEMF;
   Double_t        noise_hcal_minRBXEMF;
   Int_t           noise_hcal_noiseFilterStatus;
   Int_t           noise_hcal_noiseType;
   Int_t           noise_hcal_num10GeVHits;
   Int_t           noise_hcal_num25GeVHits;
   Int_t           noise_hcal_numFlatNoiseChannels;
   Int_t           noise_hcal_numIsolatedNoiseChannels;
   Int_t           noise_hcal_numProblematicRBXs;
   Int_t           noise_hcal_numSpikeNoiseChannels;
   Int_t           noise_hcal_numTriangleNoiseChannels;
   Int_t           noise_hcal_numTS4TS5NoiseChannels;
   Int_t           noise_hcal_passHighLevelNoiseFilter;
   Int_t           noise_hcal_passLooseNoiseFilter;
   Int_t           noise_hcal_passTightNoiseFilter;
   Double_t        noise_hcal_rms10GeVHitTime;
   Double_t        noise_hcal_rms25GeVHitTime;
   Double_t        noise_hcal_spikeNoiseSumE;
   Double_t        noise_hcal_spikeNoiseSumEt;
   Double_t        noise_hcal_triangleNoiseSumE;
   Double_t        noise_hcal_triangleNoiseSumEt;
   Double_t        noise_hcal_TS4TS5NoiseSumE;
   Double_t        noise_hcal_TS4TS5NoiseSumEt;
   Int_t           trig_HLTName[20];
   Int_t           trig_n;
   Int_t           trig_L1prescale[1000];   //[trig_n]
   Int_t           trig_HLTprescale[1000];   //[trig_n]
   Int_t           trig_name[1000][20];   //[trig_n]
   Int_t           trig_filter[1000][20];   //[trig_n]
   Double_t        trig_pt[1000];   //[trig_n]
   Double_t        trig_eta[1000];   //[trig_n]
   Double_t        trig_phi[1000];   //[trig_n]
   Int_t           truth_n;
   Int_t           truth_pdgid[100];   //[truth_n]
   Int_t           truth_bvtxid[100];   //[truth_n]
   Int_t           truth_evtxid[100];   //[truth_n]
   Double_t        truth_E[100];   //[truth_n]
   Double_t        truth_Et[100];   //[truth_n]
   Double_t        truth_p[100];   //[truth_n]
   Double_t        truth_pt[100];   //[truth_n]
   Double_t        truth_px[100];   //[truth_n]
   Double_t        truth_py[100];   //[truth_n]
   Double_t        truth_pz[100];   //[truth_n]
   Double_t        truth_eta[100];   //[truth_n]
   Double_t        truth_phi[100];   //[truth_n]
   Double_t        truth_m[100];   //[truth_n]
   Int_t           truthl_n;
   Int_t           truthl_ori[100];   //[truthl_n]
   Int_t           truthl_pdgid[100];   //[truthl_n]
   Double_t        truthl_E[100];   //[truthl_n]
   Double_t        truthl_Et[100];   //[truthl_n]
   Double_t        truthl_p[100];   //[truthl_n]
   Double_t        truthl_pt[100];   //[truthl_n]
   Double_t        truthl_px[100];   //[truthl_n]
   Double_t        truthl_py[100];   //[truthl_n]
   Double_t        truthl_pz[100];   //[truthl_n]
   Double_t        truthl_eta[100];   //[truthl_n]
   Double_t        truthl_phi[100];   //[truthl_n]
   Int_t           pdf_id1;
   Int_t           pdf_id2;
   Float_t         pdf_x1;
   Float_t         pdf_x2;
   Float_t         pdf_f1;
   Float_t         pdf_f2;
   Float_t         pdf_scale;
   Int_t           vtx_n;
   Int_t           vtx_ntr[100];   //[vtx_n]
   Int_t           vtx_fake[100];   //[vtx_n]
   Double_t        vtx_ndof[100];   //[vtx_n]
   Double_t        vtx_x[100];   //[vtx_n]
   Double_t        vtx_y[100];   //[vtx_n]
   Double_t        vtx_z[100];   //[vtx_n]
   Double_t        vtx_chi[100];   //[vtx_n]
   Double_t        bs_x;
   Double_t        bs_y;
   Double_t        bs_z;
   Int_t           tracks_n;
   Double_t        tracks_hqf;
   Double_t        met_et[10];
   Double_t        met_ex[10];
   Double_t        met_ey[10];
   Double_t        met_eta[10];
   Double_t        met_phi[10];
   Double_t        met_sumet[10];
   Double_t        met_sumetsig[10];
   Double_t        met_etsignif[10];
   Double_t        met_CaloMETInmHF[10];
   Double_t        met_CaloMETInpHF[10];
   Double_t        met_CaloMETPhiInmHF[10];
   Double_t        met_CaloMETPhiInpHF[10];
   Double_t        met_CaloSETInmHF[10];
   Double_t        met_CaloSETInpHF[10];
   Double_t        met_emEtFraction[10];
   Double_t        met_etFractionHadronic[10];
   Double_t        met_maxEtInEmTowers[10];
   Double_t        met_maxEtInHadTowers[10];
   Double_t        met_emEtInHF[10];
   Double_t        met_emEtInEE[10];
   Double_t        met_emEtInEB[10];
   Double_t        met_hadEtInHF[10];
   Double_t        met_hadEtInHE[10];
   Double_t        met_hadEtInHO[10];
   Double_t        met_hadEtInHB[10];
   Double_t        met_ChargedEMEtFraction[10];
   Double_t        met_ChargedHadEtFraction[10];
   Double_t        met_MuonEtFraction[10];
   Double_t        met_NeutralEMFraction[10];
   Double_t        met_NeutralHadEtFraction[10];
   Double_t        met_Type6EtFraction[10];
   Double_t        met_Type7EtFraction[10];
   Int_t           calojet_n;
   Double_t        calojet_E[100];   //[calojet_n]
   Double_t        calojet_Et[100];   //[calojet_n]
   Double_t        calojet_p[100];   //[calojet_n]
   Double_t        calojet_pt[100];   //[calojet_n]
   Double_t        calojet_pt_raw[100];   //[calojet_n]
   Double_t        calojet_px[100];   //[calojet_n]
   Double_t        calojet_py[100];   //[calojet_n]
   Double_t        calojet_pz[100];   //[calojet_n]
   Double_t        calojet_eta[100];   //[calojet_n]
   Double_t        calojet_phi[100];   //[calojet_n]
   Double_t        calojet_fem[100];   //[calojet_n]
   Double_t        calojet_fhad[100];   //[calojet_n]
   Double_t        calojet_btag[100];   //[calojet_n]
   Double_t        calojet_charge[100];   //[calojet_n]
   Double_t        calojet_fHPD[100];   //[calojet_n]
   Double_t        calojet_fRBX[100];   //[calojet_n]
   Int_t           calojet_n90hits[100];   //[calojet_n]
   Int_t           calojet_n90[100];   //[calojet_n]
   Int_t           calojet_flav[100];   //[calojet_n]
   Int_t           calojet_truth[100];   //[calojet_n]
   Int_t           calojet_const[100];   //[calojet_n]
   Int_t           calojet_ID[100];   //[calojet_n]
   Int_t           pfjet_n;
   Double_t        pfjet_E[100];   //[pfjet_n]
   Double_t        pfjet_Et[100];   //[pfjet_n]
   Double_t        pfjet_p[100];   //[pfjet_n]
   Double_t        pfjet_pt[100];   //[pfjet_n]
   Double_t        pfjet_pt_raw[100];   //[pfjet_n]
   Double_t        pfjet_px[100];   //[pfjet_n]
   Double_t        pfjet_py[100];   //[pfjet_n]
   Double_t        pfjet_pz[100];   //[pfjet_n]
   Double_t        pfjet_eta[100];   //[pfjet_n]
   Double_t        pfjet_phi[100];   //[pfjet_n]
   Double_t        pfjet_btag[100];   //[pfjet_n]
   Double_t        pfjet_charge[100];   //[pfjet_n]
   Int_t           pfjet_n90[100];   //[pfjet_n]
   Int_t           pfjet_flav[100];   //[pfjet_n]
   Int_t           pfjet_truth[100];   //[pfjet_n]
   Int_t           pfjet_const[100];   //[pfjet_n]
   Int_t           pfjet_PFN[100][7];   //[pfjet_n]
   Double_t        pfjet_PFF[100][7];   //[pfjet_n]
   Int_t           truthjet_n;
   Double_t        truthjet_E[100];   //[truthjet_n]
   Double_t        truthjet_Et[100];   //[truthjet_n]
   Double_t        truthjet_p[100];   //[truthjet_n]
   Double_t        truthjet_pt[100];   //[truthjet_n]
   Double_t        truthjet_px[100];   //[truthjet_n]
   Double_t        truthjet_py[100];   //[truthjet_n]
   Double_t        truthjet_pz[100];   //[truthjet_n]
   Double_t        truthjet_eta[100];   //[truthjet_n]
   Double_t        truthjet_phi[100];   //[truthjet_n]
   Int_t           fatjet_n;
   Int_t           fatjet_nsub[100];   //[fatjet_n]
   Double_t        fatjet_pt[100];   //[fatjet_n]
   Double_t        fatjet_px[100];   //[fatjet_n]
   Double_t        fatjet_py[100];   //[fatjet_n]
   Double_t        fatjet_pz[100];   //[fatjet_n]
   Double_t        fatjet_E[100];   //[fatjet_n]
   Double_t        fatjet_eta[100];   //[fatjet_n]
   Double_t        fatjet_phi[100];   //[fatjet_n]
   Double_t        fatjet_sub_pt[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_px[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_py[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_pz[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_E[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_eta[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_phi[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fem[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fhad[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_btag[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_n90[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fHPD[100][10];   //[fatjet_n]
   Double_t        fatjet_sub_fRBX[100][10];   //[fatjet_n]
   Int_t           SC_n;
   Int_t           SC_truth[200];   //[SC_n]
   Double_t        SC_E[200];   //[SC_n]
   Double_t        SC_phi[200];   //[SC_n]
   Double_t        SC_eta[200];   //[SC_n]
   Int_t           ele_n;
   Double_t        ele_E[100];   //[ele_n]
   Double_t        ele_Et[100];   //[ele_n]
   Double_t        ele_p[100];   //[ele_n]
   Double_t        ele_pt[100];   //[ele_n]
   Double_t        ele_px[100];   //[ele_n]
   Double_t        ele_py[100];   //[ele_n]
   Double_t        ele_pz[100];   //[ele_n]
   Double_t        ele_eta[100];   //[ele_n]
   Double_t        ele_phi[100];   //[ele_n]
   Double_t        ele_charge[100];   //[ele_n]
   Double_t        ele_TrkChiNorm[100];   //[ele_n]
   Double_t        ele_d0vtx[100];   //[ele_n]
   Double_t        ele_d0bs[100];   //[ele_n]
   Double_t        ele_sd0[100];   //[ele_n]
   Int_t           ele_hits[100];   //[ele_n]
   Int_t           ele_truth[100];   //[ele_n]
   Int_t           ele_isECal[100];   //[ele_n]
   Int_t           ele_isTracker[100];   //[ele_n]
   Int_t           ele_ID[100][100];   //[ele_n]
   Int_t           ele_ValidHitFirstPxlB[100];   //[ele_n]
   Int_t           ele_TrkExpHitsInner[100];   //[ele_n]
   Double_t        ele_HCalOverEm[100];   //[ele_n]
   Double_t        ele_Dr03TkSumPt[100];   //[ele_n]
   Double_t        ele_Dr04HCalSumEt[100];   //[ele_n]
   Double_t        ele_Dr03HCalSumEt[100];   //[ele_n]
   Double_t        ele_Dr04ECalSumEt[100];   //[ele_n]
   Double_t        ele_Dr03ECalSumEt[100];   //[ele_n]
   Double_t        ele_SigmaIetaIeta[100];   //[ele_n]
   Double_t        ele_dEtaSCTrackAtVtx[100];   //[ele_n]
   Double_t        ele_dPhiSCTrackAtVtx[100];   //[ele_n]
   Double_t        ele_dr03HcalDepth1[100];   //[ele_n]
   Double_t        ele_dr03HcalDepth2[100];   //[ele_n]
   Double_t        ele_e2x5Max[100];   //[ele_n]
   Double_t        ele_e5x5[100];   //[ele_n]
   Double_t        ele_e1x5[100];   //[ele_n]
   Double_t        ele_caloEt[100];   //[ele_n]
   Double_t        ele_SCeta[100];   //[ele_n]
   Double_t        ele_convdist[100];   //[ele_n]
   Double_t        ele_convdcot[100];   //[ele_n]
   Double_t        ele_convr[100];   //[ele_n]
   Double_t        ele_fbrem[100];   //[ele_n]
   Int_t           ele_trign[100];   //[ele_n]
   Int_t           ele_trig[100][500];   //[ele_n]
   Int_t           ele_SC[100];   //[ele_n]
   Int_t           ele_numberOfHits[100];   //[ele_n]
   Int_t           pfele_n;
   Double_t        pfele_p[100];   //[pfele_n]
   Double_t        pfele_E[100];   //[pfele_n]
   Double_t        pfele_Et[100];   //[pfele_n]
   Double_t        pfele_pt[100];   //[pfele_n]
   Double_t        pfele_px[100];   //[pfele_n]
   Double_t        pfele_py[100];   //[pfele_n]
   Double_t        pfele_pz[100];   //[pfele_n]
   Double_t        pfele_eta[100];   //[pfele_n]
   Double_t        pfele_phi[100];   //[pfele_n]
   Double_t        pfele_charge[100];   //[pfele_n]
   Int_t           pfele_truth[100];   //[pfele_n]
   Int_t           pfele_trign[100];   //[pfele_n]
   Int_t           pfele_trig[100][500];   //[pfele_n]
   Int_t           pfele_SC[100];   //[pfele_n]
   Int_t           muo_n;
   Double_t        muo_E[100];   //[muo_n]
   Double_t        muo_Et[100];   //[muo_n]
   Double_t        muo_p[100];   //[muo_n]
   Double_t        muo_pt[100];   //[muo_n]
   Double_t        muo_px[100];   //[muo_n]
   Double_t        muo_py[100];   //[muo_n]
   Double_t        muo_pz[100];   //[muo_n]
   Double_t        muo_eta[100];   //[muo_n]
   Double_t        muo_phi[100];   //[muo_n]
   Double_t        muo_charge[100];   //[muo_n]
   Double_t        muo_RelTrkIso[100];   //[muo_n]
   Double_t        muo_TrkIso[100];   //[muo_n]
   Double_t        muo_ECalIso[100];   //[muo_n]
   Double_t        muo_HCalIso[100];   //[muo_n]
   Double_t        muo_TrkIsoDep[100];   //[muo_n]
   Double_t        muo_ECalIsoDep[100];   //[muo_n]
   Double_t        muo_HCalIsoDep[100];   //[muo_n]
   Double_t        muo_AllIso[100];   //[muo_n]
   Double_t        muo_TrkChiNormCm[100];   //[muo_n]
   Double_t        muo_TrkChiNormTk[100];   //[muo_n]
   Double_t        muo_d0Cm[100];   //[muo_n]
   Double_t        muo_d0Tk[100];   //[muo_n]
   Double_t        muo_sd0Cm[100];   //[muo_n]
   Double_t        muo_sd0Tk[100];   //[muo_n]
   Double_t        muo_calocomp[100];   //[muo_n]
   Double_t        muo_calotower_e[100];   //[muo_n]
   Int_t           muo_prompttight[100];   //[muo_n]
   Int_t           muo_hitsCm[100];   //[muo_n]
   Int_t           muo_hitsTk[100];   //[muo_n]
   Int_t           muo_truth[100];   //[muo_n]
   Int_t           muo_trign[100];   //[muo_n]
   Int_t           muo_trig[100][500];   //[muo_n]
   Int_t           muo_ID[100][24];   //[muo_n]
   Int_t           muo_ValidMuonHitsCm[100];   //[muo_n]
   Int_t           muo_ValidTrackerHitsCm[100];   //[muo_n]
   Int_t           muo_ValidPixelHitsCm[100];   //[muo_n]
   Int_t           muo_ChambersMatched[100];   //[muo_n]
   Double_t        muo_d0bsCm[100];   //[muo_n]
   Double_t        muo_d0OriginCm[100];   //[muo_n]
   Double_t        muo_dzbsCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersMeasCm[100];   //[muo_n]
   Int_t           muo_TrackerLayersNotMeasCm[100];   //[muo_n]
   Double_t        muo_Cocktail_pt[100];   //[muo_n]
   Double_t        muo_Cocktail_phi[100];   //[muo_n]
   Double_t        muo_Cocktail_eta[100];   //[muo_n]
   Double_t        muo_Valid_fraction[100];   //[muo_n]
   Int_t           PFmuo_n;
   Double_t        PFmuo_p[100];   //[PFmuo_n]
   Double_t        PFmuo_pt[100];   //[PFmuo_n]
   Double_t        PFmuo_E[100];   //[PFmuo_n]
   Double_t        PFmuo_Et[100];   //[PFmuo_n]
   Double_t        PFmuo_px[100];   //[PFmuo_n]
   Double_t        PFmuo_py[100];   //[PFmuo_n]
   Double_t        PFmuo_pz[100];   //[PFmuo_n]
   Double_t        PFmuo_eta[100];   //[PFmuo_n]
   Double_t        PFmuo_phi[100];   //[PFmuo_n]
   Double_t        PFmuo_Charge[100];   //[PFmuo_n]
   Double_t        PFmuo_particleIso[100];   //[PFmuo_n]
   Double_t        PFmuo_chadIso[100];   //[PFmuo_n]
   Double_t        PFmuo_nhadIso[100];   //[PFmuo_n]
   Double_t        PFmuo_gamIso[100];   //[PFmuo_n]
   Double_t        PFmuo_RelTrkIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkIso[100];   //[PFmuo_n]
   Double_t        PFmuo_ECalIso[100];   //[PFmuo_n]
   Double_t        PFmuo_HCalIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_ECalIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_HCalIsoDep[100];   //[PFmuo_n]
   Double_t        PFmuo_AllIso[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkChiNormCm[100];   //[PFmuo_n]
   Double_t        PFmuo_TrkChiNormTk[100];   //[PFmuo_n]
   Double_t        PFmuo_d0Cm[100];   //[PFmuo_n]
   Double_t        PFmuo_d0Tk[100];   //[PFmuo_n]
   Double_t        PFmuo_sd0Cm[100];   //[PFmuo_n]
   Double_t        PFmuo_sd0Tk[100];   //[PFmuo_n]
   Double_t        PFmuo_calocomp[100];   //[PFmuo_n]
   Double_t        PFmuo_calotower_e[100];   //[PFmuo_n]
   Int_t           PFmuo_hitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_hitsTk[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidMuonHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidTrackerHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ValidPixelHitsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_ChambersMatched[100];   //[PFmuo_n]
   Double_t        PFmuo_d0bsCm[100];   //[PFmuo_n]
   Double_t        PFmuo_d0OriginCm[100];   //[PFmuo_n]
   Double_t        PFmuo_dzbsCm[100];   //[PFmuo_n]
   Int_t           PFmuo_TrackerLayersMeasCm[100];   //[PFmuo_n]
   Int_t           PFmuo_TrackerLayersNotMeasCm[100];   //[PFmuo_n]
   Double_t        PFmuo_Valid_fraction[100];   //[PFmuo_n]
   Int_t           tau_n;
   Double_t        tau_p[100];   //[tau_n]
   Double_t        tau_pt[100];   //[tau_n]
   Double_t        tau_E[100];   //[tau_n]
   Double_t        tau_Et[100];   //[tau_n]
   Double_t        tau_Px[100];   //[tau_n]
   Double_t        tau_Py[100];   //[tau_n]
   Double_t        tau_Pz[100];   //[tau_n]
   Double_t        tau_Eta[100];   //[tau_n]
   Double_t        tau_Phi[100];   //[tau_n]
   Int_t           tau_DecayMode[100];   //[tau_n]
   Double_t        tau_vx[100];   //[tau_n]
   Double_t        tau_vy[100];   //[tau_n]
   Double_t        tau_vz[100];   //[tau_n]
   Double_t        tau_vx2[100];   //[tau_n]
   Double_t        tau_vy2[100];   //[tau_n]
   Double_t        tau_vz2[100];   //[tau_n]
   Double_t        tau_ECalIso[100];   //[tau_n]
   Double_t        tau_HCalIso[100];   //[tau_n]
   Double_t        tau_AllIso[100];   //[tau_n]
   Double_t        tau_TrackIso[100];   //[tau_n]
   Double_t        tau_ParticleIso[100];   //[tau_n]
   Double_t        tau_ChadIso[100];   //[tau_n]
   Double_t        tau_NhadIso[100];   //[tau_n]
   Double_t        tau_GamIso[100];   //[tau_n]
   Double_t        susyScanM0;
   Double_t        susyScanM12;
   Double_t        susyScanA0;
   Double_t        susyScanCrossSection;
   Double_t        susyScanMu;
   Double_t        susyScanRun;
   Double_t        susyScantanbeta;

   // List of branches
   TBranch        *b_global_weight;   //!
   TBranch        *b_global_procID;   //!
   TBranch        *b_global_pthat;   //!
   TBranch        *b_global_bfield;   //!
   TBranch        *b_global_store;   //!
   TBranch        *b_global_run;   //!
   TBranch        *b_global_event;   //!
   TBranch        *b_global_bx;   //!
   TBranch        *b_global_orbit;   //!
   TBranch        *b_global_exp;   //!
   TBranch        *b_global_isdata;   //!
   TBranch        *b_lumi_section;   //!
   TBranch        *b_lumi_del;   //!
   TBranch        *b_lumi_rec;   //!
   TBranch        *b_lumi_delerr;   //!
   TBranch        *b_lumi_recerr;   //!
   TBranch        *b_pu_bunchx;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_pu_num_int;   //!
   TBranch        *b_pu_inst_Lumi;   //!
   TBranch        *b_pu_zPos;   //!
   TBranch        *b_pu_sumPthi;   //!
   TBranch        *b_pu_sumPtlo;   //!
   TBranch        *b_noise_pLoose;   //!
   TBranch        *b_noise_pTight;   //!
   TBranch        *b_noise_pHigh;   //!
   TBranch        *b_noise_ecal_r9;   //!
   TBranch        *b_noise_ecal_E;   //!
   TBranch        *b_noise_ecal_pt;   //!
   TBranch        *b_noise_ecal_px;   //!
   TBranch        *b_noise_ecal_py;   //!
   TBranch        *b_noise_ecal_pz;   //!
   TBranch        *b_noise_ecal_eta;   //!
   TBranch        *b_noise_ecal_phi;   //!
   TBranch        *b_noise_ecal_time;   //!
   TBranch        *b_noise_ecal_chi;   //!
   TBranch        *b_noise_ecal_flag;   //!
   TBranch        *b_noise_ecal_ieta;   //!
   TBranch        *b_noise_ecal_iphi;   //!
   TBranch        *b_noise_hcal_eventChargeFraction;   //!
   TBranch        *b_noise_hcal_eventEMEnergy;   //!
   TBranch        *b_noise_hcal_eventEMFraction;   //!
   TBranch        *b_noise_hcal_eventHadEnergy;   //!
   TBranch        *b_noise_hcal_eventTrackEnergy;   //!
   TBranch        *b_noise_hcal_flatNoiseSumE;   //!
   TBranch        *b_noise_hcal_flatNoiseSumEt;   //!
   TBranch        *b_noise_hcal_HasBadRBXTS4TS5;   //!
   TBranch        *b_noise_hcal_isolatedNoiseSumE;   //!
   TBranch        *b_noise_hcal_isolatedNoiseSumEt;   //!
   TBranch        *b_noise_hcal_max10GeVHitTime;   //!
   TBranch        *b_noise_hcal_max25GeVHitTime;   //!
   TBranch        *b_noise_hcal_maxE10TS;   //!
   TBranch        *b_noise_hcal_maxE2Over10TS;   //!
   TBranch        *b_noise_hcal_maxE2TS;   //!
   TBranch        *b_noise_hcal_maxHPDHits;   //!
   TBranch        *b_noise_hcal_maxHPDNoOtherHits;   //!
   TBranch        *b_noise_hcal_maxRBXHits;   //!
   TBranch        *b_noise_hcal_maxZeros;   //!
   TBranch        *b_noise_hcal_min10GeVHitTime;   //!
   TBranch        *b_noise_hcal_min25GeVHitTime;   //!
   TBranch        *b_noise_hcal_minE10TS;   //!
   TBranch        *b_noise_hcal_minE2Over10TS;   //!
   TBranch        *b_noise_hcal_minE2TS;   //!
   TBranch        *b_noise_hcal_minHPDEMF;   //!
   TBranch        *b_noise_hcal_minRBXEMF;   //!
   TBranch        *b_noise_hcal_noiseFilterStatus;   //!
   TBranch        *b_noise_hcal_noiseType;   //!
   TBranch        *b_noise_hcal_num10GeVHits;   //!
   TBranch        *b_noise_hcal_num25GeVHits;   //!
   TBranch        *b_noise_hcal_numFlatNoiseChannels;   //!
   TBranch        *b_noise_hcal_numIsolatedNoiseChannels;   //!
   TBranch        *b_noise_hcal_numProblematicRBXs;   //!
   TBranch        *b_noise_hcal_numSpikeNoiseChannels;   //!
   TBranch        *b_noise_hcal_numTriangleNoiseChannels;   //!
   TBranch        *b_noise_hcal_numTS4TS5NoiseChannels;   //!
   TBranch        *b_noise_hcal_passHighLevelNoiseFilter;   //!
   TBranch        *b_noise_hcal_passLooseNoiseFilter;   //!
   TBranch        *b_noise_hcal_passTightNoiseFilter;   //!
   TBranch        *b_noise_hcal_rms10GeVHitTime;   //!
   TBranch        *b_noise_hcal_rms25GeVHitTime;   //!
   TBranch        *b_noise_hcal_spikeNoiseSumE;   //!
   TBranch        *b_noise_hcal_spikeNoiseSumEt;   //!
   TBranch        *b_noise_hcal_triangleNoiseSumE;   //!
   TBranch        *b_noise_hcal_triangleNoiseSumEt;   //!
   TBranch        *b_noise_hcal_TS4TS5NoiseSumE;   //!
   TBranch        *b_noise_hcal_TS4TS5NoiseSumEt;   //!
   TBranch        *b_trig_HLTName;   //!
   TBranch        *b_trig_n;   //!
   TBranch        *b_trig_L1prescale;   //!
   TBranch        *b_trig_HLTprescale;   //!
   TBranch        *b_trig_name;   //!
   TBranch        *b_trig_filter;   //!
   TBranch        *b_trig_pt;   //!
   TBranch        *b_trig_eta;   //!
   TBranch        *b_trig_phi;   //!
   TBranch        *b_truth_n;   //!
   TBranch        *b_truth_pdgid;   //!
   TBranch        *b_truth_bvtxid;   //!
   TBranch        *b_truth_evtxid;   //!
   TBranch        *b_truth_E;   //!
   TBranch        *b_truth_Et;   //!
   TBranch        *b_truth_p;   //!
   TBranch        *b_truth_pt;   //!
   TBranch        *b_truth_px;   //!
   TBranch        *b_truth_py;   //!
   TBranch        *b_truth_pz;   //!
   TBranch        *b_truth_eta;   //!
   TBranch        *b_truth_phi;   //!
   TBranch        *b_truth_m;   //!
   TBranch        *b_truthl_n;   //!
   TBranch        *b_truthl_ori;   //!
   TBranch        *b_truthl_pdgid;   //!
   TBranch        *b_truthl_E;   //!
   TBranch        *b_truthl_Et;   //!
   TBranch        *b_truthl_p;   //!
   TBranch        *b_truthl_pt;   //!
   TBranch        *b_truthl_px;   //!
   TBranch        *b_truthl_py;   //!
   TBranch        *b_truthl_pz;   //!
   TBranch        *b_truthl_eta;   //!
   TBranch        *b_truthl_phi;   //!
   TBranch        *b_pdf_id1;   //!
   TBranch        *b_pdf_id2;   //!
   TBranch        *b_pdf_x1;   //!
   TBranch        *b_pdf_x2;   //!
   TBranch        *b_pdf_f1;   //!
   TBranch        *b_pdf_f2;   //!
   TBranch        *b_pdf_scale;   //!
   TBranch        *b_vtx_n;   //!
   TBranch        *b_vtx_ntr;   //!
   TBranch        *b_vtx_fake;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_chi;   //!
   TBranch        *b_bs_x;   //!
   TBranch        *b_bs_y;   //!
   TBranch        *b_bs_z;   //!
   TBranch        *b_tracks_n;   //!
   TBranch        *b_tracks_hqf;   //!
   TBranch        *b_met_et;   //!
   TBranch        *b_met_ex;   //!
   TBranch        *b_met_ey;   //!
   TBranch        *b_met_eta;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_sumetsig;   //!
   TBranch        *b_met_etsignif;   //!
   TBranch        *b_met_CaloMETInmHF;   //!
   TBranch        *b_met_CaloMETInpHF;   //!
   TBranch        *b_met_CaloMETPhiInmHF;   //!
   TBranch        *b_met_CaloMETPhiInpHF;   //!
   TBranch        *b_met_CaloSETInmHF;   //!
   TBranch        *b_met_CaloSETInpHF;   //!
   TBranch        *b_met_emEtFraction;   //!
   TBranch        *b_met_etFractionHadronic;   //!
   TBranch        *b_met_maxEtInEmTowers;   //!
   TBranch        *b_met_maxEtInHadTowers;   //!
   TBranch        *b_met_emEtInHF;   //!
   TBranch        *b_met_emEtInEE;   //!
   TBranch        *b_met_emEtInEB;   //!
   TBranch        *b_met_hadEtInHF;   //!
   TBranch        *b_met_hadEtInHE;   //!
   TBranch        *b_met_hadEtInHO;   //!
   TBranch        *b_met_hadEtInHB;   //!
   TBranch        *b_met_ChargedEMEtFraction;   //!
   TBranch        *b_met_ChargedHadEtFraction;   //!
   TBranch        *b_met_MuonEtFraction;   //!
   TBranch        *b_met_NeutralEMFraction;   //!
   TBranch        *b_met_NeutralHadEtFraction;   //!
   TBranch        *b_met_Type6EtFraction;   //!
   TBranch        *b_met_Type7EtFraction;   //!
   TBranch        *b_calojet_n;   //!
   TBranch        *b_calojet_E;   //!
   TBranch        *b_calojet_Et;   //!
   TBranch        *b_calojet_p;   //!
   TBranch        *b_calojet_pt;   //!
   TBranch        *b_calojet_pt_raw;   //!
   TBranch        *b_calojet_px;   //!
   TBranch        *b_calojet_py;   //!
   TBranch        *b_calojet_pz;   //!
   TBranch        *b_calojet_eta;   //!
   TBranch        *b_calojet_phi;   //!
   TBranch        *b_calojet_fem;   //!
   TBranch        *b_calojet_fhad;   //!
   TBranch        *b_calojet_btag;   //!
   TBranch        *b_calojet_charge;   //!
   TBranch        *b_calojet_fHPD;   //!
   TBranch        *b_calojet_fRBX;   //!
   TBranch        *b_calojet_n90hits;   //!
   TBranch        *b_calojet_n90;   //!
   TBranch        *b_calojet_flav;   //!
   TBranch        *b_calojet_truth;   //!
   TBranch        *b_calojet_const;   //!
   TBranch        *b_calojet_ID;   //!
   TBranch        *b_pfjet_n;   //!
   TBranch        *b_pfjet_E;   //!
   TBranch        *b_pfjet_Et;   //!
   TBranch        *b_pfjet_p;   //!
   TBranch        *b_pfjet_pt;   //!
   TBranch        *b_pfjet_pt_raw;   //!
   TBranch        *b_pfjet_px;   //!
   TBranch        *b_pfjet_py;   //!
   TBranch        *b_pfjet_pz;   //!
   TBranch        *b_pfjet_eta;   //!
   TBranch        *b_pfjet_phi;   //!
   TBranch        *b_pfjet_btag;   //!
   TBranch        *b_pfjet_charge;   //!
   TBranch        *b_pfjet_n90;   //!
   TBranch        *b_pfjet_flav;   //!
   TBranch        *b_pfjet_truth;   //!
   TBranch        *b_pfjet_const;   //!
   TBranch        *b_pfjet_PFN;   //!
   TBranch        *b_pfjet_PFF;   //!
   TBranch        *b_truthjet_n;   //!
   TBranch        *b_truthjet_E;   //!
   TBranch        *b_truthjet_Et;   //!
   TBranch        *b_truthjet_p;   //!
   TBranch        *b_truthjet_pt;   //!
   TBranch        *b_truthjet_px;   //!
   TBranch        *b_truthjet_py;   //!
   TBranch        *b_truthjet_pz;   //!
   TBranch        *b_truthjet_eta;   //!
   TBranch        *b_truthjet_phi;   //!
   TBranch        *b_fatjet_n;   //!
   TBranch        *b_fatjet_nsub;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_px;   //!
   TBranch        *b_fatjet_py;   //!
   TBranch        *b_fatjet_pz;   //!
   TBranch        *b_fatjet_E;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_sub_pt;   //!
   TBranch        *b_fatjet_sub_px;   //!
   TBranch        *b_fatjet_sub_py;   //!
   TBranch        *b_fatjet_sub_pz;   //!
   TBranch        *b_fatjet_sub_E;   //!
   TBranch        *b_fatjet_sub_eta;   //!
   TBranch        *b_fatjet_sub_phi;   //!
   TBranch        *b_fatjet_sub_fem;   //!
   TBranch        *b_fatjet_sub_fhad;   //!
   TBranch        *b_fatjet_sub_btag;   //!
   TBranch        *b_fatjet_sub_n90;   //!
   TBranch        *b_fatjet_sub_fHPD;   //!
   TBranch        *b_fatjet_sub_fRBX;   //!
   TBranch        *b_SC_n;   //!
   TBranch        *b_SC_truth;   //!
   TBranch        *b_SC_E;   //!
   TBranch        *b_SC_phi;   //!
   TBranch        *b_SC_eta;   //!
   TBranch        *b_ele_n;   //!
   TBranch        *b_ele_E;   //!
   TBranch        *b_ele_Et;   //!
   TBranch        *b_ele_p;   //!
   TBranch        *b_ele_pt;   //!
   TBranch        *b_ele_px;   //!
   TBranch        *b_ele_py;   //!
   TBranch        *b_ele_pz;   //!
   TBranch        *b_ele_eta;   //!
   TBranch        *b_ele_phi;   //!
   TBranch        *b_ele_charge;   //!
   TBranch        *b_ele_TrkChiNorm;   //!
   TBranch        *b_ele_d0vtx;   //!
   TBranch        *b_ele_d0bs;   //!
   TBranch        *b_ele_sd0;   //!
   TBranch        *b_ele_hits;   //!
   TBranch        *b_ele_truth;   //!
   TBranch        *b_ele_isECal;   //!
   TBranch        *b_ele_isTracker;   //!
   TBranch        *b_ele_ID;   //!
   TBranch        *b_ele_ValidHitFirstPxlB;   //!
   TBranch        *b_ele_TrkExpHitsInner;   //!
   TBranch        *b_ele_HCalOverEm;   //!
   TBranch        *b_ele_Dr03TkSumPt;   //!
   TBranch        *b_ele_Dr04HCalSumEt;   //!
   TBranch        *b_ele_Dr03HCalSumEt;   //!
   TBranch        *b_ele_Dr04ECalSumEt;   //!
   TBranch        *b_ele_Dr03ECalSumEt;   //!
   TBranch        *b_ele_SigmaIetaIeta;   //!
   TBranch        *b_ele_dEtaSCTrackAtVtx;   //!
   TBranch        *b_ele_dPhiSCTrackAtVtx;   //!
   TBranch        *b_ele_dr03HcalDepth1;   //!
   TBranch        *b_ele_dr03HcalDepth2;   //!
   TBranch        *b_ele_e2x5Max;   //!
   TBranch        *b_ele_e5x5;   //!
   TBranch        *b_ele_e1x5;   //!
   TBranch        *b_ele_caloEt;   //!
   TBranch        *b_ele_SCeta;   //!
   TBranch        *b_ele_convdist;   //!
   TBranch        *b_ele_convdcot;   //!
   TBranch        *b_ele_convr;   //!
   TBranch        *b_ele_fbrem;   //!
   TBranch        *b_ele_trign;   //!
   TBranch        *b_ele_trig;   //!
   TBranch        *b_ele_SC;   //!
   TBranch        *b_ele_numberOfHits;   //!
   TBranch        *b_pfele_n;   //!
   TBranch        *b_pfele_p;   //!
   TBranch        *b_pfele_E;   //!
   TBranch        *b_pfele_Et;   //!
   TBranch        *b_pfele_pt;   //!
   TBranch        *b_pfele_px;   //!
   TBranch        *b_pfele_py;   //!
   TBranch        *b_pfele_pz;   //!
   TBranch        *b_pfele_eta;   //!
   TBranch        *b_pfele_phi;   //!
   TBranch        *b_pfele_charge;   //!
   TBranch        *b_pfele_truth;   //!
   TBranch        *b_pfele_trign;   //!
   TBranch        *b_pfele_trig;   //!
   TBranch        *b_pfele_SC;   //!
   TBranch        *b_muo_n;   //!
   TBranch        *b_muo_E;   //!
   TBranch        *b_muo_Et;   //!
   TBranch        *b_muo_p;   //!
   TBranch        *b_muo_pt;   //!
   TBranch        *b_muo_px;   //!
   TBranch        *b_muo_py;   //!
   TBranch        *b_muo_pz;   //!
   TBranch        *b_muo_eta;   //!
   TBranch        *b_muo_phi;   //!
   TBranch        *b_muo_charge;   //!
   TBranch        *b_muo_RelTrkIso;   //!
   TBranch        *b_muo_TrkIso;   //!
   TBranch        *b_muo_ECalIso;   //!
   TBranch        *b_muo_HCalIso;   //!
   TBranch        *b_muo_TrkIsoDep;   //!
   TBranch        *b_muo_ECalIsoDep;   //!
   TBranch        *b_muo_HCalIsoDep;   //!
   TBranch        *b_muo_AllIso;   //!
   TBranch        *b_muo_TrkChiNormCm;   //!
   TBranch        *b_muo_TrkChiNormTk;   //!
   TBranch        *b_muo_d0Cm;   //!
   TBranch        *b_muo_d0Tk;   //!
   TBranch        *b_muo_sd0Cm;   //!
   TBranch        *b_muo_sd0Tk;   //!
   TBranch        *b_muo_calocomp;   //!
   TBranch        *b_muo_calotower_e;   //!
   TBranch        *b_muo_prompttight;   //!
   TBranch        *b_muo_hitsCm;   //!
   TBranch        *b_muo_hitsTk;   //!
   TBranch        *b_muo_truth;   //!
   TBranch        *b_muo_trign;   //!
   TBranch        *b_muo_trig;   //!
   TBranch        *b_muo_ID;   //!
   TBranch        *b_muo_ValidMuonHitsCm;   //!
   TBranch        *b_muo_ValidTrackerHitsCm;   //!
   TBranch        *b_muo_ValidPixelHitsCm;   //!
   TBranch        *b_muo_ChambersMatched;   //!
   TBranch        *b_muo_d0bsCm;   //!
   TBranch        *b_muo_d0OriginCm;   //!
   TBranch        *b_muo_dzbsCm;   //!
   TBranch        *b_muo_TrackerLayersMeasCm;   //!
   TBranch        *b_muo_TrackerLayersNotMeasCm;   //!
   TBranch        *b_muo_Cocktail_pt;   //!
   TBranch        *b_muo_Cocktail_phi;   //!
   TBranch        *b_muo_Cocktail_eta;   //!
   TBranch        *b_muo_Valid_fraction;   //!
   TBranch        *b_PFmuo_n;   //!
   TBranch        *b_PFmuo_p;   //!
   TBranch        *b_PFmuo_pt;   //!
   TBranch        *b_PFmuo_E;   //!
   TBranch        *b_PFmuo_Et;   //!
   TBranch        *b_PFmuo_px;   //!
   TBranch        *b_PFmuo_py;   //!
   TBranch        *b_PFmuo_pz;   //!
   TBranch        *b_PFmuo_eta;   //!
   TBranch        *b_PFmuo_phi;   //!
   TBranch        *b_PFmuo_Charge;   //!
   TBranch        *b_PFmuo_particleIso;   //!
   TBranch        *b_PFmuo_chadIso;   //!
   TBranch        *b_PFmuo_nhadIso;   //!
   TBranch        *b_PFmuo_gamIso;   //!
   TBranch        *b_PFmuo_RelTrkIso;   //!
   TBranch        *b_PFmuo_TrkIso;   //!
   TBranch        *b_PFmuo_ECalIso;   //!
   TBranch        *b_PFmuo_HCalIso;   //!
   TBranch        *b_PFmuo_TrkIsoDep;   //!
   TBranch        *b_PFmuo_ECalIsoDep;   //!
   TBranch        *b_PFmuo_HCalIsoDep;   //!
   TBranch        *b_PFmuo_AllIso;   //!
   TBranch        *b_PFmuo_TrkChiNormCm;   //!
   TBranch        *b_PFmuo_TrkChiNormTk;   //!
   TBranch        *b_PFmuo_d0Cm;   //!
   TBranch        *b_PFmuo_d0Tk;   //!
   TBranch        *b_PFmuo_sd0Cm;   //!
   TBranch        *b_PFmuo_sd0Tk;   //!
   TBranch        *b_PFmuo_calocomp;   //!
   TBranch        *b_PFmuo_calotower_e;   //!
   TBranch        *b_PFmuo_hitsCm;   //!
   TBranch        *b_PFmuo_hitsTk;   //!
   TBranch        *b_PFmuo_ValidMuonHitsCm;   //!
   TBranch        *b_PFmuo_ValidTrackerHitsCm;   //!
   TBranch        *b_PFmuo_ValidPixelHitsCm;   //!
   TBranch        *b_PFmuo_ChambersMatched;   //!
   TBranch        *b_PFmuo_d0bsCm;   //!
   TBranch        *b_PFmuo_d0OriginCm;   //!
   TBranch        *b_PFmuo_dzbsCm;   //!
   TBranch        *b_PFmuo_TrackerLayersMeasCm;   //!
   TBranch        *b_PFmuo_TrackerLayersNotMeasCm;   //!
   TBranch        *b_PFmuo_Valid_fraction;   //!
   TBranch        *b_tau_n;   //!
   TBranch        *b_tau_p;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_E;   //!
   TBranch        *b_tau_Et;   //!
   TBranch        *b_tau_Px;   //!
   TBranch        *b_tau_Py;   //!
   TBranch        *b_tau_Pz;   //!
   TBranch        *b_tau_Eta;   //!
   TBranch        *b_tau_Phi;   //!
   TBranch        *b_tau_DecayMode;   //!
   TBranch        *b_tau_vx;   //!
   TBranch        *b_tau_vy;   //!
   TBranch        *b_tau_vz;   //!
   TBranch        *b_tau_vx2;   //!
   TBranch        *b_tau_vy2;   //!
   TBranch        *b_tau_vz2;   //!
   TBranch        *b_tau_ECalIso;   //!
   TBranch        *b_tau_HCalIso;   //!
   TBranch        *b_tau_AllIso;   //!
   TBranch        *b_tau_TrackIso;   //!
   TBranch        *b_tau_ParticleIso;   //!
   TBranch        *b_tau_ChadIso;   //!
   TBranch        *b_tau_NhadIso;   //!
   TBranch        *b_tau_GamIso;   //!
   TBranch        *b_susyScanM0;   //!
   TBranch        *b_susyScanM12;   //!
   TBranch        *b_susyScanA0;   //!
   TBranch        *b_susyScanCrossSection;   //!
   TBranch        *b_susyScanMu;   //!
   TBranch        *b_susyScanRun;   //!
   TBranch        *b_susyScantanbeta;   //!

   allData(TTree *tree=0);
   virtual ~allData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TreeContent_cxx
TreeContent::TreeContent(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out.root");
      if (!f) {
         f = new TFile("out.root");
         f->cd("out.root:/ACSkimAnalysis");
      }
      tree = (TTree*)gDirectory->Get("allData");

   }
   Init(tree);
}

TreeContent::~TreeContent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeContent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeContent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeContent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("global_weight", &global_weight, &b_global_weight);
   fChain->SetBranchAddress("global_procID", &global_procID, &b_global_procID);
   fChain->SetBranchAddress("global_pthat", &global_pthat, &b_global_pthat);
   fChain->SetBranchAddress("global_bfield", &global_bfield, &b_global_bfield);
   fChain->SetBranchAddress("global_store", &global_store, &b_global_store);
   fChain->SetBranchAddress("global_run", &global_run, &b_global_run);
   fChain->SetBranchAddress("global_event", &global_event, &b_global_event);
   fChain->SetBranchAddress("global_bx", &global_bx, &b_global_bx);
   fChain->SetBranchAddress("global_orbit", &global_orbit, &b_global_orbit);
   fChain->SetBranchAddress("global_exp", &global_exp, &b_global_exp);
   fChain->SetBranchAddress("global_isdata", &global_isdata, &b_global_isdata);
   fChain->SetBranchAddress("lumi_section", &lumi_section, &b_lumi_section);
   fChain->SetBranchAddress("lumi_del", &lumi_del, &b_lumi_del);
   fChain->SetBranchAddress("lumi_rec", &lumi_rec, &b_lumi_rec);
   fChain->SetBranchAddress("lumi_delerr", &lumi_delerr, &b_lumi_delerr);
   fChain->SetBranchAddress("lumi_recerr", &lumi_recerr, &b_lumi_recerr);
   fChain->SetBranchAddress("noise_pLoose", &noise_pLoose, &b_noise_pLoose);
   fChain->SetBranchAddress("noise_pTight", &noise_pTight, &b_noise_pTight);
   fChain->SetBranchAddress("noise_pHigh", &noise_pHigh, &b_noise_pHigh);
   fChain->SetBranchAddress("noise_ecal_r9", &noise_ecal_r9, &b_noise_ecal_r9);
   fChain->SetBranchAddress("noise_ecal_E", &noise_ecal_E, &b_noise_ecal_E);
   fChain->SetBranchAddress("noise_ecal_pt", &noise_ecal_pt, &b_noise_ecal_pt);
   fChain->SetBranchAddress("noise_ecal_px", &noise_ecal_px, &b_noise_ecal_px);
   fChain->SetBranchAddress("noise_ecal_py", &noise_ecal_py, &b_noise_ecal_py);
   fChain->SetBranchAddress("noise_ecal_pz", &noise_ecal_pz, &b_noise_ecal_pz);
   fChain->SetBranchAddress("noise_ecal_eta", &noise_ecal_eta, &b_noise_ecal_eta);
   fChain->SetBranchAddress("noise_ecal_phi", &noise_ecal_phi, &b_noise_ecal_phi);
   fChain->SetBranchAddress("noise_ecal_time", &noise_ecal_time, &b_noise_ecal_time);
   fChain->SetBranchAddress("noise_ecal_chi", &noise_ecal_chi, &b_noise_ecal_chi);
   fChain->SetBranchAddress("noise_ecal_flag", &noise_ecal_flag, &b_noise_ecal_flag);
   fChain->SetBranchAddress("noise_ecal_ieta", &noise_ecal_ieta, &b_noise_ecal_ieta);
   fChain->SetBranchAddress("noise_ecal_iphi", &noise_ecal_iphi, &b_noise_ecal_iphi);
   fChain->SetBranchAddress("trig_HLTName", trig_HLTName, &b_trig_HLTName);
   fChain->SetBranchAddress("trig_n", &trig_n, &b_trig_n);
   fChain->SetBranchAddress("trig_L1prescale", trig_L1prescale, &b_trig_L1prescale);
   fChain->SetBranchAddress("trig_HLTprescale", trig_HLTprescale, &b_trig_HLTprescale);
   fChain->SetBranchAddress("trig_name", trig_name, &b_trig_name);
   fChain->SetBranchAddress("trig_filter", trig_filter, &b_trig_filter);
   fChain->SetBranchAddress("trig_pt", trig_pt, &b_trig_pt);
   fChain->SetBranchAddress("trig_eta", trig_eta, &b_trig_eta);
   fChain->SetBranchAddress("trig_phi", trig_phi, &b_trig_phi);
   fChain->SetBranchAddress("truth_n", &truth_n, &b_truth_n);
   fChain->SetBranchAddress("truth_pdgid", truth_pdgid, &b_truth_pdgid);
   fChain->SetBranchAddress("truth_bvtxid", truth_bvtxid, &b_truth_bvtxid);
   fChain->SetBranchAddress("truth_evtxid", truth_evtxid, &b_truth_evtxid);
   fChain->SetBranchAddress("truth_E", truth_E, &b_truth_E);
   fChain->SetBranchAddress("truth_Et", truth_Et, &b_truth_Et);
   fChain->SetBranchAddress("truth_p", truth_p, &b_truth_p);
   fChain->SetBranchAddress("truth_pt", truth_pt, &b_truth_pt);
   fChain->SetBranchAddress("truth_px", truth_px, &b_truth_px);
   fChain->SetBranchAddress("truth_py", truth_py, &b_truth_py);
   fChain->SetBranchAddress("truth_pz", truth_pz, &b_truth_pz);
   fChain->SetBranchAddress("truth_eta", truth_eta, &b_truth_eta);
   fChain->SetBranchAddress("truth_phi", truth_phi, &b_truth_phi);
   fChain->SetBranchAddress("truth_m", truth_m, &b_truth_m);
   fChain->SetBranchAddress("truthl_n", &truthl_n, &b_truthl_n);
   fChain->SetBranchAddress("truthl_ori", truthl_ori, &b_truthl_ori);
   fChain->SetBranchAddress("truthl_pdgid", truthl_pdgid, &b_truthl_pdgid);
   fChain->SetBranchAddress("truthl_E", truthl_E, &b_truthl_E);
   fChain->SetBranchAddress("truthl_Et", truthl_Et, &b_truthl_Et);
   fChain->SetBranchAddress("truthl_p", truthl_p, &b_truthl_p);
   fChain->SetBranchAddress("truthl_pt", truthl_pt, &b_truthl_pt);
   fChain->SetBranchAddress("truthl_px", truthl_px, &b_truthl_px);
   fChain->SetBranchAddress("truthl_py", truthl_py, &b_truthl_py);
   fChain->SetBranchAddress("truthl_pz", truthl_pz, &b_truthl_pz);
   fChain->SetBranchAddress("truthl_eta", truthl_eta, &b_truthl_eta);
   fChain->SetBranchAddress("truthl_phi", truthl_phi, &b_truthl_phi);
   fChain->SetBranchAddress("pdf_id1", &pdf_id1, &b_pdf_id1);
   fChain->SetBranchAddress("pdf_id2", &pdf_id2, &b_pdf_id2);
   fChain->SetBranchAddress("pdf_x1", &pdf_x1, &b_pdf_x1);
   fChain->SetBranchAddress("pdf_x2", &pdf_x2, &b_pdf_x2);
   fChain->SetBranchAddress("pdf_f1", &pdf_f1, &b_pdf_f1);
   fChain->SetBranchAddress("pdf_f2", &pdf_f2, &b_pdf_f2);
   fChain->SetBranchAddress("pdf_scale", &pdf_scale, &b_pdf_scale);
   fChain->SetBranchAddress("vtx_n", &vtx_n, &b_vtx_n);
   fChain->SetBranchAddress("vtx_ntr", vtx_ntr, &b_vtx_ntr);
   fChain->SetBranchAddress("vtx_fake", vtx_fake, &b_vtx_fake);
   fChain->SetBranchAddress("vtx_ndof", vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_chi", vtx_chi, &b_vtx_chi);
   fChain->SetBranchAddress("bs_x", &bs_x, &b_bs_x);
   fChain->SetBranchAddress("bs_y", &bs_y, &b_bs_y);
   fChain->SetBranchAddress("bs_z", &bs_z, &b_bs_z);
   fChain->SetBranchAddress("tracks_n", &tracks_n, &b_tracks_n);
   fChain->SetBranchAddress("tracks_hqf", &tracks_hqf, &b_tracks_hqf);
   fChain->SetBranchAddress("met_et", met_et, &b_met_et);
   fChain->SetBranchAddress("met_ex", met_ex, &b_met_ex);
   fChain->SetBranchAddress("met_ey", met_ey, &b_met_ey);
   fChain->SetBranchAddress("met_eta", met_eta, &b_met_eta);
   fChain->SetBranchAddress("met_phi", met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_sumet", met_sumet, &b_met_sumet);
   fChain->SetBranchAddress("met_sumetsig", met_sumetsig, &b_met_sumetsig);
   fChain->SetBranchAddress("met_etsignif", met_etsignif, &b_met_etsignif);
   fChain->SetBranchAddress("calojet_n", &calojet_n, &b_calojet_n);
   fChain->SetBranchAddress("calojet_E", calojet_E, &b_calojet_E);
   fChain->SetBranchAddress("calojet_Et", calojet_Et, &b_calojet_Et);
   fChain->SetBranchAddress("calojet_p", calojet_p, &b_calojet_p);
   fChain->SetBranchAddress("calojet_pt", calojet_pt, &b_calojet_pt);
   fChain->SetBranchAddress("calojet_pt_raw", calojet_pt_raw, &b_calojet_pt_raw);
   fChain->SetBranchAddress("calojet_px", calojet_px, &b_calojet_px);
   fChain->SetBranchAddress("calojet_py", calojet_py, &b_calojet_py);
   fChain->SetBranchAddress("calojet_pz", calojet_pz, &b_calojet_pz);
   fChain->SetBranchAddress("calojet_eta", calojet_eta, &b_calojet_eta);
   fChain->SetBranchAddress("calojet_phi", calojet_phi, &b_calojet_phi);
   fChain->SetBranchAddress("calojet_fem", calojet_fem, &b_calojet_fem);
   fChain->SetBranchAddress("calojet_fhad", calojet_fhad, &b_calojet_fhad);
   fChain->SetBranchAddress("calojet_btag", calojet_btag, &b_calojet_btag);
   fChain->SetBranchAddress("calojet_charge", calojet_charge, &b_calojet_charge);
   fChain->SetBranchAddress("calojet_fHPD", calojet_fHPD, &b_calojet_fHPD);
   fChain->SetBranchAddress("calojet_fRBX", calojet_fRBX, &b_calojet_fRBX);
   fChain->SetBranchAddress("calojet_n90hits", calojet_n90hits, &b_calojet_n90hits);
   fChain->SetBranchAddress("calojet_n90", calojet_n90, &b_calojet_n90);
   fChain->SetBranchAddress("calojet_flav", calojet_flav, &b_calojet_flav);
   fChain->SetBranchAddress("calojet_truth", calojet_truth, &b_calojet_truth);
   fChain->SetBranchAddress("calojet_const", calojet_const, &b_calojet_const);
   fChain->SetBranchAddress("calojet_ID", calojet_ID, &b_calojet_ID);
   fChain->SetBranchAddress("pfjet_n", &pfjet_n, &b_pfjet_n);
   fChain->SetBranchAddress("pfjet_E", pfjet_E, &b_pfjet_E);
   fChain->SetBranchAddress("pfjet_Et", pfjet_Et, &b_pfjet_Et);
   fChain->SetBranchAddress("pfjet_p", pfjet_p, &b_pfjet_p);
   fChain->SetBranchAddress("pfjet_pt", pfjet_pt, &b_pfjet_pt);
   fChain->SetBranchAddress("pfjet_pt_raw", pfjet_pt_raw, &b_pfjet_pt_raw);
   fChain->SetBranchAddress("pfjet_px", pfjet_px, &b_pfjet_px);
   fChain->SetBranchAddress("pfjet_py", pfjet_py, &b_pfjet_py);
   fChain->SetBranchAddress("pfjet_pz", pfjet_pz, &b_pfjet_pz);
   fChain->SetBranchAddress("pfjet_eta", pfjet_eta, &b_pfjet_eta);
   fChain->SetBranchAddress("pfjet_phi", pfjet_phi, &b_pfjet_phi);
   fChain->SetBranchAddress("pfjet_btag", pfjet_btag, &b_pfjet_btag);
   fChain->SetBranchAddress("pfjet_charge", pfjet_charge, &b_pfjet_charge);
   fChain->SetBranchAddress("pfjet_n90", pfjet_n90, &b_pfjet_n90);
   fChain->SetBranchAddress("pfjet_flav", pfjet_flav, &b_pfjet_flav);
   fChain->SetBranchAddress("pfjet_truth", pfjet_truth, &b_pfjet_truth);
   fChain->SetBranchAddress("pfjet_const", pfjet_const, &b_pfjet_const);
   fChain->SetBranchAddress("pfjet_PFN", pfjet_PFN, &b_pfjet_PFN);
   fChain->SetBranchAddress("pfjet_PFF", pfjet_PFF, &b_pfjet_PFF);
   fChain->SetBranchAddress("truthjet_n", &truthjet_n, &b_truthjet_n);
   fChain->SetBranchAddress("truthjet_E", truthjet_E, &b_truthjet_E);
   fChain->SetBranchAddress("truthjet_Et", truthjet_Et, &b_truthjet_Et);
   fChain->SetBranchAddress("truthjet_p", truthjet_p, &b_truthjet_p);
   fChain->SetBranchAddress("truthjet_pt", truthjet_pt, &b_truthjet_pt);
   fChain->SetBranchAddress("truthjet_px", truthjet_px, &b_truthjet_px);
   fChain->SetBranchAddress("truthjet_py", truthjet_py, &b_truthjet_py);
   fChain->SetBranchAddress("truthjet_pz", truthjet_pz, &b_truthjet_pz);
   fChain->SetBranchAddress("truthjet_eta", truthjet_eta, &b_truthjet_eta);
   fChain->SetBranchAddress("truthjet_phi", truthjet_phi, &b_truthjet_phi);
   fChain->SetBranchAddress("fatjet_n", &fatjet_n, &b_fatjet_n);
   fChain->SetBranchAddress("fatjet_nsub", fatjet_nsub, &b_fatjet_nsub);
   fChain->SetBranchAddress("fatjet_pt", fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_px", fatjet_px, &b_fatjet_px);
   fChain->SetBranchAddress("fatjet_py", fatjet_py, &b_fatjet_py);
   fChain->SetBranchAddress("fatjet_pz", fatjet_pz, &b_fatjet_pz);
   fChain->SetBranchAddress("fatjet_E", fatjet_E, &b_fatjet_E);
   fChain->SetBranchAddress("fatjet_eta", fatjet_eta, &b_fatjet_eta);
   fChain->SetBranchAddress("fatjet_phi", fatjet_phi, &b_fatjet_phi);
   fChain->SetBranchAddress("fatjet_sub_pt", fatjet_sub_pt, &b_fatjet_sub_pt);
   fChain->SetBranchAddress("fatjet_sub_px", fatjet_sub_px, &b_fatjet_sub_px);
   fChain->SetBranchAddress("fatjet_sub_py", fatjet_sub_py, &b_fatjet_sub_py);
   fChain->SetBranchAddress("fatjet_sub_pz", fatjet_sub_pz, &b_fatjet_sub_pz);
   fChain->SetBranchAddress("fatjet_sub_E", fatjet_sub_E, &b_fatjet_sub_E);
   fChain->SetBranchAddress("fatjet_sub_eta", fatjet_sub_eta, &b_fatjet_sub_eta);
   fChain->SetBranchAddress("fatjet_sub_phi", fatjet_sub_phi, &b_fatjet_sub_phi);
   fChain->SetBranchAddress("fatjet_sub_fem", fatjet_sub_fem, &b_fatjet_sub_fem);
   fChain->SetBranchAddress("fatjet_sub_fhad", fatjet_sub_fhad, &b_fatjet_sub_fhad);
   fChain->SetBranchAddress("fatjet_sub_btag", fatjet_sub_btag, &b_fatjet_sub_btag);
   fChain->SetBranchAddress("fatjet_sub_n90", fatjet_sub_n90, &b_fatjet_sub_n90);
   fChain->SetBranchAddress("fatjet_sub_fHPD", fatjet_sub_fHPD, &b_fatjet_sub_fHPD);
   fChain->SetBranchAddress("fatjet_sub_fRBX", fatjet_sub_fRBX, &b_fatjet_sub_fRBX);
   fChain->SetBranchAddress("SC_n", &SC_n, &b_SC_n);
   fChain->SetBranchAddress("SC_truth", SC_truth, &b_SC_truth);
   fChain->SetBranchAddress("SC_E", SC_E, &b_SC_E);
   fChain->SetBranchAddress("SC_phi", SC_phi, &b_SC_phi);
   fChain->SetBranchAddress("SC_eta", SC_eta, &b_SC_eta);
   fChain->SetBranchAddress("ele_n", &ele_n, &b_ele_n);
   fChain->SetBranchAddress("ele_E", ele_E, &b_ele_E);
   fChain->SetBranchAddress("ele_Et", ele_Et, &b_ele_Et);
   fChain->SetBranchAddress("ele_p", ele_p, &b_ele_p);
   fChain->SetBranchAddress("ele_pt", ele_pt, &b_ele_pt);
   fChain->SetBranchAddress("ele_px", ele_px, &b_ele_px);
   fChain->SetBranchAddress("ele_py", ele_py, &b_ele_py);
   fChain->SetBranchAddress("ele_pz", ele_pz, &b_ele_pz);
   fChain->SetBranchAddress("ele_eta", ele_eta, &b_ele_eta);
   fChain->SetBranchAddress("ele_phi", ele_phi, &b_ele_phi);
   fChain->SetBranchAddress("ele_charge", ele_charge, &b_ele_charge);
   fChain->SetBranchAddress("ele_TrkChiNorm", ele_TrkChiNorm, &b_ele_TrkChiNorm);
   fChain->SetBranchAddress("ele_d0vtx", ele_d0vtx, &b_ele_d0vtx);
   fChain->SetBranchAddress("ele_d0bs", ele_d0bs, &b_ele_d0bs);
   fChain->SetBranchAddress("ele_sd0", ele_sd0, &b_ele_sd0);
   fChain->SetBranchAddress("ele_hits", ele_hits, &b_ele_hits);
   fChain->SetBranchAddress("ele_truth", ele_truth, &b_ele_truth);
   fChain->SetBranchAddress("ele_isECal", ele_isECal, &b_ele_isECal);
   fChain->SetBranchAddress("ele_isTracker", ele_isTracker, &b_ele_isTracker);
   fChain->SetBranchAddress("ele_ID", ele_ID, &b_ele_ID);
   fChain->SetBranchAddress("ele_ValidHitFirstPxlB", ele_ValidHitFirstPxlB, &b_ele_ValidHitFirstPxlB);
   fChain->SetBranchAddress("ele_TrkExpHitsInner", ele_TrkExpHitsInner, &b_ele_TrkExpHitsInner);
   fChain->SetBranchAddress("ele_HCalOverEm", ele_HCalOverEm, &b_ele_HCalOverEm);
   fChain->SetBranchAddress("ele_Dr03TkSumPt", ele_Dr03TkSumPt, &b_ele_Dr03TkSumPt);
   fChain->SetBranchAddress("ele_Dr04HCalSumEt", ele_Dr04HCalSumEt, &b_ele_Dr04HCalSumEt);
   fChain->SetBranchAddress("ele_Dr03HCalSumEt", ele_Dr03HCalSumEt, &b_ele_Dr03HCalSumEt);
   fChain->SetBranchAddress("ele_Dr04ECalSumEt", ele_Dr04ECalSumEt, &b_ele_Dr04ECalSumEt);
   fChain->SetBranchAddress("ele_Dr03ECalSumEt", ele_Dr03ECalSumEt, &b_ele_Dr03ECalSumEt);
   fChain->SetBranchAddress("ele_SigmaIetaIeta", ele_SigmaIetaIeta, &b_ele_SigmaIetaIeta);
   fChain->SetBranchAddress("ele_dEtaSCTrackAtVtx", ele_dEtaSCTrackAtVtx, &b_ele_dEtaSCTrackAtVtx);
   fChain->SetBranchAddress("ele_dPhiSCTrackAtVtx", ele_dPhiSCTrackAtVtx, &b_ele_dPhiSCTrackAtVtx);
   fChain->SetBranchAddress("ele_dr03HcalDepth1", ele_dr03HcalDepth1, &b_ele_dr03HcalDepth1);
   fChain->SetBranchAddress("ele_dr03HcalDepth2", ele_dr03HcalDepth2, &b_ele_dr03HcalDepth2);
   fChain->SetBranchAddress("ele_e2x5Max", ele_e2x5Max, &b_ele_e2x5Max);
   fChain->SetBranchAddress("ele_e5x5", ele_e5x5, &b_ele_e5x5);
   fChain->SetBranchAddress("ele_e1x5", ele_e1x5, &b_ele_e1x5);
   fChain->SetBranchAddress("ele_caloEt", ele_caloEt, &b_ele_caloEt);
   fChain->SetBranchAddress("ele_SCeta", ele_SCeta, &b_ele_SCeta);
   fChain->SetBranchAddress("ele_convdist", ele_convdist, &b_ele_convdist);
   fChain->SetBranchAddress("ele_convdcot", ele_convdcot, &b_ele_convdcot);
   fChain->SetBranchAddress("ele_convr", ele_convr, &b_ele_convr);
   fChain->SetBranchAddress("ele_fbrem", ele_fbrem, &b_ele_fbrem);
   fChain->SetBranchAddress("ele_trign", ele_trign, &b_ele_trign);
   fChain->SetBranchAddress("ele_trig", ele_trig, &b_ele_trig);
   fChain->SetBranchAddress("ele_SC", ele_SC, &b_ele_SC);
   fChain->SetBranchAddress("pfele_n", &pfele_n, &b_pfele_n);
   fChain->SetBranchAddress("pfele_p", pfele_p, &b_pfele_p);
   fChain->SetBranchAddress("pfele_E", pfele_E, &b_pfele_E);
   fChain->SetBranchAddress("pfele_Et", pfele_Et, &b_pfele_Et);
   fChain->SetBranchAddress("pfele_pt", pfele_pt, &b_pfele_pt);
   fChain->SetBranchAddress("pfele_px", pfele_px, &b_pfele_px);
   fChain->SetBranchAddress("pfele_py", pfele_py, &b_pfele_py);
   fChain->SetBranchAddress("pfele_pz", pfele_pz, &b_pfele_pz);
   fChain->SetBranchAddress("pfele_eta", pfele_eta, &b_pfele_eta);
   fChain->SetBranchAddress("pfele_phi", pfele_phi, &b_pfele_phi);
   fChain->SetBranchAddress("pfele_charge", pfele_charge, &b_pfele_charge);
   fChain->SetBranchAddress("pfele_truth", pfele_truth, &b_pfele_truth);
   fChain->SetBranchAddress("pfele_trign", pfele_trign, &b_pfele_trign);
   fChain->SetBranchAddress("pfele_trig", pfele_trig, &b_pfele_trig);
   fChain->SetBranchAddress("pfele_SC", pfele_SC, &b_pfele_SC);
   fChain->SetBranchAddress("muo_n", &muo_n, &b_muo_n);
   fChain->SetBranchAddress("muo_E", muo_E, &b_muo_E);
   fChain->SetBranchAddress("muo_Et", muo_Et, &b_muo_Et);
   fChain->SetBranchAddress("muo_p", muo_p, &b_muo_p);
   fChain->SetBranchAddress("muo_pt", muo_pt, &b_muo_pt);
   fChain->SetBranchAddress("muo_px", muo_px, &b_muo_px);
   fChain->SetBranchAddress("muo_py", muo_py, &b_muo_py);
   fChain->SetBranchAddress("muo_pz", muo_pz, &b_muo_pz);
   fChain->SetBranchAddress("muo_eta", muo_eta, &b_muo_eta);
   fChain->SetBranchAddress("muo_phi", muo_phi, &b_muo_phi);
   fChain->SetBranchAddress("muo_charge", muo_charge, &b_muo_charge);
   fChain->SetBranchAddress("muo_RelTrkIso", muo_RelTrkIso, &b_muo_RelTrkIso);
   fChain->SetBranchAddress("muo_TrkIso", muo_TrkIso, &b_muo_TrkIso);
   fChain->SetBranchAddress("muo_ECalIso", muo_ECalIso, &b_muo_ECalIso);
   fChain->SetBranchAddress("muo_HCalIso", muo_HCalIso, &b_muo_HCalIso);
   fChain->SetBranchAddress("muo_TrkIsoDep", muo_TrkIsoDep, &b_muo_TrkIsoDep);
   fChain->SetBranchAddress("muo_ECalIsoDep", muo_ECalIsoDep, &b_muo_ECalIsoDep);
   fChain->SetBranchAddress("muo_HCalIsoDep", muo_HCalIsoDep, &b_muo_HCalIsoDep);
   fChain->SetBranchAddress("muo_AllIso", muo_AllIso, &b_muo_AllIso);
   fChain->SetBranchAddress("muo_TrkChiNormCm", muo_TrkChiNormCm, &b_muo_TrkChiNormCm);
   fChain->SetBranchAddress("muo_TrkChiNormTk", muo_TrkChiNormTk, &b_muo_TrkChiNormTk);
   fChain->SetBranchAddress("muo_d0Cm", muo_d0Cm, &b_muo_d0Cm);
   fChain->SetBranchAddress("muo_d0Tk", muo_d0Tk, &b_muo_d0Tk);
   fChain->SetBranchAddress("muo_sd0Cm", muo_sd0Cm, &b_muo_sd0Cm);
   fChain->SetBranchAddress("muo_sd0Tk", muo_sd0Tk, &b_muo_sd0Tk);
   fChain->SetBranchAddress("muo_calocomp", muo_calocomp, &b_muo_calocomp);
   fChain->SetBranchAddress("muo_calotower_e", muo_calotower_e, &b_muo_calotower_e);
   fChain->SetBranchAddress("muo_prompttight", muo_prompttight, &b_muo_prompttight);
   fChain->SetBranchAddress("muo_hitsCm", muo_hitsCm, &b_muo_hitsCm);
   fChain->SetBranchAddress("muo_hitsTk", muo_hitsTk, &b_muo_hitsTk);
   fChain->SetBranchAddress("muo_truth", muo_truth, &b_muo_truth);
   fChain->SetBranchAddress("muo_trign", muo_trign, &b_muo_trign);
   fChain->SetBranchAddress("muo_trig", muo_trig, &b_muo_trig);
   fChain->SetBranchAddress("muo_ID", muo_ID, &b_muo_ID);
   fChain->SetBranchAddress("muo_ValidMuonHitsCm", muo_ValidMuonHitsCm, &b_muo_ValidMuonHitsCm);
   fChain->SetBranchAddress("muo_ValidTrackerHitsCm", muo_ValidTrackerHitsCm, &b_muo_ValidTrackerHitsCm);
   fChain->SetBranchAddress("muo_ValidPixelHitsCm", muo_ValidPixelHitsCm, &b_muo_ValidPixelHitsCm);
   fChain->SetBranchAddress("muo_ChambersMatched", muo_ChambersMatched, &b_muo_ChambersMatched);
   fChain->SetBranchAddress("muo_d0bsCm", muo_d0bsCm, &b_muo_d0bsCm);
   fChain->SetBranchAddress("muo_d0OriginCm", muo_d0OriginCm, &b_muo_d0OriginCm);
   fChain->SetBranchAddress("muo_dzbsCm", muo_dzbsCm, &b_muo_dzbsCm);
   fChain->SetBranchAddress("muo_TrackerLayersMeasCm", muo_TrackerLayersMeasCm, &b_muo_TrackerLayersMeasCm);
   fChain->SetBranchAddress("muo_TrackerLayersNotMeasCm", muo_TrackerLayersNotMeasCm, &b_muo_TrackerLayersNotMeasCm);
   fChain->SetBranchAddress("muo_Cocktail_pt", muo_Cocktail_pt, &b_muo_Cocktail_pt);
   fChain->SetBranchAddress("muo_Cocktail_phi", muo_Cocktail_phi, &b_muo_Cocktail_phi);
   fChain->SetBranchAddress("muo_Cocktail_eta", muo_Cocktail_eta, &b_muo_Cocktail_eta);
   Notify();
}

Bool_t TreeContent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeContent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeContent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeContent_cxx
