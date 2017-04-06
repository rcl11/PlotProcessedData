#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/BitManip.hh>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TColor.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <string>
#include <fstream>
#include <bitset>

#include "styles.h"


std::pair<std::string, std::string> ParseRunInfo(const std::string filename){
    std::string run_info = filename;
    run_info = run_info.erase(0,run_info.find("_r")+1);
    std::string run_number = run_info;
    run_number.erase(run_number.find("_"), run_number.size()+1);
    std::string subrun_number = run_info;
    subrun_number.erase(0,subrun_number.find("_s")+1);
    subrun_number.erase(subrun_number.find("_p"),subrun_number.size()-1);
    std::pair<std::string, std::string> info;
    info.first = run_number;
    info.second = subrun_number;
    return info;
}

void CreateRunPlots( const std::vector<std::string>& files, bool ntuple=true, std::string output_dir="plots/" )
{
    
  TH1D* hNHits_vs_run = new TH1D( "hNHits_vs_run", "Mean number of hits per run", files.size(), 0.0, files.size() );
  TH1D* hTotalQ_vs_run = new TH1D( "hTotalQ_vs_run", "Mean total charge per run", files.size(), 0.0, files.size() );
  TH1D* hNEvents_vs_run = new TH1D( "hNEvents_vs_run", "Number of events per run", files.size(), 0.0, files.size() );
  std::string start_run;
  std::string end_run;
  std::string postfix = "";
  if(ntuple) postfix = "_ntuple";
  for(unsigned i=0; i<files.size(); ++i){
      
      SetStyle();
      TCanvas* c1 = new TCanvas();
      
      //=============================================================================
      //Create a collection of histograns
      TH1D* hNHits = new TH1D( "hNHits", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hNHits_log = new TH1D( "hNHits_log", "Number of hits per event on log scale", 10000, 0.0, 10000.0 );
      TH1D* hTotalQ = new TH1D( "hTotalQ", "hTotalQ", 1000, 0.0, 10000.0 );
      TH1D* hNHits_pulseGT = new TH1D( "hNHits_pulseGT", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_pulseGT = new TH1D( "hTotalQ_pulseGT", "hTotalQ_pulseGT", 1000, 0.0, 10000.0 );
      TH1D* hNHits_N100M = new TH1D( "hNHits_N100M", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_N100M = new TH1D( "hTotalQ_N100M", "hTotalQ_N100M", 1000, 0.0, 10000.0 );
      TH1D* hNHits_N100H = new TH1D( "hNHits_N100H", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_N100H = new TH1D( "hTotalQ_N100H", "hTotalQ_N100H", 1000, 0.0, 10000.0 );
      TH1D* hNHits_N20 = new TH1D( "hNHits_N20", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_N20 = new TH1D( "hTotalQ_N20", "hTotalQ_N20", 1000, 0.0, 10000.0 );
      TH1D* hNHits_ESUMH = new TH1D( "hNHits_ESUMH", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_ESUMH = new TH1D( "hTotalQ_ESUMH", "hTotalQ_ESUMH", 1000, 0.0, 10000.0 );
      TH1D* hNHits_OWLEH = new TH1D( "hNHits_OWLEH", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ_OWLEH = new TH1D( "hTotalQ_OWLEH", "hTotalQ_OWLEH", 1000, 0.0, 10000.0 );
      TH1D* hfitValid = new TH1D( "hfitValid", "fitValid", 2, 0, 2 );
      TH1D* hitr = new TH1D( "hitr", "ITR", 100, 0, 1 );
      TH1D* hposx = new TH1D( "hposx", "Fitted x position", 160, -8000.0, 8000.0 );
      TH1D* hposy = new TH1D( "hposy", "Fitted y position", 160, -8000.0, 8000.0 );
      TH1D* hposz = new TH1D( "hposz", "Fitted z position", 160, -8000.0, 8000.0 );
      TH1D* hposR = new TH1D( "hposR", "Fitted R position", 80, 0.0, 8000.0 );
      TH2D* hposxy = new TH2D( "hposxy", "Fitted X vs Y", 160, -8000, 8000.0, 160, -8000, 8000.0 );
      TH2D* hposrz = new TH2D( "hposrz", "Fitted r vs z", 80, 0, 8000.0, 160, -8000, 8000.0 );
      TH1D* hrpmt = new TH1D( "hrpmt", "Hit PMT R", 100, 0.0, 10000.0 );
      TH1D* hxpmt = new TH1D( "hxpmt", "Hit PMT X", 100, -10000.0, 10000.0 );
      TH1D* hypmt = new TH1D( "hypmt", "Hit PMT Y", 100, -10000.0, 10000.0 );
      TH1D* hzpmt = new TH1D( "hzpmt", "Hit PMT Z", 100, -10000.0, 10000.0 );
      TH1D* htpmt = new TH1D( "htpmt", "Hit PMT raw time", 1000, -500.0, 500.0 );
      
      TFile *f = new TFile(files[i].c_str());
      
      int n_events=0;
      
      //=============================================================================
      //Fill histograms using info from ratds or ntuple
      if(ntuple) {
          TTree *t1 = (TTree*)f->Get("output");
          Int_t nhits;
          Double_t charge;
          ULong64_t flag;
          ULong64_t applied_flag;
          Double_t posx, posy, posz;
          bool fit_valid;
          Double_t itr;
          Int_t triggerWord;
          t1->SetBranchAddress("nhits",&nhits);
          t1->SetBranchAddress("q",&charge);
          t1->SetBranchAddress("dcFlagged",&flag);
          t1->SetBranchAddress("dcApplied",&applied_flag);
          t1->SetBranchAddress("posx",&posx);
          t1->SetBranchAddress("posy",&posy);
          t1->SetBranchAddress("posz",&posz);
          t1->SetBranchAddress("fitValid",&fit_valid);
          t1->SetBranchAddress("triggerWord",&triggerWord);
          t1->SetBranchAddress("itr",&itr);

          Long64_t nentries = t1->GetEntries();
          for (Long64_t j=0;j<nentries;j++) {
             t1->GetEntry(j);
             //analysis_mask
             //bool dataclean = !(flag & 0b111111111111110);
             //analysis_mask without tpmuonfollowercut-short
             bool dataclean = !(flag & 0b11111111111110);
             //analysis_mask
             //bool compatibility_cut = (applied_flag & 0b111111111111110) == 0b111111111111110;
             //analysis_mask without tpmuonfollowercut-short
             bool compatibility_cut = (applied_flag & 0b11111111111110) == 0b11111111111110;
             if(dataclean && compatibility_cut) {
                hNHits->Fill(nhits);
                hNHits_log->Fill(nhits);
                hTotalQ->Fill(charge);
                //Fill some nhits and total q plots for different triggers fired
                std::bitset<32> bits = std::bitset<32>(triggerWord);
                //Numbers to test taken from RAT documentation
                if(bits.test(10)){
                    hNHits_pulseGT->Fill(nhits);
                    hTotalQ_pulseGT->Fill(charge);
                }
                if(bits.test(1)){
                    hNHits_N100M->Fill(nhits);
                    hTotalQ_N100M->Fill(charge);
                }
                if(bits.test(2)){
                    hNHits_N100H->Fill(nhits);
                    hTotalQ_N100H->Fill(charge);
                }
                if(bits.test(3)){
                    hNHits_N20->Fill(nhits);
                    hTotalQ_N20->Fill(charge);
                }
                if(bits.test(6)){
                    hNHits_ESUMH->Fill(nhits);
                    hTotalQ_ESUMH->Fill(charge);
                }
                if(bits.test(9)){
                    hNHits_OWLEH->Fill(nhits);
                    hTotalQ_OWLEH->Fill(charge);
                }
                hfitValid->Fill(fit_valid); 
                if(fit_valid){
                    hposxy->Fill(posx,posy);
                    hposrz->Fill(sqrt(posx*posx + posy*posy), posz);
                    hposx->Fill(posx);
                    hposy->Fill(posy);
                    hposz->Fill(posz);
                    hposR->Fill(sqrt(posx*posx + posy*posy + posz*posz));
                    hitr->Fill(itr);
                }
                n_events++;
             }
          }
      } else {
          RAT::DU::DSReader dsReader( files[i] );
          const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
          //Add data cleaning cut. Note should possibly be changed to "analysis_mask" once bug is fixed with tpmuonfollowercut
          //ULong64_t rDataCleaningWord = RAT::GetDataCleaningWord( "analysis_mask" );
          //Temporary mask removing tpmuonfollowercut-short which has a bug for this production
          ULong64_t rDataCleaningWord = RAT::GetDataCleaningWord( "analysis_mask_temp" );
          for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
          {
              const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
              for( size_t iEv = 0; iEv < rDS.GetEVCount(); iEv++ )
              {   
                  RAT::DS::EV rEV = rDS.GetEV(iEv);
                  const RAT::DS::Meta& rMeta = dsReader.GetMeta();
                  if(RAT::EventIsClean( rEV, rMeta, rDataCleaningWord )){
                      n_events++;
                      hNHits->Fill( rEV.GetNhits() );
                      hNHits_log->Fill( rEV.GetNhits() );
                      hTotalQ->Fill( rEV.GetTotalCharge() );
                      //Fill some nhits and total q plots for different triggers fired
                      //std::cout << std::bitset<32>(rEV.GetTrigType())/*.to_string()*/ << std::endl;
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::PulseGT)){
                          hNHits_pulseGT->Fill(rEV.GetNhits()); 
                          hTotalQ_pulseGT->Fill( rEV.GetTotalCharge() );
                      }
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N100Med)){
                          hNHits_N100M->Fill(rEV.GetNhits()); 
                          hTotalQ_N100M->Fill( rEV.GetTotalCharge() );
                      }
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N100High)){
                          hNHits_N100H->Fill(rEV.GetNhits()); 
                          hTotalQ_N100H->Fill( rEV.GetTotalCharge() );
                      }
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N20)){
                          hNHits_N20->Fill(rEV.GetNhits()); 
                          hTotalQ_N20->Fill( rEV.GetTotalCharge() );
                      }
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::ESHigh)){
                          hNHits_ESUMH->Fill(rEV.GetNhits()); 
                          hTotalQ_ESUMH->Fill( rEV.GetTotalCharge() );
                      }
                      if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::OWLESHigh)){
                          hNHits_OWLEH->Fill(rEV.GetNhits()); 
                          hTotalQ_OWLEH->Fill( rEV.GetTotalCharge() );
                      }
                      
                      if(rEV.FitResultExists("waterFitter") && rEV.GetFitResult("waterFitter").GetValid()){
                          RAT::DS::FitVertex rvertex = rEV.GetFitResult("waterFitter").GetVertex(0);
                          if( rvertex.ContainsPosition() && rvertex.ValidPosition() ) {
                              hfitValid->Fill(1.);
                              hposxy->Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                              hposrz->Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y()), 
                                    rvertex.GetPosition().Z());
                              hposx->Fill(rvertex.GetPosition().X());
                              hposy->Fill(rvertex.GetPosition().Y());
                              hposz->Fill(rvertex.GetPosition().Z());
                              hposR->Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y() + rvertex.GetPosition().Z()*rvertex.GetPosition().Z()) );
                              hitr->Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"));
                          } else hfitValid->Fill(0.);
                      } else hfitValid->Fill(0.);
                      RAT::DS::CalPMTs& calpmts = rEV.GetCalPMTs();
                      //NOTE: i think we should return to GetCount() and GetPMT() after next reprocess (assuming we want to include the "inward" PMTs (normal + HQE) not just "normal" here)
                      for(unsigned int ipmt=0;ipmt<calpmts.GetNormalCount();ipmt++){
                          TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetNormalPMT(ipmt).GetID());
                          double pmt_r = pmtpos.Mag();
                          hrpmt->Fill(pmt_r);
                          hxpmt->Fill(pmtpos.X());
                          hypmt->Fill(pmtpos.Y());
                          hzpmt->Fill(pmtpos.Z());
                          htpmt->Fill((calpmts.GetNormalPMT(ipmt)).GetTime());
                      }
                  }
              }
          }
      }
      
      std::pair<std::string,std::string> run_info = ParseRunInfo(files[i]);
      std::string outname = run_info.first + "_" + run_info.second; 
      
      //=============================================================================
      //Make some plots
      hNHits->SetFillStyle(1001);
      hNHits->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits->GetYaxis()->SetTitle( "Events" );
      hNHits->GetXaxis()->SetTitle( "Number of hits" );
      hNHits->Draw("hist");
      TLatex *title_latex = new TLatex();
      title_latex->SetNDC();
      title_latex->SetTextSize(0.04);
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_"+outname+postfix+".pdf").c_str());
      c1->SetLogx();
      c1->SetLogy();
      c1->SaveAs((output_dir+"nhits_log_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_log_"+outname+postfix+".pdf").c_str());
      c1->SetLogx(0);
      c1->SetLogy(0);
      c1->Clear();
      
      hTotalQ->SetFillStyle(1001);
      hTotalQ->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ->GetYaxis()->SetTitle( "Events" );
      hTotalQ->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_pulseGT->SetFillStyle(1001);
      hNHits_pulseGT->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_pulseGT->GetYaxis()->SetTitle( "Events" );
      hNHits_pulseGT->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_pulseGT->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_pulseGT_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_pulseGT_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_pulseGT->SetFillStyle(1001);
      hTotalQ_pulseGT->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_pulseGT->GetYaxis()->SetTitle( "Events" );
      hTotalQ_pulseGT->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_pulseGT->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_pulseGT_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_pulseGT_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_N100M->SetFillStyle(1001);
      hNHits_N100M->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_N100M->GetYaxis()->SetTitle( "Events" );
      hNHits_N100M->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_N100M->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_N100M_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_N100M_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_N100M->SetFillStyle(1001);
      hTotalQ_N100M->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_N100M->GetYaxis()->SetTitle( "Events" );
      hTotalQ_N100M->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_N100M->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_N100M_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_N100M_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_N100H->SetFillStyle(1001);
      hNHits_N100H->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_N100H->GetYaxis()->SetTitle( "Events" );
      hNHits_N100H->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_N100H->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_N100H_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_N100H_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_N100H->SetFillStyle(1001);
      hTotalQ_N100H->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_N100H->GetYaxis()->SetTitle( "Events" );
      hTotalQ_N100H->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_N100H->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_N100H_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_N100H_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_N20->SetFillStyle(1001);
      hNHits_N20->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_N20->GetYaxis()->SetTitle( "Events" );
      hNHits_N20->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_N20->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_N20_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_N20_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_N20->SetFillStyle(1001);
      hTotalQ_N20->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_N20->GetYaxis()->SetTitle( "Events" );
      hTotalQ_N20->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_N20->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_N20_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_N20_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_ESUMH->SetFillStyle(1001);
      hNHits_ESUMH->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_ESUMH->GetYaxis()->SetTitle( "Events" );
      hNHits_ESUMH->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_ESUMH->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_ESUMH_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_ESUMH_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_ESUMH->SetFillStyle(1001);
      hTotalQ_ESUMH->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_ESUMH->GetYaxis()->SetTitle( "Events" );
      hTotalQ_ESUMH->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_ESUMH->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_ESUMH_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_ESUMH_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_OWLEH->SetFillStyle(1001);
      hNHits_OWLEH->SetFillColor(TColor::GetColor(220, 24, 24));
      hNHits_OWLEH->GetYaxis()->SetTitle( "Events" );
      hNHits_OWLEH->GetXaxis()->SetTitle( "Number of hits" );
      hNHits_OWLEH->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"nhits_OWLEH_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"nhits_OWLEH_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hTotalQ_OWLEH->SetFillStyle(1001);
      hTotalQ_OWLEH->SetFillColor(TColor::GetColor(142, 24, 220));
      hTotalQ_OWLEH->GetYaxis()->SetTitle( "Events" );
      hTotalQ_OWLEH->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ_OWLEH->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"totalQ_OWLEH_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"totalQ_OWLEH_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hNHits_vs_run->Fill(i, hNHits_log->GetMean());
      hNHits_vs_run->SetBinError(i+1,hNHits_log->GetMeanError());
      std::string bin_label = run_info.first;
      if(i==0) start_run = bin_label;
      if(i==files.size()-1) end_run = bin_label;
      hNHits_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
      hTotalQ_vs_run->Fill(i, hTotalQ->GetMean());
      hTotalQ_vs_run->SetBinError(i+1,hTotalQ->GetMeanError());
      hTotalQ_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
      hNEvents_vs_run->Fill(i, n_events);
      hNEvents_vs_run->SetBinError(i+1,sqrt(n_events));
      hNEvents_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
      
      hposx->GetYaxis()->SetTitle( "Events" );
      hposx->GetXaxis()->SetTitle( "X Position (mm)" );
      hposx->SetFillStyle(1001);
      hposx->SetFillColor(TColor::GetColor(220, 24, 70));
      hposx->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"posx_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"posx_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hposy->GetYaxis()->SetTitle( "Events" );
      hposy->GetXaxis()->SetTitle( "Y Position (mm)" );
      hposy->SetFillStyle(1001);
      hposy->SetFillColor(TColor::GetColor(24,220,57));
      hposy->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"posy_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"posy_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hposz->GetYaxis()->SetTitle( "Events" );
      hposz->GetXaxis()->SetTitle( "Z Position (mm)" );
      hposz->SetFillStyle(1001);
      hposz->SetFillColor(TColor::GetColor(24,113,220));
      hposz->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"posz_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"posz_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hposR->GetYaxis()->SetTitle( "Events" );
      hposR->GetXaxis()->SetTitle( "R Position (mm)" );
      hposR->SetFillStyle(1001);
      hposR->SetFillColor(kYellow);
      hposR->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"posR_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"posR_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hfitValid->GetYaxis()->SetTitle( "Events" );
      hfitValid->GetXaxis()->SetTitle( "Fit is valid" );
      hfitValid->SetFillStyle(1001);
      hfitValid->SetFillColor(TColor::GetColor(255, 117, 250));
      hfitValid->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"fitValid_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"fitValid_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      hitr->GetYaxis()->SetTitle( "Events" );
      hitr->GetXaxis()->SetTitle( "ITR" );
      hitr->SetFillStyle(1001);
      hitr->SetFillColor(TColor::GetColor(222, 157, 59));
      hitr->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c1->Update();
      c1->SaveAs((output_dir+"itr_"+outname+postfix+".png").c_str());
      c1->SaveAs((output_dir+"itr_"+outname+postfix+".pdf").c_str());
      c1->Clear();
      
      
      if(!ntuple){
          hrpmt->GetYaxis()->SetTitle( "Events" );
          hrpmt->GetXaxis()->SetTitle( "hit PMT R Position (mm)" );
          hrpmt->SetFillStyle(1001);
          hrpmt->SetFillColor(kOrange);
          hrpmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c1->Update();
          c1->SaveAs((output_dir+"rpmt_"+outname+postfix+".png").c_str());
          c1->SaveAs((output_dir+"rpmt_"+outname+postfix+".pdf").c_str());
          c1->Clear();
          
          hxpmt->GetYaxis()->SetTitle( "Events" );
          hxpmt->GetXaxis()->SetTitle( "hit PMT X Position (mm)" );
          hxpmt->SetFillStyle(1001);
          hxpmt->SetFillColor(TColor::GetColor(220, 24, 70));
          hxpmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c1->Update();
          c1->SaveAs((output_dir+"xpmt_"+outname+postfix+".png").c_str());
          c1->SaveAs((output_dir+"xpmt_"+outname+postfix+".pdf").c_str());
          c1->Clear();
          
          hypmt->GetYaxis()->SetTitle( "Events" );
          hypmt->GetXaxis()->SetTitle( "hit PMT Y Position (mm)" );
          hypmt->SetFillStyle(1001);
          hypmt->SetFillColor(TColor::GetColor(24,220,57));
          hypmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c1->Update();
          c1->SaveAs((output_dir+"ypmt_"+outname+postfix+".png").c_str());
          c1->SaveAs((output_dir+"ypmt_"+outname+postfix+".pdf").c_str());
          c1->Clear();
          
          hzpmt->GetYaxis()->SetTitle( "Events" );
          hzpmt->GetXaxis()->SetTitle( "hit PMT Z Position (mm)" );
          hzpmt->SetFillStyle(1001);
          hzpmt->SetFillColor(TColor::GetColor(24,113,220));
          hzpmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c1->Update();
          c1->SaveAs((output_dir+"zpmt_"+outname+postfix+".png").c_str());
          c1->SaveAs((output_dir+"zpmt_"+outname+postfix+".pdf").c_str());
          c1->Clear();
          
          htpmt->GetYaxis()->SetTitle( "Events" );
          htpmt->GetXaxis()->SetTitle( "hit PMT raw time (ns)" );
          htpmt->SetFillStyle(1001);
          htpmt->SetFillColor(TColor::GetColor(76, 220, 215));
          htpmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c1->Update();
          c1->SaveAs((output_dir+"tpmt_"+outname+postfix+".png").c_str());
          c1->SaveAs((output_dir+"tpmt_"+outname+postfix+".pdf").c_str());
          c1->Clear();
      }
      
      TCanvas* c2 = new TCanvas("c2","c2",600,500);
      hposxy->GetYaxis()->SetTitle( "Y position (mm)" );
      hposxy->GetXaxis()->SetTitle( "X position (mm)" );
      hposxy->SetMarkerStyle(20);
      hposxy->Draw("colzsame");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c2->SetRightMargin(0.15);
      c2->Update();
      c2->SaveAs((output_dir+"posxy_"+outname+postfix+".png").c_str());
      c2->SaveAs((output_dir+"posxy_"+outname+postfix+".pdf").c_str());
      c2->Clear();
      
      TCanvas* c3 = new TCanvas("c3","c3",600,500);
      hposrz->GetYaxis()->SetTitle( "Z position (mm)" );
      hposrz->GetXaxis()->SetTitle( "r = #sqrt{(x^2 + y^2)} position (mm)" );
      hposrz->SetMarkerStyle(20);
      hposrz->Draw("colzsame");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c3->SetRightMargin(0.15);
      c3->Update();
      c3->SaveAs((output_dir+"posrz_"+outname+postfix+".png").c_str());
      c3->SaveAs((output_dir+"posrz_"+outname+postfix+".pdf").c_str());
      c3->Clear();
      
      delete hNHits, hNHits_pulseGT, hNHits_N100M, hNHits_N100H, hNHits_N20, hNHits_ESUMH, hNHits_OWLEH;
      delete hTotalQ, hTotalQ_pulseGT, hTotalQ_N100M, hTotalQ_N100H, hTotalQ_N20, hTotalQ_ESUMH, hTotalQ_OWLEH;
      delete hposxy, hposrz, hposx, hposy, hposz, hposR, hrpmt, hxpmt, hypmt, hzpmt, htpmt, hfitValid;
      delete c1;
      delete c2;
      delete c3;
  }
  SetStyle();
  TCanvas* c100 = new TCanvas("c100","c100",1000,400);
  hNHits_vs_run->SetMarkerColor(TColor::GetColor(220, 24, 24));
  hNHits_vs_run->SetMarkerStyle(20);
  hNHits_vs_run->SetLineColor(TColor::GetColor(220, 24, 24));
  hNHits_vs_run->GetYaxis()->SetTitle( "Mean #Hits" );
  hNHits_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNHits_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNHits_vs_run->Draw("PE1");
  c100->SetLeftMargin(0.1);
  c100->SetBottomMargin(0.18);
  c100->Update();
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hTotalQ_vs_run->SetMarkerColor(TColor::GetColor(142, 24, 220));
  hTotalQ_vs_run->SetMarkerStyle(20);
  hTotalQ_vs_run->SetLineColor(TColor::GetColor(142, 24, 220));
  hTotalQ_vs_run->GetYaxis()->SetTitle( "Mean TotalQ" );
  hTotalQ_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hTotalQ_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hTotalQ_vs_run->Draw("PE1");
  c100->SetLeftMargin(0.1);
  c100->SetBottomMargin(0.18);
  c100->Update();
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNEvents_vs_run->SetMarkerColor(kGreen);
  hNEvents_vs_run->SetMarkerStyle(20);
  hNEvents_vs_run->SetLineColor(kGreen);
  hNEvents_vs_run->GetYaxis()->SetTitle( "# Events" );
  hNEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNEvents_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNEvents_vs_run->Draw("PE1");
  c100->SetLeftMargin(0.1);
  c100->SetBottomMargin(0.18);
  c100->Update();
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  delete c100;
  delete hNHits_vs_run;
  delete hTotalQ_vs_run;
  delete hNEvents_vs_run;
 
}

void makeRunPlots(std::string filelist, std::string output_dir = "plots/")
{
  std::ifstream infile(filelist.c_str());
  std::vector<std::string> files;
  std::string line;
  while (std::getline(infile, line)){
      std::istringstream iss(line);
      files.push_back(iss.str());
  }
  //Possible to use both ntuple files or ratds files. Assumes that the txt file containing the filelist follows my naming convention
  bool ntuple=false;
  if(files[0].find("ntuple")!=std::string::npos) ntuple=true;
  CreateRunPlots(files,ntuple,output_dir);
  return;
}
