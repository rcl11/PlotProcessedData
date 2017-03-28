#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Meta.hh>

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
      TH1D* hNHits = new TH1D( "hNHits", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ = new TH1D( "hTotalQ", "hTotalQ", 1000, 0.0, 10000.0 );
      TH1D* hposx = new TH1D( "hposx", "hposx", 80, 0.0, 8000.0 );
      TH1D* hposy = new TH1D( "hposy", "hposy", 80, 0.0, 8000.0 );
      TH1D* hposz = new TH1D( "hposz", "hposz", 80, 0.0, 8000.0 );
      TH1D* hposR = new TH1D( "hposR", "hposR", 80, 0.0, 8000.0 );
      TH2D* hposxy = new TH2D( "hposxy", "hposxy", 160, -8000, 8000.0, 160, -8000, 8000.0 );
      TH2D* hposrz = new TH2D( "hposrz", "hposrz", 80, 0, 8000.0, 160, -8000, 8000.0 );
      TH1D* hrpmt = new TH1D( "hrpmt", "hrpmt", 80, 0.0, 8000.0 );
      
      TFile *f = new TFile(files[i].c_str());
      TTree *t1 = (TTree*)f->Get("output");
      
      int n_events=0;
      if(ntuple) {
          Int_t nhits;
          Double_t charge;
          ULong64_t flag;
          ULong64_t applied_flag;
          Double_t posx, posy, posz;
          bool fit_valid;
          t1->SetBranchAddress("nhits",&nhits);
          t1->SetBranchAddress("q",&charge);
          t1->SetBranchAddress("dcFlagged",&flag);
          t1->SetBranchAddress("dcApplied",&applied_flag);
          t1->SetBranchAddress("posx",&posx);
          t1->SetBranchAddress("posy",&posy);
          t1->SetBranchAddress("posz",&posz);
          t1->SetBranchAddress("fitValid",&fit_valid);

          Long64_t nentries = t1->GetEntries();
          for (Long64_t j=0;j<nentries;j++) {
             t1->GetEntry(j);
             //Playing around with data cleaning cut here. I believe we need analysis_mask once new processing is ready
             //Add data cleaning cut. Corresponds to "default" but without prescale, calculated using Morgan's dcflags.py tool 
             //bool dataclean = !(flag & 0b10011111111111110);
             //analysis_mask
             //bool dataclean = !(flag & 0b111111111111110);
             //analysis_mask without tpmuonfollowercut-short
             bool dataclean = !(flag & 0b11111111111110);
             //bool compatibility_cut = (applied_flag & 0b10011111111111110) == 0b10011111111111110;
             //analysis_mask
             //bool compatibility_cut = (applied_flag & 0b111111111111110) == 0b111111111111110;
             //analysis_mask without tpmuonfollowercut-short
             bool compatibility_cut = (applied_flag & 0b11111111111110) == 0b11111111111110;
             if(dataclean && compatibility_cut) {
                hNHits->Fill(nhits);
                hTotalQ->Fill(charge);
                if(fit_valid){
                    hposxy->Fill(posx,posy);
                    hposrz->Fill(sqrt(posx*posx + posy*posy), posz);
                    hposx->Fill(posx);
                    hposy->Fill(posy);
                    hposz->Fill(posz);
                    hposR->Fill(sqrt(posx*posx + posy*posy + posz*posz));
                }
                n_events++;
             }
          }
      } else {
          RAT::DU::DSReader dsReader( files[i] );
          const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
          //Add data cleaning cut. Note should possibly be changed to "analysis_mask", however I believe this requires reprocessing first?
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
                      hTotalQ->Fill( rEV.GetTotalCharge() );
                      if(rEV.FitResultExists("partialWaterFitter") && rEV.GetFitResult("partialWaterFitter").GetValid()){
                          RAT::DS::FitVertex rvertex = rEV.GetFitResult("partialWaterFitter").GetVertex(0);
                          if( rvertex.ContainsPosition() && rvertex.ValidPosition() ) {
                              hposxy->Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                              hposrz->Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y()), 
                                    rvertex.GetPosition().Z());
                              hposx->Fill(rvertex.GetPosition().X());
                              hposy->Fill(rvertex.GetPosition().Y());
                              hposz->Fill(rvertex.GetPosition().Y());
                              hposR->Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y() + rvertex.GetPosition().Z()*rvertex.GetPosition().Z()) );
                          }
                      }
                      RAT::DS::CalPMTs& calpmts = rEV.GetCalPMTs();
                      for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++){
                          TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                          double pmt_r = pmtpos.Mag();
                          hrpmt->Fill(pmt_r);
                      }
                  }
              }
          }
      }
      
      std::pair<std::string,std::string> run_info = ParseRunInfo(files[i]);
      std::string outname = run_info.first + "_" + run_info.second; 
      //Make some plots
      hNHits->SetFillStyle(1001);
      hNHits->SetFillColor(TColor::GetColor(9, 114, 251));
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
      c1->Clear();
      SetStyle();
      
      TCanvas* c2 = new TCanvas();
      hTotalQ->SetFillStyle(1001);
      hTotalQ->SetFillColor(kMagenta);
      hTotalQ->GetYaxis()->SetTitle( "Events" );
      hTotalQ->GetXaxis()->SetTitle( "Total Charge" );
      hTotalQ->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c2->Update();
      c2->SaveAs((output_dir+"totalQ_"+outname+postfix+".png").c_str());
      c2->SaveAs((output_dir+"totalQ_"+outname+postfix+".pdf").c_str());
      c2->Clear();
      
      hNHits_vs_run->Fill(i, hNHits->GetMean());
      hNHits_vs_run->SetBinError(i+1,hNHits->GetMeanError());
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
      
      TCanvas* c3 = new TCanvas("c3","c3",600,500);
      hposxy->GetYaxis()->SetTitle( "Y position (mm)" );
      hposxy->GetXaxis()->SetTitle( "X position (mm)" );
      hposxy->SetMarkerStyle(20);
      hposxy->Draw("colzsame");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c3->SetRightMargin(0.15);
      c3->Update();
      c3->SaveAs((output_dir+"posxy_"+outname+postfix+".png").c_str());
      c3->SaveAs((output_dir+"posxy_"+outname+postfix+".pdf").c_str());
      c3->Clear();
      
      TCanvas* c4 = new TCanvas("c4","c4",600,500);
      hposrz->GetYaxis()->SetTitle( "Z position (mm)" );
      hposrz->GetXaxis()->SetTitle( "R position (mm)" );
      hposrz->SetMarkerStyle(20);
      hposrz->Draw("colzsame");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c4->SetRightMargin(0.15);
      c4->Update();
      c4->SaveAs((output_dir+"posrz_"+outname+postfix+".png").c_str());
      c4->SaveAs((output_dir+"posrz_"+outname+postfix+".pdf").c_str());
      c4->Clear();
      
      TCanvas* c5 = new TCanvas("c5","c5",600,500);
      hposx->GetYaxis()->SetTitle( "Events" );
      hposx->GetXaxis()->SetTitle( "X Position (mm)" );
      hposx->SetFillStyle(1001);
      hposx->SetFillColor(kGreen);
      hposx->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c5->Update();
      c5->SaveAs((output_dir+"posx_"+outname+postfix+".png").c_str());
      c5->SaveAs((output_dir+"posx_"+outname+postfix+".pdf").c_str());
      c5->Clear();
      
      TCanvas* c6 = new TCanvas("c6","c6",600,500);
      hposy->GetYaxis()->SetTitle( "Events" );
      hposy->GetXaxis()->SetTitle( "Y Position (mm)" );
      hposy->SetFillStyle(1001);
      hposy->SetFillColor(kGreen);
      hposy->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c6->Update();
      c6->SaveAs((output_dir+"posy_"+outname+postfix+".png").c_str());
      c6->SaveAs((output_dir+"posy_"+outname+postfix+".pdf").c_str());
      c6->Clear();
      
      TCanvas* c7 = new TCanvas("c7","c7",600,500);
      hposz->GetYaxis()->SetTitle( "Events" );
      hposz->GetXaxis()->SetTitle( "Z Position (mm)" );
      hposz->SetFillStyle(1001);
      hposz->SetFillColor(kCyan);
      hposz->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c7->Update();
      c7->SaveAs((output_dir+"posz_"+outname+postfix+".png").c_str());
      c7->SaveAs((output_dir+"posz_"+outname+postfix+".pdf").c_str());
      c7->Clear();
      
      TCanvas* c8 = new TCanvas("c8","c8",600,500);
      hposR->GetYaxis()->SetTitle( "Events" );
      hposR->GetXaxis()->SetTitle( "R Position (mm)" );
      hposR->SetFillStyle(1001);
      hposR->SetFillColor(kYellow);
      hposR->Draw("hist");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c8->Update();
      c8->SaveAs((output_dir+"posR_"+outname+postfix+".png").c_str());
      c8->SaveAs((output_dir+"posR_"+outname+postfix+".pdf").c_str());
      c8->Clear();
      
      TCanvas* c9 = new TCanvas("c9","c9",600,500);
      if(!ntuple){
          hrpmt->GetYaxis()->SetTitle( "Events" );
          hrpmt->GetXaxis()->SetTitle( "PMT R Position (mm)" );
          hrpmt->SetFillStyle(1001);
          hrpmt->SetFillColor(kOrange);
          hrpmt->Draw("hist");
          title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
          c9->Update();
          c9->SaveAs((output_dir+"rpmt_"+outname+postfix+".png").c_str());
          c9->SaveAs((output_dir+"rpmt_"+outname+postfix+".pdf").c_str());
          c9->Clear();
      }
      
      delete hNHits;
      delete hTotalQ;
      delete hposxy;
      delete hposrz;
      delete hposx;
      delete hposy;
      delete hposz;
      delete hposR;
      delete hrpmt;
      delete c1;
      delete c2;
      delete c3;
      delete c4;
      delete c5;
      delete c6;
      delete c7;
      delete c8;
      delete c9;
  }
  SetStyle();
  TCanvas* c10 = new TCanvas("c10","c10",1000,400);
  hNHits_vs_run->SetMarkerColor(TColor::GetColor(9, 114, 251));
  hNHits_vs_run->SetMarkerStyle(20);
  hNHits_vs_run->SetLineColor(TColor::GetColor(9, 114, 251));
  hNHits_vs_run->GetYaxis()->SetTitle( "Mean #Hits" );
  hNHits_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNHits_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNHits_vs_run->Draw("PE1");
  c10->SetLeftMargin(0.1);
  c10->SetBottomMargin(0.18);
  c10->Update();
  c10->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c10->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  SetStyle();
  TCanvas* c11 = new TCanvas("c11","c11",1000,400);
  hTotalQ_vs_run->SetMarkerColor(kMagenta);
  hTotalQ_vs_run->SetMarkerStyle(20);
  hTotalQ_vs_run->SetLineColor(kMagenta);
  hTotalQ_vs_run->GetYaxis()->SetTitle( "Mean TotalQ" );
  hTotalQ_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hTotalQ_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hTotalQ_vs_run->Draw("PE1");
  c11->SetLeftMargin(0.1);
  c11->SetBottomMargin(0.18);
  c11->Update();
  c11->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c11->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  SetStyle();
  TCanvas* c12 = new TCanvas("c12","c12",1000,400);
  hNEvents_vs_run->SetMarkerColor(kRed);
  hNEvents_vs_run->SetMarkerStyle(20);
  hNEvents_vs_run->SetLineColor(kRed);
  hNEvents_vs_run->GetYaxis()->SetTitle( "# Events" );
  hNEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNEvents_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNEvents_vs_run->Draw("PE1");
  c12->SetLeftMargin(0.1);
  c12->SetBottomMargin(0.18);
  c12->Update();
  c12->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c12->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  delete c10;
  delete c11;
  delete c12;
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
