#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DataCleaningUtility.hh>
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

void CreateRunPlots( const std::vector<std::string>& files, bool ntuple=true )
{
    
  TH1D* hNHits_vs_run = new TH1D( "hNHits_vs_run", "Mean number of hits per run", files.size(), 0.0, files.size() );
  TH1D* hTotalQ_vs_run = new TH1D( "hTotalQ_vs_run", "Mean total charge per run", files.size(), 0.0, files.size() );
  std::string start_run;
  std::string end_run;
  std::string postfix = "";
  if(ntuple) postfix = "_ntuple";
  for(unsigned i=0; i<files.size(); ++i){
      SetStyle();
      TCanvas* c1 = new TCanvas();
      TH1D* hNHits = new TH1D( "hNHits", "Number of hits per event", 100, 0.0, 100.0 );
      TH1D* hTotalQ = new TH1D( "hTotalQ", "hTotalQ", 500, 0.0, 5000.0 );
      TH2D* hposxy = new TH2D( "hposxy", "hposxy", 200, 0.0, 1000.0, 200, 0.0, 1000.0 );
      
      TFile *f = new TFile(files[i].c_str());
      TTree *t1 = (TTree*)f->Get("output");
      
      if(ntuple) {
          Int_t nhits;
          Double_t charge;
          ULong64_t flag;
          ULong64_t applied_flag;
          Double_t posx, posy;
          bool fit_valid;
          t1->SetBranchAddress("nhits",&nhits);
          t1->SetBranchAddress("q",&charge);
          t1->SetBranchAddress("dcFlagged",&flag);
          t1->SetBranchAddress("dcApplied",&applied_flag);
          t1->SetBranchAddress("posx",&posx);
          t1->SetBranchAddress("posy",&posy);
          t1->SetBranchAddress("fitValid",&fit_valid);

          Long64_t nentries = t1->GetEntries();
          for (Long64_t j=0;j<nentries;j++) {
             t1->GetEntry(j);
             //Add data cleaning cut. Corresponds to "default_apply", calculated using Morgan's dcflags.py tool 
             bool dataclean = !(flag & 0b10011111111111111);
             //Add data cleaning cut. Corresponds to "default", calculated using Morgan's dcflags.py tool 
             //bool dataclean = !(flag & 0b10000011111111111);
             bool compatibility_cut = applied_flag & 0b10000011111111111; 
             if(dataclean && compatibility_cut) {
                hNHits->Fill(nhits);
                hTotalQ->Fill(charge);
                if(fit_valid){
                    hposxy->Fill(posx,posy);
                }
             }
          }
      } else {
          RAT::DU::DSReader dsReader( files[0] );
          //Add data cleaning cut. Note should possibly be changed to "analysis_mask", however I believe this requires reprocessing first?
          ULong64_t rDataCleaningWord = RAT::GetDataCleaningWord( "default_apply" );
          for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
          {
              const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
              for( size_t iEv = 0; iEv < rDS.GetEVCount(); iEv++ )
              {
                  RAT::DS::EV rEV = rDS.GetEV(iEv);
                  const RAT::DS::Meta& rMeta = dsReader.GetMeta();
                  if(RAT::EventIsClean( rEV, rMeta, rDataCleaningWord )){
                      hNHits->Fill( rEV.GetNhits() );
                      hTotalQ->Fill( rEV.GetTotalCharge() );
                      if(rEV.FitResultExists("partialWaterFitter")){
                          RAT::DS::FitVertex rvertex = rEV.GetFitResult("partialWaterFitter").GetVertex(0);
                          if( rvertex.ContainsPosition() && rvertex.ValidPosition() ) {
                              hposxy->Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                          }
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
      c1->SaveAs(("nhits_"+outname+postfix+".png").c_str());
      c1->SaveAs(("nhits_"+outname+postfix+".pdf").c_str());
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
      c2->SaveAs(("totalQ_"+outname+postfix+".png").c_str());
      c2->SaveAs(("totalQ_"+outname+postfix+".pdf").c_str());
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
      
      TCanvas* c3 = new TCanvas("c3","c3",600,500);
      //gStyle->SetPalette(51,0);
      //gStyle->SetNumberContours(999);
      hposxy->GetYaxis()->SetTitle( "Y position" );
      hposxy->GetXaxis()->SetTitle( "X position" );
      hposxy->SetMarkerStyle(20);
      hposxy->Draw("colzsame");
      title_latex->DrawLatex(0.55, 0.94, ("Run: "+run_info.first + " Subrun: " + run_info.second).c_str()  );
      c3->SetRightMargin(0.15);
      c3->Update();
      c3->SaveAs(("posxy_"+outname+postfix+".png").c_str());
      c3->SaveAs(("posxy_"+outname+postfix+".pdf").c_str());
      c3->Clear();
      
      delete hNHits;
      delete hTotalQ;
      delete hposxy;
      delete c1;
      delete c2;
      delete c3;
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
  c10->SaveAs(("nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c10->SaveAs(("nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
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
  c11->SaveAs(("totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c11->SaveAs(("totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  delete c10;
  delete c11;
  delete hNHits_vs_run;
  delete hTotalQ_vs_run;
 
}

void makeRunPlots(std::string filelist)
{
  //gROOT->SetBatch(kTRUE);
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
  CreateRunPlots(files,ntuple);
  return;
}
