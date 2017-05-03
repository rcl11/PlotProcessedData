#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/DataCleaningUtility.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Meta.hh>
#include <RAT/BitManip.hh>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TColor.h>
#include <TFile.h>
#include <TTree.h>
#include <TCut.h>
#include <string>
#include <bitset>
#include <map>

#include <boost/version.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

#include "../interface/THPlot.hh"
#include "../interface/styles.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;


//Function to pull run info from the filename
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

//Create map to store all plots requested in config files and set style choices
std::map<std::string,THPlot::THPlot> ParseConfigs(std::string config_dir="config/"){
  std::map<std::string,THPlot::THPlot> plot_map;

  fs::path p(config_dir.c_str());
  fs::directory_iterator end_itr;

  for (fs::directory_iterator itr(p); itr != end_itr; ++itr)
  {
    if (fs::is_regular_file(itr->path())) {
        std::string current_file = itr->path().string();
        //Create plot object setting styles from the config file
        THPlot::THPlot plot(current_file.c_str());
        current_file = current_file.erase(current_file.find(".cfg"),current_file.size()-1);
        current_file = current_file.erase(0,current_file.find("/")+1);
        plot_map[current_file] = plot;
     }
  }
  return plot_map;
}

void CreateRunPlots( const std::vector<std::string>& files, bool ntuple=true, std::string output_dir="plots/" )
{
    
  TH1D* hnhits_vs_run = new TH1D( "hnhits_vs_run", "Mean number of hits per run", files.size(), 0.0, files.size() );
  TH1D* htotalQ_vs_run = new TH1D( "htotalQ_vs_run", "Mean total charge per run", files.size(), 0.0, files.size() );
  TH1D* hNEvents_vs_run = new TH1D( "hNEvents_vs_run", "Number of events per run", files.size(), 0.0, files.size() );
  TH1D* hNCleanEvents_vs_run = new TH1D( "hNCleanEvents_vs_run", "Number of clean events per run", files.size(), 0.0, files.size() );
  TH1D* hNCleanEventsnorm_vs_run = new TH1D( "hNCleanEventsnorm_vs_run", "Normalised number of clean events per run", files.size(), 0.0, files.size() );
  TH1D* hFracCleanEvents_vs_run = new TH1D( "hFracCleanEvents_vs_run", "Fraction of clean events per run", files.size(), 0.0, files.size() );
  std::string start_run;
  std::string end_run;
  std::string postfix = "";
  if(ntuple) postfix = "_ntuple";
  double run_duration = 0;
  
  for(unsigned i=0; i<files.size(); ++i){
      
    //Create a collection of histograms
    std::map<std::string,THPlot::THPlot> plot_map = ParseConfigs();
    TH1D hnhits = plot_map["nhits"].GetHist(); 
    TH1D hnhits_pulseGT = plot_map["nhits_pulseGT"].GetHist(); 
    TH1D hnhits_N100M = plot_map["nhits_N100M"].GetHist(); 
    TH1D hnhits_N100H = plot_map["nhits_N100H"].GetHist(); 
    TH1D hnhits_N20 = plot_map["nhits_N20"].GetHist(); 
    TH1D hnhits_ESUMH = plot_map["nhits_ESUMH"].GetHist(); 
    TH1D hnhits_OWLEH = plot_map["nhits_OWLEH"].GetHist(); 
    TH1D htotalQ = plot_map["totalQ"].GetHist(); 
    TH1D htotalQ_pulseGT = plot_map["totalQ_pulseGT"].GetHist(); 
    TH1D htotalQ_N100M = plot_map["totalQ_N100M"].GetHist(); 
    TH1D htotalQ_N100H = plot_map["totalQ_N100H"].GetHist(); 
    TH1D htotalQ_N20 = plot_map["totalQ_N20"].GetHist(); 
    TH1D htotalQ_ESUMH = plot_map["totalQ_ESUMH"].GetHist(); 
    TH1D htotalQ_OWLEH = plot_map["totalQ_OWLEH"].GetHist(); 
    TH1D hfitValid = plot_map["fitValid"].GetHist(); 
    TH1D hitr = plot_map["itr"].GetHist(); 
    TH1D hposx = plot_map["posx"].GetHist(); 
    TH1D hposy = plot_map["posy"].GetHist(); 
    TH1D hposz = plot_map["posz"].GetHist(); 
    TH1D hposR = plot_map["posR"].GetHist(); 
    TH1D hposR3 = plot_map["posR3"].GetHist(); 
    TH1D hrpmt = plot_map["rpmt"].GetHist(); 
    TH1D hxpmt = plot_map["xpmt"].GetHist(); 
    TH1D hypmt = plot_map["ypmt"].GetHist(); 
    TH1D hzpmt = plot_map["zpmt"].GetHist(); 
    TH1D htpmt = plot_map["tpmt"].GetHist(); 
    TH2D hposxy = plot_map["posxy"].Get2DHist(); 
    TH2D hposrz = plot_map["posrz"].Get2DHist(); 
    TH2D hposRz = plot_map["posRz"].Get2DHist(); 
    TH2D hnhitsz = plot_map["nhitsz"].Get2DHist(); 
    TH1D hduration = plot_map["duration"].GetHist(); 
    TH1D htrigger = plot_map["trigger"].GetHist(); 
  
    TFile *f = new TFile(files[i].c_str());
    
    std::vector<std::string> trig_names = {"N100L","N100M","N100H","N20","N20LB","ESUML","ESUMH","OWLN","OWLEL","OWLEH","PULGT"}; 
    for(unsigned h=0;h<trig_names.size();h++){
      htrigger.GetXaxis()->SetBinLabel(h+1,trig_names[h].c_str());
    }
  
    int n_cleanevents=0;
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
      Int_t uTDays, uTSecs, uTNSecs; 
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
      t1->SetBranchAddress("uTDays",&uTDays);
      t1->SetBranchAddress("uTSecs",&uTSecs);
      t1->SetBranchAddress("uTNSecs",&uTNSecs);


      Long64_t nentries = t1->GetEntries();
      int start_days = 0;
      int start_secs = 0;
      int start_nsecs = 0;
      int end_days = 0;
      int end_secs = 0;
      int end_nsecs = 0;
      
      for (Long64_t j=0;j<nentries;j++) {
        t1->GetEntry(j);
        if(j==0) {
          start_days = uTDays;
          start_secs = uTSecs;
          start_nsecs = uTNSecs;
        }
        if(j==nentries-1) {
          end_days = uTDays;
          end_secs = uTSecs;
          end_nsecs = uTNSecs;
        }
        //analysis_mask
        bool dataclean = ( (flag & 0b111111111111110) == 0b111111111111110);
        //analysis_mask
        bool compatibility_cut = (applied_flag & 0b111111111111110) == 0b111111111111110;
        if(dataclean && compatibility_cut) {
          hnhits.Fill(nhits);
          htotalQ.Fill(charge);
          //Fill some nhits and total q plots for different triggers fired
          std::bitset<32> bits = std::bitset<32>(triggerWord);
          //Numbers to test taken from RAT documentation
          if(bits.test(10)){
            hnhits_pulseGT.Fill(nhits);
            htotalQ_pulseGT.Fill(charge);
            htrigger.Fill("PULGT", 1);
          }
          if(bits.test(0)){
            htrigger.Fill("N100L", 1);
          }
          if(bits.test(1)){
            hnhits_N100M.Fill(nhits);
            htotalQ_N100M.Fill(charge);
            htrigger.Fill("N100M", 1);
          }
          if(bits.test(2)){
            hnhits_N100H.Fill(nhits);
            htotalQ_N100H.Fill(charge);
            htrigger.Fill("N100H", 1);
          }
          if(bits.test(3)){
            hnhits_N20.Fill(nhits);
            htotalQ_N20.Fill(charge);
            htrigger.Fill("N20", 1);
          }
          if(bits.test(4)){
            htrigger.Fill("N20LB", 1);
          }
          if(bits.test(5)){
            htrigger.Fill("ESUML", 1);
          }
          if(bits.test(6)){
            hnhits_ESUMH.Fill(nhits);
            htotalQ_ESUMH.Fill(charge);
            htrigger.Fill("ESUMH", 1);
          }
          if(bits.test(7)){
            htrigger.Fill("OWLN", 1);
          }
          if(bits.test(8)){
            htrigger.Fill("OWLEL", 1);
          }
          if(bits.test(9)){
            hnhits_OWLEH.Fill(nhits);
            htotalQ_OWLEH.Fill(charge);
            htrigger.Fill("OWLEH", 1);
          }
          //std::cout << "testing pulseGT trigger " << bits.test(10) << " " << htrigger->GetBinContent(11) << std::endl; 
          hfitValid.Fill(fit_valid); 
          if(fit_valid){
            hposrz.Fill(sqrt(posx*posx + posy*posy), posz);
            hposx.Fill(posx);
            hposy.Fill(posy);
            hposz.Fill(posz);
            hposxy.Fill(posx,posy);
            hnhitsz.Fill(nhits,posz);
            double R = sqrt(posx*posx + posy*posy + posz*posz);
            hposR.Fill(R);
            hposR3.Fill(pow(R,3)/pow(6005.3,3));
            hposRz.Fill(R, posz);
            hitr.Fill(itr); 
          }
          n_cleanevents++;
        }
        n_events++;
      }
      run_duration = (end_days-start_days)*60*60*24 + (end_secs-start_secs) + ((end_nsecs-start_nsecs) * 1E-9);
      hduration.Fill(run_duration);
    } else {
      RAT::DU::DSReader dsReader( files[i] );
      const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
      //Add data cleaning cut.
      ULong64_t rDataCleaningWord = RAT::GetDataCleaningWord( "analysis_mask" );
      RAT::DS::UniversalTime start_time;
      RAT::DS::UniversalTime end_time;
      for( size_t iEntry = 0; iEntry < dsReader.GetEntryCount(); iEntry++ )
      {
        const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry );
        for( size_t iEv = 0; iEv < rDS.GetEVCount(); iEv++ )
        {   
          RAT::DS::EV rEV = rDS.GetEV(iEv);
          if(iEntry == 0 && iEv == 0 ) start_time = rEV.GetUniversalTime();
          if(iEntry == dsReader.GetEntryCount()-1 && iEv == rDS.GetEVCount()-1 ) end_time = rEV.GetUniversalTime();
          const RAT::DS::Meta& rMeta = dsReader.GetMeta();
          n_events++;
          if(RAT::EventIsClean( rEV, rMeta, rDataCleaningWord )){
            n_cleanevents++;
            hnhits.Fill( rEV.GetNhits() );
            htotalQ.Fill( rEV.GetTotalCharge() );
            //Fill some nhits and total q plots for different triggers fired
            //std::cout << std::bitset<32>(rEV.GetTrigType())/*.to_string()*/ << std::endl;
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N100Low)){
                htrigger.Fill("N100L", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N100Med)){
                hnhits_N100M.Fill(rEV.GetNhits()); 
                htotalQ_N100M.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("N100M", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N100High)){
                hnhits_N100H.Fill(rEV.GetNhits()); 
                htotalQ_N100H.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("N100H", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N20)){
                hnhits_N20.Fill(rEV.GetNhits()); 
                htotalQ_N20.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("N20", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::N20LB)){
                htrigger.Fill("N20LB", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::ESLow)){
                htrigger.Fill("ESUML", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::ESHigh)){
                hnhits_ESUMH.Fill(rEV.GetNhits()); 
                htotalQ_ESUMH.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("ESUMH", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::OWLN)){
                htrigger.Fill("OWLN", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::OWLESLow)){
                htrigger.Fill("OWLEL", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::OWLESHigh)){
                hnhits_OWLEH.Fill(rEV.GetNhits()); 
                htotalQ_OWLEH.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("OWLEH", 1);
            }
            if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::PulseGT)){
                hnhits_pulseGT.Fill(rEV.GetNhits()); 
                htotalQ_pulseGT.Fill( rEV.GetTotalCharge() );
                htrigger.Fill("PULGT", 1);
            }
                  
            if(rEV.FitResultExists("waterFitter") && rEV.GetFitResult("waterFitter").GetValid()){
                RAT::DS::FitVertex rvertex = rEV.GetFitResult("waterFitter").GetVertex(0);
              if( rvertex.ContainsPosition() && rvertex.ValidPosition() ) {
                hfitValid.Fill(1.);
                double R = sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y() + rvertex.GetPosition().Z()*rvertex.GetPosition().Z());
                hposx.Fill(rvertex.GetPosition().X());
                hposy.Fill(rvertex.GetPosition().Y());
                hposz.Fill(rvertex.GetPosition().Z());
                hposxy.Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                hnhitsz.Fill(rEV.GetNhits(),rvertex.GetPosition().Z());
                hposRz.Fill(R,rvertex.GetPosition().Z());
                hposrz.Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X() + rvertex.GetPosition().Y()*rvertex.GetPosition().Y()),rvertex.GetPosition().Z());
                hposR.Fill(R );
                hposR3.Fill(pow(R,3)/pow(6005.3,3));
                hitr.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"));
              } else hfitValid.Fill(0.);
            } else hfitValid.Fill(0.);
            RAT::DS::CalPMTs& calpmts = rEV.GetCalPMTs();
            //NOTE: i think we should return to GetCount() and GetPMT() after next reprocess (assuming we want to include the "inward" PMTs (normal + HQE) not just "normal" here)
            for(unsigned int ipmt=0;ipmt<calpmts.GetNormalCount();ipmt++){
              TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetNormalPMT(ipmt).GetID());
              double pmt_r = pmtpos.Mag();
              hrpmt.Fill(pmt_r);
              hxpmt.Fill(pmtpos.X());
              hypmt.Fill(pmtpos.Y());
              hzpmt.Fill(pmtpos.Z());
              htpmt.Fill((calpmts.GetNormalPMT(ipmt)).GetTime());
            }
          }
        }
      }
      run_duration = ((end_time-start_time).GetDays())*60*60*24 + ((end_time-start_time).GetSeconds()) + ((end_time-start_time).GetNanoSeconds() * 1E-9);
      hduration.Fill(run_duration);
    }
    plot_map["nhits"].SetHist(hnhits);
    plot_map["nhits_pulseGT"].SetHist(hnhits_pulseGT);
    plot_map["nhits_N100M"].SetHist(hnhits_N100M);
    plot_map["nhits_N100H"].SetHist(hnhits_N100H);
    plot_map["nhits_N20"].SetHist(hnhits_N20);
    plot_map["nhits_ESUMH"].SetHist(hnhits_ESUMH);
    plot_map["nhits_OWLEH"].SetHist(hnhits_OWLEH);
    plot_map["totalQ"].SetHist(htotalQ);
    plot_map["totalQ_pulseGT"].SetHist(htotalQ_pulseGT);
    plot_map["totalQ_N100M"].SetHist(htotalQ_N100M);
    plot_map["totalQ_N100H"].SetHist(htotalQ_N100H);
    plot_map["totalQ_N20"].SetHist(htotalQ_N20);
    plot_map["totalQ_ESUMH"].SetHist(htotalQ_ESUMH);
    plot_map["totalQ_OWLEH"].SetHist(htotalQ_OWLEH);
    plot_map["fitValid"].SetHist(hfitValid);
    plot_map["itr"].SetHist(hitr);
    plot_map["posx"].SetHist(hposx);
    plot_map["posy"].SetHist(hposy);
    plot_map["posz"].SetHist(hposz);
    plot_map["posR"].SetHist(hposR);
    plot_map["posR3"].SetHist(hposR3);
    plot_map["rpmt"].SetHist(hrpmt);
    plot_map["xpmt"].SetHist(hxpmt);
    plot_map["ypmt"].SetHist(hypmt);
    plot_map["zpmt"].SetHist(hzpmt);
    plot_map["tpmt"].SetHist(htpmt);
    plot_map["duration"].SetHist(hduration);
    plot_map["trigger"].SetHist(htrigger);
    plot_map["posxy"].Set2DHist(hposxy);
    plot_map["nhitsz"].Set2DHist(hnhitsz);
    plot_map["posrz"].Set2DHist(hposrz);
    plot_map["posRz"].Set2DHist(hposRz);
  
    //=============================================================================
    //Make some plots
   
    std::pair<std::string,std::string> run_info = ParseRunInfo(files[i]);
    std::string outname = run_info.first + "_" + run_info.second; 
    std::map<std::string, THPlot::THPlot>::iterator it;
    for ( it = plot_map.begin(); it != plot_map.end(); it++ )
    {
        THPlot::THPlot plot = it->second;
        plot.SetOutFilename(output_dir + plot.GetOutFilename() + "_" + outname + postfix);
        plot.SetRunInfo(run_info);
        plot.GeneratePlot();
    }

    //Code for the vs run plots. For now dont use the class for this
    hnhits_vs_run->Fill(i, hnhits.GetMean());
    hnhits_vs_run->SetBinError(i+1,hnhits.GetMeanError());
    std::string bin_label = run_info.first;
    if(i==0) start_run = bin_label;
    if(i==files.size()-1) end_run = bin_label;
    hnhits_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    htotalQ_vs_run->Fill(i, htotalQ.GetMean());
    htotalQ_vs_run->SetBinError(i+1,htotalQ.GetMeanError());
    htotalQ_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    hNCleanEvents_vs_run->Fill(i, n_cleanevents);
    hNCleanEvents_vs_run->SetBinError(i+1,sqrt(n_cleanevents));
    hNCleanEvents_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    hNEvents_vs_run->Fill(i, n_events);
    hNEvents_vs_run->SetBinError(i+1,sqrt(n_events));
    hNEvents_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    hNCleanEventsnorm_vs_run->Fill(i, n_cleanevents/run_duration);
    hNCleanEventsnorm_vs_run->SetBinError(i+1,sqrt(n_cleanevents/run_duration));
    hNCleanEventsnorm_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    hFracCleanEvents_vs_run->Fill(i, float(n_cleanevents)/float(n_events));
    hFracCleanEvents_vs_run->SetBinError(i+1,sqrt( pow((sqrt(n_cleanevents)/n_cleanevents),2) + pow((sqrt(n_events)/n_events),2)));
    hFracCleanEvents_vs_run->GetXaxis()->SetBinLabel(i+1,bin_label.c_str());
    delete f;
  }
  

  SetStyle();
  TCanvas* c100 = new TCanvas("c100","c100",1000,400);
  hnhits_vs_run->SetMarkerColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_run->SetMarkerStyle(20);
  hnhits_vs_run->SetLineColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_run->GetYaxis()->SetTitle( "Mean #Hits" );
  hnhits_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hnhits_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hnhits_vs_run->Draw("PE1");
  c100->SetLeftMargin(0.1);
  c100->SetBottomMargin(0.18);
  c100->Update();
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  htotalQ_vs_run->SetMarkerColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_run->SetMarkerStyle(20);
  htotalQ_vs_run->SetLineColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_run->GetYaxis()->SetTitle( "Mean totalQ" );
  htotalQ_vs_run->GetXaxis()->SetTitle( "Run ID" );
  htotalQ_vs_run->GetXaxis()->SetTitleOffset(1.9);
  htotalQ_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNCleanEvents_vs_run->SetMarkerColor(kGreen);
  hNCleanEvents_vs_run->SetMarkerStyle(20);
  hNCleanEvents_vs_run->SetLineColor(kGreen);
  hNCleanEvents_vs_run->GetYaxis()->SetTitle( "# Clean Events" );
  hNCleanEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNCleanEvents_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNCleanEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"ncleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"ncleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNEvents_vs_run->SetMarkerColor(kGreen);
  hNEvents_vs_run->SetMarkerStyle(20);
  hNEvents_vs_run->SetLineColor(kGreen);
  hNEvents_vs_run->GetYaxis()->SetTitle( "# Events" );
  hNEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNEvents_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNCleanEventsnorm_vs_run->SetMarkerColor(kGreen+3);
  hNCleanEventsnorm_vs_run->SetMarkerStyle(20);
  hNCleanEventsnorm_vs_run->SetLineColor(kGreen+3);
  hNCleanEventsnorm_vs_run->GetYaxis()->SetTitle( "# Clean events per second" );
  hNCleanEventsnorm_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNCleanEventsnorm_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hNCleanEventsnorm_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"ncleaneventsnorm_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"ncleaneventsnorm_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hFracCleanEvents_vs_run->SetMarkerColor(kMagenta);
  hFracCleanEvents_vs_run->SetMarkerStyle(20);
  hFracCleanEvents_vs_run->SetLineColor(kMagenta);
  hFracCleanEvents_vs_run->GetYaxis()->SetTitle( "Fraction of clean Events" );
  hFracCleanEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hFracCleanEvents_vs_run->GetXaxis()->SetTitleOffset(1.9);
  hFracCleanEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"fraccleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"fraccleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  delete c100;
  delete hnhits_vs_run;
  delete htotalQ_vs_run;
  delete hNEvents_vs_run;
  delete hNCleanEvents_vs_run;
  delete hFracCleanEvents_vs_run;
  delete hNCleanEventsnorm_vs_run;
 
}


int main(int argc, char** argv){

   po::options_description desc("Options"); 
   std::string filelist;
   std::string directory = "plots/";
   desc.add_options() 
    ("filelist", po::value<std::string>(&filelist)) 
    ("directory", po::value<std::string>(&directory)->default_value("plots/"));

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc),  
                  vm);
  po::notify(vm);

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
  CreateRunPlots(files,ntuple,directory);


  return 0;
}
