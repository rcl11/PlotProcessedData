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
#include <fstream>
#include <bitset>
#include <map>

#include <boost/version.hpp>

#include "../interface/TH1DPlot.hh"


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

/*map<std::string,1DPlot> ParseConfigs(std::string config_dir="configs/"){




}*/

void CreateRunPlots( const std::vector<std::string>& files, bool ntuple=true, std::string output_dir="plots/" )
{
    
  std::string postfix = "";
  if(ntuple) postfix = "_ntuple";
  for(unsigned i=0; i<files.size(); ++i){
      
      TCanvas* c1 = new TCanvas();
      
      //=============================================================================
      //Create a collection of histograms
      TH1D* hNHits = new TH1D( "hNHits", "Number of hits per event", 100, 0.0, 100.0 );
      
      TFile *f = new TFile(files[i].c_str());
      
      int n_events=0;
      
      //=============================================================================
      //Fill histograms using info from ratds or ntuple
      if(ntuple) {
          TTree *t1 = (TTree*)f->Get("output");
          Int_t nhits;
          ULong64_t flag;
          ULong64_t applied_flag;
          t1->SetBranchAddress("nhits",&nhits);
          t1->SetBranchAddress("dcFlagged",&flag);
          t1->SetBranchAddress("dcApplied",&applied_flag);

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
                n_events++;
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
      
      
  }
 
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

int main(){

makeRunPlots("filelists/filelist_rat625_waterFit_ntuple.dat","dump/");

return 0;

}
