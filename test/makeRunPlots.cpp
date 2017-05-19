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
#include <vector>
#include <bitset>
#include <map>
#include <stdlib.h>
#include <algorithm>

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
std::map<std::string,THPlot::THPlot> ParseConfigs(std::string config_dir="config/", std::string totals=""){
  std::map<std::string,THPlot::THPlot> plot_map;

  fs::path p(config_dir.c_str());
  fs::directory_iterator end_itr;

  for (fs::directory_iterator itr(p); itr != end_itr; ++itr)
  {
    if (fs::is_regular_file(itr->path())) {
        std::string current_file = itr->path().string();
        //Create plot object setting styles from the config file
        THPlot::THPlot plot(current_file.c_str(), totals);
        current_file = current_file.erase(current_file.find(".cfg"),current_file.size()-1);
        current_file = current_file.erase(0,current_file.find("/")+1);
        plot_map[current_file] = plot;
     }
  }
  return plot_map;
}

//Finds cavity temperature for a given universal time. Based on information from deltaV, input by hand from the output of Teal's script
TGraph* MakeCavityTempHist(){
  fstream tempfile("data/temps.dat");
  std::string reftime_str, temp_str;
  std::getline(tempfile,reftime_str);
  std::getline(tempfile,temp_str);
  std::vector<double> reftimes, temps;
  std::stringstream reftime_ss(reftime_str);
  std::stringstream temp_ss(temp_str);
  std::string temp, reftime;
  while(std::getline(reftime_ss,reftime,',')) {
    double corrected_time =  atof(reftime.c_str())/(60*60) - (40*24*365);
    reftimes.push_back(corrected_time);
  }
  while(std::getline(temp_ss,temp,',')) {
    temps.push_back(atof(temp.c_str()));
  }
  std::reverse(reftimes.begin(),reftimes.end());
  std::reverse(temps.begin(),temps.end());
  //Remove any points where temp = 0, this happens when the database returns N/A for whatever reason 
  std::vector<double> corrreftimes, corrtemps; 
  for(unsigned k = 0; k<reftimes.size(); k++){
    if(!(temps[k]==0)){
      corrreftimes.push_back(reftimes[k]);
      corrtemps.push_back(temps[k]);
    }
  }
  SetStyle();
  TGraph* timetemp = new TGraph(corrreftimes.size(), &corrreftimes[0], &corrtemps[0] );
  TCanvas* ctest = new TCanvas("ctest","ctest",1200,500);
  timetemp->GetXaxis()->SetTitle("Universal time (hours)");
  timetemp->GetYaxis()->SetTitle("Temperature (C)");
  timetemp->Draw();
  ctest->SaveAs("data/temp_dist.png");
  ctest->SaveAs("data/temp_dist.pdf");
  delete ctest;
  return timetemp;
}

double GetCavityTemp(double time, TGraph* temphist){
    return temphist->Eval(time); 
}

void CreateRunPlots( const std::vector<std::vector<std::string> >& files, bool ntuple=true, std::string output_dir="plots/" )
{
  SetStyle();
  TH1D* hnhits_vs_run = new TH1D( "hnhits_vs_run", "Mean number of hits per run", files.size(), 0.0, files.size() );
  hnhits_vs_run->Sumw2();
  TH1D* hnhits_vs_time = new TH1D( "hnhits_vs_time", "Number of hits per 5 minutes", 3600, 64300, 64600 );
  hnhits_vs_time->Sumw2();
  TH1D* htotalQ_vs_run = new TH1D( "htotalQ_vs_run", "Mean total charge per run", files.size(), 0.0, files.size() );
  htotalQ_vs_run->Sumw2();
  TH1D* htotalQ_vs_time = new TH1D( "htotalQ_vs_time", "Total charge per 5 minutes", 3600, 64300, 64600 );
  htotalQ_vs_time->Sumw2();
  TH1D* hNEvents_vs_run = new TH1D( "hNEvents_vs_run", "Number of events per run", files.size(), 0.0, files.size() );
  hNEvents_vs_run->Sumw2();
  TH1D* hNEvents_vs_time = new TH1D( "hNEvents_vs_rtime", "Number of events per 5 minutes",3600,64300,64600 );
  hNEvents_vs_time->Sumw2();
  TH1D* hNCleanEvents_vs_run = new TH1D( "hNCleanEvents_vs_run", "Number of clean events per run", files.size(), 0.0, files.size() );
  hNCleanEvents_vs_run->Sumw2();
  TH1D* hNCleanEvents_vs_time = new TH1D( "hNCleanEvents_vs_time", "Number of clean events per 5 minutes",3600,64300,64600   );
  hNCleanEvents_vs_time->Sumw2();
  TH1D* hNCleanEventsnorm_vs_run = new TH1D( "hNCleanEventsnorm_vs_run", "Normalised number of clean events per run", files.size(), 0.0, files.size() );
  hNCleanEventsnorm_vs_run->Sumw2();
  TH1D* hFracCleanEvents_vs_run = new TH1D( "hFracCleanEvents_vs_run", "Fraction of clean events per run", files.size(), 0.0, files.size() );
  hFracCleanEvents_vs_run->Sumw2();
  TH1D* hFracCleanEvents_vs_time = new TH1D( "hFracCleanEvents_vs_time", "Fraction of clean events per 5 minutes",3600,64300,64600  );
  hFracCleanEvents_vs_time->Sumw2();
  TH1D* hTemp_vs_run = new TH1D( "hTemp_vs_run", "Mean cavity temperature per run", files.size(), 0.0, files.size() );
  hTemp_vs_run->Sumw2();

  std::string start_run = "";
  std::string end_run = "";
  double start_run_time=0;
  double end_run_time=0;
  std::string postfix = "";
  if(ntuple) postfix = "_ntuple";
  double run_duration = 0;
  int file_count=0;
  TGraph * temphist = MakeCavityTempHist();
  std::map<std::string,THPlot::THPlot> plot_map_totals = ParseConfigs("config/","total");
  std::vector<std::map<std::string,THPlot::THPlot> > plot_maps;
  std::vector<std::string> trig_names = {"N100L","N100M","N100H","N20","N20LB","ESUML","ESUMH","OWLN","OWLEL","OWLEH","PULGT","Prescale","Pedestal"}; 
  std::vector<std::string> dataclean_names = {"prescale","zerozerocut","crateisotropy","ftscut","flashergeocut","icttimespread","junkcut","muontag","neckcut","owlcut","qcluster","qvnhit","qvt","ringoffire","tpmuonshort"}; 
  //This should really go inside the class but im feeling hacky today
  TH1D htrigger_total = plot_map_totals["trigger"].GetHist();
  TH1D hdataclean_total = plot_map_totals["dataclean"].GetHist();
  for(unsigned h=0;h<trig_names.size();h++){
    htrigger_total.GetXaxis()->SetBinLabel(h+1,trig_names[h].c_str());
  }
  for(unsigned h=0;h<dataclean_names.size();h++){
    hdataclean_total.GetXaxis()->SetBinLabel(h+1,dataclean_names[h].c_str());
  }
  hdataclean_total.GetXaxis()->SetBinLabel(dataclean_names.size()+1,"allevents");
  plot_map_totals["trigger"].SetHist(htrigger_total);
  plot_map_totals["dataclean"].SetHist(hdataclean_total);

  std::vector<std::string> allruns;
  for(unsigned i=0; i<files.size(); ++i){
    std::cout << "Creating plots for run number " << i << std::endl;
      
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
    //Fit related
    TH1D hfitValid = plot_map["fitValid"].GetHist(); 
    TH1D hitr = plot_map["itr"].GetHist(); 
    TH1D htime = plot_map["time"].GetHist(); 
    TH1D hposx = plot_map["posx"].GetHist(); 
    TH1D hposy = plot_map["posy"].GetHist(); 
    TH1D hposz = plot_map["posz"].GetHist(); 
    TH1D hposR = plot_map["posR"].GetHist(); 
    TH1D hposR3 = plot_map["posR3"].GetHist(); 
    TH2D hposxy = plot_map["posxy"].Get2DHist(); 
    TH2D htimeposx = plot_map["timeposx"].Get2DHist(); 
    TH2D htimeposy = plot_map["timeposy"].Get2DHist(); 
    TH2D htimeposz = plot_map["timeposz"].Get2DHist(); 
    TH2D herrposx = plot_map["errposx"].Get2DHist(); 
    TH2D herrposy = plot_map["errposy"].Get2DHist(); 
    TH2D herrposz = plot_map["errposz"].Get2DHist(); 
    TH2D herrposxnhits = plot_map["errposxnhits"].Get2DHist(); 
    TH2D herrposxitr = plot_map["errposxitr"].Get2DHist(); 
    TH2D herrposynhits = plot_map["errposynhits"].Get2DHist(); 
    TH2D herrposyitr = plot_map["errposyitr"].Get2DHist(); 
    TH2D herrposznhits = plot_map["errposznhits"].Get2DHist(); 
    TH2D herrposzitr = plot_map["errposzitr"].Get2DHist(); 
    TH2D herrtimex = plot_map["errtimex"].Get2DHist(); 
    TH2D herrtimey = plot_map["errtimey"].Get2DHist(); 
    TH2D herrtimez = plot_map["errtimez"].Get2DHist(); 
    TH2D herrenergy = plot_map["errenergy"].Get2DHist(); 
    TH2D hposrhoz = plot_map["posrhoz"].Get2DHist(); 
    TH2D hposRz = plot_map["posRz"].Get2DHist();
    TH1D hrpmt = plot_map["rpmt"].GetHist(); 
    TH1D hxpmt = plot_map["xpmt"].GetHist(); 
    TH1D hypmt = plot_map["ypmt"].GetHist(); 
    TH1D hzpmt = plot_map["zpmt"].GetHist(); 
    TH1D htpmt = plot_map["tpmt"].GetHist(); 
    TH1D hbeta14 = plot_map["beta14"].GetHist(); 
    TH1D henergy = plot_map["energy"].GetHist(); 
    TH2D hnhitsz = plot_map["nhitsz"].Get2DHist(); 
    //Second copy of fit related for some good fit conditions
    TH1D hitrgoodfit = plot_map["itrgoodfit"].GetHist(); 
    TH1D htimegoodfit = plot_map["timegoodfit"].GetHist(); 
    TH1D hposxgoodfit = plot_map["posxgoodfit"].GetHist(); 
    TH1D hposygoodfit = plot_map["posygoodfit"].GetHist(); 
    TH1D hposzgoodfit = plot_map["poszgoodfit"].GetHist(); 
    TH1D hposRgoodfit = plot_map["posRgoodfit"].GetHist(); 
    TH1D hposR3goodfit = plot_map["posR3goodfit"].GetHist(); 
    TH2D hposxygoodfit = plot_map["posxygoodfit"].Get2DHist(); 
    TH2D htimeposxgoodfit = plot_map["timeposxgoodfit"].Get2DHist(); 
    TH2D htimeposygoodfit = plot_map["timeposygoodfit"].Get2DHist(); 
    TH2D htimeposzgoodfit = plot_map["timeposzgoodfit"].Get2DHist(); 
    TH2D herrposxgoodfit = plot_map["errposxgoodfit"].Get2DHist(); 
    TH2D herrposygoodfit = plot_map["errposygoodfit"].Get2DHist(); 
    TH2D herrposzgoodfit = plot_map["errposzgoodfit"].Get2DHist(); 
    TH2D herrposxnhitsgoodfit = plot_map["errposxnhitsgoodfit"].Get2DHist(); 
    TH2D herrposxitrgoodfit = plot_map["errposxitrgoodfit"].Get2DHist(); 
    TH2D herrposynhitsgoodfit = plot_map["errposynhitsgoodfit"].Get2DHist(); 
    TH2D herrposyitrgoodfit = plot_map["errposyitrgoodfit"].Get2DHist(); 
    TH2D herrposznhitsgoodfit = plot_map["errposznhitsgoodfit"].Get2DHist(); 
    TH2D herrposzitrgoodfit = plot_map["errposzitrgoodfit"].Get2DHist(); 
    TH2D herrtimexgoodfit = plot_map["errtimexgoodfit"].Get2DHist(); 
    TH2D herrtimeygoodfit = plot_map["errtimeygoodfit"].Get2DHist(); 
    TH2D herrtimezgoodfit = plot_map["errtimezgoodfit"].Get2DHist(); 
    TH2D herrenergygoodfit = plot_map["errenergygoodfit"].Get2DHist(); 
    TH2D hposrhozgoodfit = plot_map["posrhozgoodfit"].Get2DHist(); 
    TH2D hposRzgoodfit = plot_map["posRzgoodfit"].Get2DHist();
    TH1D hrpmtgoodfit = plot_map["rpmtgoodfit"].GetHist(); 
    TH1D hxpmtgoodfit = plot_map["xpmtgoodfit"].GetHist(); 
    TH1D hypmtgoodfit = plot_map["ypmtgoodfit"].GetHist(); 
    TH1D hzpmtgoodfit = plot_map["zpmtgoodfit"].GetHist(); 
    TH1D htpmtgoodfit = plot_map["tpmtgoodfit"].GetHist(); 
    TH1D hbeta14goodfit = plot_map["beta14goodfit"].GetHist(); 
    TH1D henergygoodfit = plot_map["energygoodfit"].GetHist(); 
    TH2D hnhitszgoodfit = plot_map["nhitszgoodfit"].Get2DHist(); 
    
    TH2D hnhitstemp = plot_map["nhitstemp"].Get2DHist(); 
    TH1D hduration = plot_map["duration"].GetHist(); 
    TH1D htrigger = plot_map["trigger"].GetHist(); 
    TH1D hdataclean = plot_map["dataclean"].GetHist(); 
    TH1D htemp = plot_map["temp"].GetHist(); 
    
    for(unsigned h=0;h<trig_names.size();h++){
      htrigger.GetXaxis()->SetBinLabel(h+1,trig_names[h].c_str());
    }
    for(unsigned h=0;h<dataclean_names.size();h++){
      hdataclean.GetXaxis()->SetBinLabel(h+1,dataclean_names[h].c_str());
    }
    hdataclean.GetXaxis()->SetBinLabel(dataclean_names.size()+1,"allevents");
  
    int n_cleanevents=0;
    int n_events=0;
    std::vector<std::string> runfiles = files[i];
    int start_days = 0;
    int start_secs = 0;
    int start_nsecs = 0;
    int end_days = 0;
    int end_secs = 0;
    int end_nsecs = 0;
          
    for(unsigned subrun=0; subrun<runfiles.size(); subrun++){ 
        std::cout << "Filling plots for subrun number " << subrun << std::endl;
        TFile *f = new TFile((runfiles[subrun]).c_str());

        //=============================================================================
        //Fill histograms using info from ratds or ntuple
        if(ntuple) {
          TTree *t1 = (TTree*)f->Get("output");
          Int_t nhits;
          Double_t charge;
          ULong64_t flag;
          ULong64_t applied_flag;
          Double_t posx, posy, posz;
          Double_t posxPosError, posyPosError, poszPosError;
          bool fit_valid;
          Double_t itr;
          Int_t triggerWord;
          Double_t beta14;
          Double_t energy;
          Double_t energyPosError;
          Double_t time;
          Double_t timePosError;
          Int_t uTDays, uTSecs, uTNSecs; 
          t1->SetBranchAddress("nhits",&nhits);
          t1->SetBranchAddress("q",&charge);
          t1->SetBranchAddress("dcFlagged",&flag);
          t1->SetBranchAddress("dcApplied",&applied_flag);
          t1->SetBranchAddress("posx",&posx);
          t1->SetBranchAddress("posy",&posy);
          t1->SetBranchAddress("posz",&posz);
          t1->SetBranchAddress("posxPosError",&posxPosError);
          t1->SetBranchAddress("posyPosError",&posyPosError);
          t1->SetBranchAddress("poszPosError",&poszPosError);
          t1->SetBranchAddress("fitValid",&fit_valid);
          t1->SetBranchAddress("triggerWord",&triggerWord);
          t1->SetBranchAddress("itr",&itr);
          t1->SetBranchAddress("uTDays",&uTDays);
          t1->SetBranchAddress("uTSecs",&uTSecs);
          t1->SetBranchAddress("uTNSecs",&uTNSecs);
          t1->SetBranchAddress("beta14",&beta14);
          t1->SetBranchAddress("energy",&energy);
          t1->SetBranchAddress("energyPosError",&energyPosError);
          t1->SetBranchAddress("time",&time);
          t1->SetBranchAddress("timePosError",&timePosError);


          Long64_t nentries = t1->GetEntries();
          for (Long64_t j=0;j<nentries;j++) {
            t1->GetEntry(j);
            if(j==0 && subrun==0) {
              start_days = uTDays;
              start_secs = uTSecs;
              start_nsecs = uTNSecs;
            }
            if(j==nentries-1 && subrun==runfiles.size()-1) {
              end_days = uTDays;
              end_secs = uTSecs;
              end_nsecs = uTNSecs;
            }
            
            double event_time_secs = (uTDays)*60*60*24 + (uTSecs) + ((uTNSecs) * 1E-9);
            if(j==0 && i==0 && subrun==0) start_run_time = event_time_secs;
            if(j==nentries-1 && i==files.size()-1 && subrun==runfiles.size()-1) end_run_time = event_time_secs;
            //analysis_mask
            bool dataclean = ( (flag & 0b111111111111110) == 0b111111111111110);
            //analysis_mask
            bool compatibility_cut = (applied_flag & 0b111111111111110) == 0b111111111111110;
            std::bitset<32> cleanbits = std::bitset<32>(flag);
            for(unsigned g=0; g<dataclean_names.size(); g++){
              if(cleanbits.test(g)) hdataclean.Fill(dataclean_names[g].c_str(),1);
            }
            hdataclean.Fill("allevents",1);
            hNEvents_vs_time->Fill(event_time_secs/(60*60),1);
            if(dataclean && compatibility_cut) {
              hNCleanEvents_vs_time->Fill(event_time_secs/(60*60),1);
              hnhits.Fill(nhits);
              double cavity_temp = GetCavityTemp(event_time_secs/(60*60), temphist);
              hnhitstemp.Fill(nhits,cavity_temp);
              htemp.Fill(cavity_temp);
              hnhits_vs_time->Fill(event_time_secs/(60*60),nhits);
              htotalQ.Fill(charge);
              htotalQ_vs_time->Fill(event_time_secs/(60*60),charge);
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
              if(bits.test(11)){
                htrigger.Fill("Prescale", 1);
              }
              if(bits.test(12)){
                htrigger.Fill("Pedestal", 1);
              }
              hfitValid.Fill(fit_valid); 
              if(fit_valid){
                hposrhoz.Fill(sqrt(posx*posx + posy*posy), posz);
                hposx.Fill(posx);
                hposy.Fill(posy);
                hposz.Fill(posz);
                hposxy.Fill(posx,posy);
                hnhitsz.Fill(nhits,posz);
                hbeta14.Fill(beta14);
                henergy.Fill(energy);
                double R = sqrt(posx*posx + posy*posy + posz*posz);
                hposR.Fill(R);
                hposR3.Fill(pow(R,3)/pow(6005.3,3));
                hposRz.Fill(R, posz);
                hitr.Fill(itr);
                herrposx.Fill(posx,posxPosError);
                herrposy.Fill(posy,posyPosError);
                herrposz.Fill(posz,poszPosError);
                herrposxnhits.Fill(nhits,posxPosError);
                herrposxitr.Fill(itr,posxPosError);
                herrposynhits.Fill(nhits,posyPosError);
                herrposyitr.Fill(itr,posyPosError);
                herrposznhits.Fill(nhits,poszPosError);
                herrposzitr.Fill(itr,poszPosError);
                herrtimex.Fill(posx,timePosError);
                herrtimey.Fill(posy,timePosError);
                herrtimez.Fill(posz,timePosError);
                htime.Fill(time);
                htimeposx.Fill(posx,time);
                htimeposy.Fill(posy,time);
                htimeposz.Fill(posz,time);
                if(itr>0.55){
                  hposrhozgoodfit.Fill(sqrt(posx*posx + posy*posy), posz);
                  hposxgoodfit.Fill(posx);
                  hposygoodfit.Fill(posy);
                  hposzgoodfit.Fill(posz);
                  hposxygoodfit.Fill(posx,posy);
                  hnhitszgoodfit.Fill(nhits,posz);
                  hbeta14goodfit.Fill(beta14);
                  henergygoodfit.Fill(energy);
                  double R = sqrt(posx*posx + posy*posy + posz*posz);
                  hposRgoodfit.Fill(R);
                  hposR3goodfit.Fill(pow(R,3)/pow(6005.3,3));
                  hposRzgoodfit.Fill(R, posz);
                  hitrgoodfit.Fill(itr);
                  herrposxgoodfit.Fill(posx,posxPosError);
                  herrposygoodfit.Fill(posy,posyPosError);
                  herrposzgoodfit.Fill(posz,poszPosError);
                  herrposxnhitsgoodfit.Fill(nhits,posxPosError);
                  herrposxitrgoodfit.Fill(itr,posxPosError);
                  herrposynhitsgoodfit.Fill(nhits,posyPosError);
                  herrposyitrgoodfit.Fill(itr,posyPosError);
                  herrposznhitsgoodfit.Fill(nhits,poszPosError);
                  herrposzitrgoodfit.Fill(itr,poszPosError);
                  herrtimexgoodfit.Fill(posx,timePosError);
                  herrtimeygoodfit.Fill(posy,timePosError);
                  herrtimezgoodfit.Fill(posz,timePosError);
                  htimegoodfit.Fill(time);
                  htimeposxgoodfit.Fill(posx,time);
                  htimeposygoodfit.Fill(posy,time);
                  htimeposzgoodfit.Fill(posz,time);
                }
              }
              n_cleanevents++;
            }
            n_events++;
          }
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
              if(iEntry == 0 && iEv == 0 && subrun==0) start_time = rEV.GetUniversalTime();
              if(iEntry == dsReader.GetEntryCount()-1 && iEv == rDS.GetEVCount()-1 && subrun==runfiles.size()-1 ) end_time = rEV.GetUniversalTime();
              double event_time_secs = ((rEV.GetUniversalTime()).GetDays())*60*60*24 + ((rEV.GetUniversalTime()).GetSeconds()) + ((rEV.GetUniversalTime()).GetNanoSeconds() * 1E-9);
              if(iEntry == 0 && iEv == 0 && i == 0 && subrun==0) start_run_time = event_time_secs;
              if(iEntry == dsReader.GetEntryCount()-1 && iEv == rDS.GetEVCount()-1 && i == files.size()-1 && subrun==runfiles.size()-1) end_run_time = event_time_secs;
              const RAT::DS::Meta& rMeta = dsReader.GetMeta();
              n_events++;
              //Below line may need changing with new PR fixing which pass this works for
              std::bitset<32> cleanbits = std::bitset<32> (rEV.GetDataCleaningFlags().GetFlags(0).GetULong64_t(0));
              for(unsigned g=0; g<dataclean_names.size(); g++){
                if(cleanbits.test(g)) hdataclean.Fill(dataclean_names[g].c_str(),1);
              }
              hdataclean.Fill("allevents",1);
              
              hNEvents_vs_time->Fill(event_time_secs/(60*60),1);
              if(RAT::EventIsClean( rEV, rMeta, rDataCleaningWord )){
              //Will need changing in next RAT release
              //if(RAT::EventIsClean( rEV, rDataCleaningWord ))
                n_cleanevents++;
                hnhits.Fill( rEV.GetNhits() );
                hNCleanEvents_vs_time->Fill(event_time_secs/(60*60),1);
                hnhits_vs_time->Fill(event_time_secs/(60*60),rEV.GetNhits());
                htotalQ.Fill( rEV.GetTotalCharge() );
                htotalQ_vs_time->Fill(event_time_secs/(60*60),rEV.GetTotalCharge());
                double cavity_temp = GetCavityTemp(event_time_secs/(60*60), temphist);
                hnhitstemp.Fill(rEV.GetNhits(),cavity_temp);
                htemp.Fill(cavity_temp);
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
                if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::Prescale)){
                    htrigger.Fill("Prescale", 1);
                }
                if(RAT::BitManip::TestBit(rEV.GetTrigType(), RAT::DU::TrigBits::Pedestal)){
                    htrigger.Fill("Pedestal", 1);
                }
                RAT::DS::CalPMTs& calpmts = rEV.GetCalPMTs();
                for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++){
                  TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                  double pmt_r = pmtpos.Mag();
                  hrpmt.Fill(pmt_r);
                  hxpmt.Fill(pmtpos.X());
                  hypmt.Fill(pmtpos.Y());
                  hzpmt.Fill(pmtpos.Z());
                  htpmt.Fill((calpmts.GetPMT(ipmt)).GetTime());
                }
                      
                if(rEV.FitResultExists("waterFitter") && rEV.GetFitResult("waterFitter").GetValid()){
                    RAT::DS::FitVertex rvertex = rEV.GetFitResult("waterFitter").GetVertex(0);
                  if( rvertex.ContainsPosition() && rvertex.ValidPosition() ) {
                    hfitValid.Fill(1.);
                    double R = sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y() + rvertex.GetPosition().Z()*rvertex.GetPosition().Z());
                    //The below are all from fOptimiser inside PositionTimeLikelihood
                    if(rvertex.ValidPositivePositionError()) {
                      herrposx.Fill(rvertex.GetPosition().X(), rvertex.GetPositivePositionError().x());
                      herrposy.Fill(rvertex.GetPosition().Y(), rvertex.GetPositivePositionError().y());
                      herrposz.Fill(rvertex.GetPosition().Z(), rvertex.GetPositivePositionError().z());
                      herrposxnhits.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().x());
                      herrposxitr.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().x());
                      herrposynhits.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().y());
                      herrposyitr.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().y());
                      herrposznhits.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().z());
                      herrposzitr.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().z());
                    }
                    if(rvertex.ValidPositiveTimeError()){
                      herrtimex.Fill(rvertex.GetPosition().X(), rvertex.GetPositiveTimeError());
                      herrtimey.Fill(rvertex.GetPosition().Y(), rvertex.GetPositiveTimeError());
                      herrtimez.Fill(rvertex.GetPosition().Z(), rvertex.GetPositiveTimeError());
                    }
                    if(rvertex.ValidTime()){
                      htime.Fill(rvertex.GetTime());
                      htimeposx.Fill(rvertex.GetPosition().X(), rvertex.GetTime());
                      htimeposy.Fill(rvertex.GetPosition().Y(), rvertex.GetTime());
                      htimeposz.Fill(rvertex.GetPosition().Z(), rvertex.GetTime());
                    }
                    hposx.Fill(rvertex.GetPosition().X());
                    hposy.Fill(rvertex.GetPosition().Y());
                    hposz.Fill(rvertex.GetPosition().Z());
                    hposxy.Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                    hnhitsz.Fill(rEV.GetNhits(),rvertex.GetPosition().Z());
                    hposRz.Fill(R,rvertex.GetPosition().Z());
                    hposrhoz.Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X() + rvertex.GetPosition().Y()*rvertex.GetPosition().Y()),rvertex.GetPosition().Z());
                    hposR.Fill(R );
                    hposR3.Fill(pow(R,3)/pow(6005.3,3));
                    hitr.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"));
                    if(rvertex.ValidEnergy() && rvertex.ContainsEnergy()) {
                        henergy.Fill(rvertex.GetEnergy());
                        if(rvertex.ValidPositiveEnergyError()){
                          //Seems to be just filled with 1s for EnergyPromptLookup
                          herrenergy.Fill(rvertex.GetEnergy(), rvertex.GetPositiveEnergyError());
                        }
                    }
                    if(rEV.GetClassifierResult( "isotropy:waterFitter" ).GetValid()) hbeta14.Fill(rEV.GetClassifierResult("isotropy:waterFitter").GetClassification("snobeta14"));
                    
                    if(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR") > 0.55){  
                        double R = sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X()  + rvertex.GetPosition().Y()*rvertex.GetPosition().Y() + rvertex.GetPosition().Z()*rvertex.GetPosition().Z());
                        //The below are all from fOptimiser inside PositionTimeLikelihood
                        if(rvertex.ValidPositivePositionError()) {
                          herrposxgoodfit.Fill(rvertex.GetPosition().X(), rvertex.GetPositivePositionError().x());
                          herrposygoodfit.Fill(rvertex.GetPosition().Y(), rvertex.GetPositivePositionError().y());
                          herrposzgoodfit.Fill(rvertex.GetPosition().Z(), rvertex.GetPositivePositionError().z());
                          herrposxnhitsgoodfit.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().x());
                          herrposxitrgoodfit.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().x());
                          herrposynhitsgoodfit.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().y());
                          herrposyitrgoodfit.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().y());
                          herrposznhitsgoodfit.Fill(rEV.GetNhits(), rvertex.GetPositivePositionError().z());
                          herrposzitrgoodfit.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"), rvertex.GetPositivePositionError().z());
                        }
                        if(rvertex.ValidPositiveTimeError()){
                          herrtimexgoodfit.Fill(rvertex.GetPosition().X(), rvertex.GetPositiveTimeError());
                          herrtimeygoodfit.Fill(rvertex.GetPosition().Y(), rvertex.GetPositiveTimeError());
                          herrtimezgoodfit.Fill(rvertex.GetPosition().Z(), rvertex.GetPositiveTimeError());
                        }
                        if(rvertex.ValidTime()){
                          htime.Fill(rvertex.GetTime());
                          htimeposxgoodfit.Fill(rvertex.GetPosition().X(), rvertex.GetTime());
                          htimeposygoodfit.Fill(rvertex.GetPosition().Y(), rvertex.GetTime());
                          htimeposzgoodfit.Fill(rvertex.GetPosition().Z(), rvertex.GetTime());
                        }
                        hposxgoodfit.Fill(rvertex.GetPosition().X());
                        hposygoodfit.Fill(rvertex.GetPosition().Y());
                        hposzgoodfit.Fill(rvertex.GetPosition().Z());
                        hposxygoodfit.Fill(rvertex.GetPosition().X(),rvertex.GetPosition().Y());
                        hnhitszgoodfit.Fill(rEV.GetNhits(),rvertex.GetPosition().Z());
                        hposRzgoodfit.Fill(R,rvertex.GetPosition().Z());
                        hposrhozgoodfit.Fill(sqrt(rvertex.GetPosition().X()*rvertex.GetPosition().X() + rvertex.GetPosition().Y()*rvertex.GetPosition().Y()),rvertex.GetPosition().Z());
                        hposRgoodfit.Fill(R );
                        hposR3goodfit.Fill(pow(R,3)/pow(6005.3,3));
                        hitrgoodfit.Fill(rEV.GetClassifierResult("ITR:waterFitter").GetClassification("ITR"));
                        if(rvertex.ValidEnergy() && rvertex.ContainsEnergy()) {
                            henergy.Fill(rvertex.GetEnergy());
                            if(rvertex.ValidPositiveEnergyError()){
                              //Seems to be just filled with 1s for EnergyPromptLookup
                              herrenergygoodfit.Fill(rvertex.GetEnergy(), rvertex.GetPositiveEnergyError());
                            }
                        }
                        if(rEV.GetClassifierResult( "isotropy:waterFitter" ).GetValid()) hbeta14goodfit.Fill(rEV.GetClassifierResult("isotropy:waterFitter").GetClassification("snobeta14"));
                        RAT::DS::CalPMTs& calpmts = rEV.GetCalPMTs();
                        for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++){
                          TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                          double pmt_r = pmtpos.Mag();
                          hrpmtgoodfit.Fill(pmt_r);
                          hxpmtgoodfit.Fill(pmtpos.X());
                          hypmtgoodfit.Fill(pmtpos.Y());
                          hzpmtgoodfit.Fill(pmtpos.Z());
                          htpmtgoodfit.Fill((calpmts.GetPMT(ipmt)).GetTime());
                        }
                      } 
                  } else hfitValid.Fill(0.);
                } else hfitValid.Fill(0.);
              }
            }
          }
          start_days = start_time.GetDays();
          start_secs = start_time.GetSeconds();
          start_nsecs = start_time.GetNanoSeconds();
          end_days = end_time.GetDays();
          end_secs = end_time.GetSeconds();
          end_nsecs = end_time.GetNanoSeconds();
        }
    }
    run_duration = (end_days-start_days)*60*60*24 + (end_secs-start_secs) + ((end_nsecs-start_nsecs) * 1E-9);
    hduration.Fill(run_duration);
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
    plot_map["time"].SetHist(htime);
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
    plot_map["beta14"].SetHist(hbeta14);
    plot_map["energy"].SetHist(henergy);
    plot_map["posxy"].Set2DHist(hposxy);
    plot_map["timeposx"].Set2DHist(htimeposx);
    plot_map["timeposy"].Set2DHist(htimeposy);
    plot_map["timeposz"].Set2DHist(htimeposz);
    plot_map["errposx"].Set2DHist(herrposx);
    plot_map["errposy"].Set2DHist(herrposy);
    plot_map["errposz"].Set2DHist(herrposz);
    plot_map["errposxnhits"].Set2DHist(herrposxnhits);
    plot_map["errposxitr"].Set2DHist(herrposxitr);
    plot_map["errposynhits"].Set2DHist(herrposynhits);
    plot_map["errposyitr"].Set2DHist(herrposyitr);
    plot_map["errposznhits"].Set2DHist(herrposznhits);
    plot_map["errposzitr"].Set2DHist(herrposzitr);
    plot_map["errtimex"].Set2DHist(herrtimex);
    plot_map["errtimey"].Set2DHist(herrtimey);
    plot_map["errtimez"].Set2DHist(herrtimez);
    plot_map["errenergy"].Set2DHist(herrenergy);
    plot_map["nhitsz"].Set2DHist(hnhitsz);
    plot_map["posrhoz"].Set2DHist(hposrhoz);
    plot_map["posRz"].Set2DHist(hposRz);
    
    plot_map["itrgoodfit"].SetHist(hitrgoodfit);
    plot_map["timegoodfit"].SetHist(htimegoodfit);
    plot_map["posxgoodfit"].SetHist(hposxgoodfit);
    plot_map["posygoodfit"].SetHist(hposygoodfit);
    plot_map["poszgoodfit"].SetHist(hposzgoodfit);
    plot_map["posRgoodfit"].SetHist(hposRgoodfit);
    plot_map["posR3goodfit"].SetHist(hposR3goodfit);
    plot_map["rpmtgoodfit"].SetHist(hrpmtgoodfit);
    plot_map["xpmtgoodfit"].SetHist(hxpmtgoodfit);
    plot_map["ypmtgoodfit"].SetHist(hypmtgoodfit);
    plot_map["zpmtgoodfit"].SetHist(hzpmtgoodfit);
    plot_map["tpmtgoodfit"].SetHist(htpmtgoodfit);
    plot_map["beta14goodfit"].SetHist(hbeta14goodfit);
    plot_map["energygoodfit"].SetHist(henergygoodfit);
    plot_map["posxygoodfit"].Set2DHist(hposxygoodfit);
    plot_map["timeposxgoodfit"].Set2DHist(htimeposxgoodfit);
    plot_map["timeposygoodfit"].Set2DHist(htimeposygoodfit);
    plot_map["timeposzgoodfit"].Set2DHist(htimeposzgoodfit);
    plot_map["errposxgoodfit"].Set2DHist(herrposxgoodfit);
    plot_map["errposygoodfit"].Set2DHist(herrposygoodfit);
    plot_map["errposzgoodfit"].Set2DHist(herrposzgoodfit);
    plot_map["errposxnhitsgoodfit"].Set2DHist(herrposxnhitsgoodfit);
    plot_map["errposxitrgoodfit"].Set2DHist(herrposxitrgoodfit);
    plot_map["errposynhitsgoodfit"].Set2DHist(herrposynhitsgoodfit);
    plot_map["errposyitrgoodfit"].Set2DHist(herrposyitrgoodfit);
    plot_map["errposznhitsgoodfit"].Set2DHist(herrposznhitsgoodfit);
    plot_map["errposzitrgoodfit"].Set2DHist(herrposzitrgoodfit);
    plot_map["errtimexgoodfit"].Set2DHist(herrtimexgoodfit);
    plot_map["errtimeygoodfit"].Set2DHist(herrtimeygoodfit);
    plot_map["errtimezgoodfit"].Set2DHist(herrtimezgoodfit);
    plot_map["errenergygoodfit"].Set2DHist(herrenergygoodfit);
    plot_map["nhitszgoodfit"].Set2DHist(hnhitszgoodfit);
    plot_map["posrhozgoodfit"].Set2DHist(hposrhozgoodfit);
    plot_map["posRzgoodfit"].Set2DHist(hposRzgoodfit);
    
    plot_map["duration"].SetHist(hduration);
    plot_map["trigger"].SetHist(htrigger);
    plot_map["dataclean"].SetHist(hdataclean);
    plot_map["temp"].SetHist(htemp);
    plot_map["nhitstemp"].Set2DHist(hnhitstemp);
  
    //=============================================================================
    //Make some plots
   
    std::pair<std::string,std::string> run_info_subrun0 = ParseRunInfo(files[i][0]);
    std::vector<std::string> runs;
    std::vector<std::string> subruns;
    runs.push_back(run_info_subrun0.first);
    allruns.push_back(run_info_subrun0.first);
    for(unsigned k=0; k<files[i].size();k++){
       std::pair<std::string,std::string> run_info_subrunx = ParseRunInfo(files[i][k]);
       subruns.push_back(run_info_subrunx.second);
    }
    
    std::string outname = run_info_subrun0.first; 
    std::map<std::string, THPlot::THPlot>::iterator it;
    std::pair<std::vector<std::string>,std::vector<std::string> > run_info;
    run_info.first = runs;
    run_info.second = subruns;
    for ( it = plot_map.begin(); it != plot_map.end(); it++ )
    { 
      std::string label = it->first;  
      THPlot::THPlot plot = it->second;
      plot.SetOutFilename(output_dir + plot.GetOutFilename() + "_" + outname + postfix);
      plot.SetRunInfo(run_info);
      plot.GeneratePlot();
      TH1D totalhist = plot_map_totals[label].GetHist();
      TH2D totaltwoDhist = plot_map_totals[label].Get2DHist();
      TH1D temphist = plot.GetHist(); 
      TH2D temptwoDhist = plot.Get2DHist();
      totalhist.Add(&temphist);
      totaltwoDhist.Add(&temptwoDhist);
      plot_map_totals[label].SetHist(totalhist);
      plot_map_totals[label].Set2DHist(totaltwoDhist);
    }
    
       
    //Dont bother plotting any runs which last less than 1 minute in the comparison plots
    //if(run_duration < 60) continue;
    //Count the files/runs which are deemed good by whatever conditions we want to set
    file_count++;

    //Code for the vs run plots. For now dont use the class for this
    hnhits_vs_run->Fill(file_count-1, hnhits.GetMean());
    hnhits_vs_run->SetBinError(file_count,hnhits.GetMeanError());
    std::string bin_label = run_info_subrun0.first;
    if(file_count==1) start_run = bin_label;
    end_run = bin_label;
    hnhits_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    hTemp_vs_run->Fill(file_count-1, htemp.GetMean());
    hTemp_vs_run->SetBinError(file_count,htemp.GetMeanError());
    hTemp_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    htotalQ_vs_run->Fill(file_count-1, htotalQ.GetMean());
    htotalQ_vs_run->SetBinError(file_count,htotalQ.GetMeanError());
    htotalQ_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    hNCleanEvents_vs_run->Fill(file_count-1, n_cleanevents);
    hNCleanEvents_vs_run->SetBinError(file_count,sqrt(n_cleanevents));
    hNCleanEvents_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    hNEvents_vs_run->Fill(file_count-1, n_events);
    hNEvents_vs_run->SetBinError(file_count,sqrt(n_events));
    hNEvents_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    hNCleanEventsnorm_vs_run->Fill(file_count-1, n_cleanevents/run_duration);
    hNCleanEventsnorm_vs_run->SetBinError(file_count,sqrt(n_cleanevents/run_duration));
    hNCleanEventsnorm_vs_run->GetXaxis()->SetBinLabel(file_count,bin_label.c_str());
    //delete f;
  }
  
  //Create a version of all per-run plots which is a sum over all runs
  std::map<std::string, THPlot::THPlot>::iterator it;
  std::pair<std::vector<std::string>,std::vector<std::string> > run_info;
  run_info.first = allruns;
  run_info.second = allruns;
  for ( it = plot_map_totals.begin(); it != plot_map_totals.end(); it++ ){
    THPlot::THPlot plot = it->second;
    plot.SetRunInfo(run_info);
    plot.SetOutFilename(output_dir + plot.GetOutFilename() + "_" +start_run+"_to_"+end_run + postfix);
    plot.GeneratePlot();
  }
  
  SetStyle();
  TCanvas* c100 = new TCanvas("c100","c100",1000,400);
  hnhits_vs_run->SetMarkerColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_run->SetMarkerStyle(20);
  hnhits_vs_run->SetLineColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_run->GetYaxis()->SetTitle( "Mean #Hits" );
  hnhits_vs_run->GetXaxis()->SetTitle( "Run ID" );
  //Necessary rescale of axis to get rid of empty bins where we didnt fill because the run wasn't good
  hnhits_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hnhits_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hnhits_vs_run->Draw("PE1");
  hnhits_vs_run->SetMinimum(0);
  c100->SetLeftMargin(0.12);
  c100->SetBottomMargin(0.20);
  c100->Update();
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nhits_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hnhits_vs_time->SetMarkerColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_time->SetMarkerStyle(20);
  hnhits_vs_time->SetLineColor(TColor::GetColor(220, 24, 24));
  hnhits_vs_time->GetYaxis()->SetTitle( "Nhits/5 minutes" );
  hnhits_vs_time->GetXaxis()->SetTitle( "Time (hours)" );
  hnhits_vs_time->GetXaxis()->SetRangeUser(start_run_time/(60*60),end_run_time/(60*60));
  hnhits_vs_time->GetYaxis()->SetTitleOffset(0.8);
  hnhits_vs_time->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"nhits_vs_time_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nhits_vs_time_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hTemp_vs_run->SetMarkerColor(TColor::GetColor(250, 155, 107));
  hTemp_vs_run->SetMarkerStyle(20);
  hTemp_vs_run->SetLineColor(TColor::GetColor(250, 155, 107));
  hTemp_vs_run->GetYaxis()->SetTitle( "Mean cavity Temp (C)" );
  hTemp_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hTemp_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hTemp_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hTemp_vs_run->SetMinimum(0);
  hTemp_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"Temp_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"Temp_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  htotalQ_vs_run->SetMarkerColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_run->SetMarkerStyle(20);
  htotalQ_vs_run->SetLineColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_run->GetYaxis()->SetTitle( "Mean totalQ" );
  htotalQ_vs_run->GetXaxis()->SetTitle( "Run ID" );
  htotalQ_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  htotalQ_vs_run->GetYaxis()->SetTitleOffset(0.8);
  htotalQ_vs_run->SetMinimum(0);
  htotalQ_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"totalQ_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  htotalQ_vs_time->SetMarkerColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_time->SetMarkerStyle(20);
  htotalQ_vs_time->SetLineColor(TColor::GetColor(142, 24, 220));
  htotalQ_vs_time->GetYaxis()->SetTitle( "Total charge/5 minutes" );
  htotalQ_vs_time->GetXaxis()->SetTitle( "Time (hours)" );
  htotalQ_vs_time->GetXaxis()->SetRangeUser(start_run_time/(60*60),end_run_time/(60*60));
  htotalQ_vs_time->GetYaxis()->SetTitleOffset(0.8);
  htotalQ_vs_time->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"totalQ_vs_time_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"totalQ_vs_time_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNCleanEvents_vs_run->SetMarkerColor(kGreen);
  hNCleanEvents_vs_run->SetMarkerStyle(20);
  hNCleanEvents_vs_run->SetLineColor(kGreen);
  hNCleanEvents_vs_run->GetYaxis()->SetTitle( "# Clean Events" );
  hNCleanEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNCleanEvents_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hNCleanEvents_vs_run->SetMinimum(0);
  hNCleanEvents_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hNCleanEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"ncleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"ncleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNEvents_vs_run->SetMarkerColor(kGreen);
  hNEvents_vs_run->SetMarkerStyle(20);
  hNEvents_vs_run->SetLineColor(kGreen);
  hNEvents_vs_run->GetYaxis()->SetTitle( "# Events" );
  hNEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNEvents_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hNEvents_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hNEvents_vs_run->SetMinimum(0);
  hNEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNCleanEvents_vs_time->SetMarkerColor(kGreen+3);
  hNCleanEvents_vs_time->SetMarkerStyle(20);
  hNCleanEvents_vs_time->SetLineColor(kGreen+3);
  hNCleanEvents_vs_time->GetYaxis()->SetTitle( "# Clean Events/5 minutes" );
  hNCleanEvents_vs_time->GetXaxis()->SetTitle( "Time (hours)" );
  hNCleanEvents_vs_time->GetXaxis()->SetRangeUser(start_run_time/(60*60),end_run_time/(60*60));
  hNCleanEvents_vs_time->GetYaxis()->SetTitleOffset(0.8);
  hNCleanEvents_vs_time->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"ncleanevents_vs_time_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"ncleanevents_vs_time_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNEvents_vs_time->SetMarkerColor(kGreen);
  hNEvents_vs_time->SetMarkerStyle(20);
  hNEvents_vs_time->SetLineColor(kGreen);
  hNEvents_vs_time->GetYaxis()->SetTitle( "# Events/5 minutes" );
  hNEvents_vs_time->GetXaxis()->SetTitle( "Time (hours)" );
  hNEvents_vs_time->GetXaxis()->SetRangeUser(start_run_time/(60*60),end_run_time/(60*60));
  hNEvents_vs_time->GetYaxis()->SetTitleOffset(0.8);
  hNEvents_vs_time->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"nevents_vs_time_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"nevents_vs_time_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hNCleanEventsnorm_vs_run->SetMarkerColor(kGreen+3);
  hNCleanEventsnorm_vs_run->SetMarkerStyle(20);
  hNCleanEventsnorm_vs_run->SetLineColor(kGreen+3);
  hNCleanEventsnorm_vs_run->GetYaxis()->SetTitle( "# Clean events per second" );
  hNCleanEventsnorm_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hNCleanEventsnorm_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hNCleanEventsnorm_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hNCleanEventsnorm_vs_run->SetMinimum(0);
  hNCleanEventsnorm_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"ncleaneventsnorm_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"ncleaneventsnorm_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hFracCleanEvents_vs_run = (TH1D*)hNCleanEvents_vs_run->Clone("hFracCleanEvents_vs_run");
  hFracCleanEvents_vs_run->Divide(hNEvents_vs_run);
  hFracCleanEvents_vs_run->SetMarkerColor(kMagenta);
  hFracCleanEvents_vs_run->SetMarkerStyle(20);
  hFracCleanEvents_vs_run->SetLineColor(kMagenta);
  hFracCleanEvents_vs_run->GetYaxis()->SetTitle( "Fraction of clean Events" );
  hFracCleanEvents_vs_run->GetXaxis()->SetTitle( "Run ID" );
  hFracCleanEvents_vs_run->GetXaxis()->SetRangeUser(0,file_count);
  hFracCleanEvents_vs_run->GetYaxis()->SetTitleOffset(0.8);
  hFracCleanEvents_vs_run->SetMinimum(0);
  hFracCleanEvents_vs_run->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"fraccleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"fraccleanevents_vs_run_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  hFracCleanEvents_vs_time = (TH1D*)hNCleanEvents_vs_time->Clone("hFracCleanEvents_vs_time");
  hFracCleanEvents_vs_time->Divide(hNEvents_vs_time);
  hFracCleanEvents_vs_time->SetMarkerColor(kMagenta);
  hFracCleanEvents_vs_time->SetMarkerStyle(20);
  hFracCleanEvents_vs_time->SetLineColor(kMagenta);
  hFracCleanEvents_vs_time->GetYaxis()->SetTitle( "Fraction of clean Events/5 minutes" );
  hFracCleanEvents_vs_time->GetXaxis()->SetTitle( "Time (hours)" );
  hFracCleanEvents_vs_time->GetXaxis()->SetRangeUser(start_run_time/(60*60),end_run_time/(60*60));
  hFracCleanEvents_vs_time->GetYaxis()->SetTitleOffset(0.8);
  hFracCleanEvents_vs_time->Draw("PE1");
  c100->Update();
  c100->SaveAs((output_dir+"fraccleanevents_vs_time_"+start_run+"_to_"+end_run+postfix+".png").c_str());
  c100->SaveAs((output_dir+"fraccleanevents_vs_time_"+start_run+"_to_"+end_run+postfix+".pdf").c_str());
  
  delete c100;
  delete hnhits_vs_run;
  delete hnhits_vs_time;
  delete htotalQ_vs_run;
  delete htotalQ_vs_time;
  delete hNEvents_vs_run;
  delete hNEvents_vs_time;
  delete hNCleanEvents_vs_run;
  delete hNCleanEvents_vs_time;
  delete hFracCleanEvents_vs_run;
  delete hFracCleanEvents_vs_time;
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
  
  
  //Create a vector of runfile vectors - i.e. runfile is a vector of subrun files
  std::map<std::string, std::vector<std::string> > runmap;
  for(unsigned i=0; i<files.size(); i++){ 
    std::pair<std::string,std::string> run_info = ParseRunInfo(files[i]);
    std::string run = run_info.first;
    std::string subrun = run_info.second;
    runmap[run].push_back(files[i]);
  }
  //What is this monstrosity
  std::vector<std::vector<std::vector<std::string> > > runvecs; 
  std::vector<std::vector<std::string> > runvec;
  int exe_on = 20; 
  int count = 0; 
  //Split up the filelists into 10 runs at a time in order to make manageable webpages
  std::map<std::string, std::vector<std::string> >::iterator it;
  for ( it = runmap.begin(); it != runmap.end(); it++ ){
      std::vector<std::string> subrun_vec = it->second;
      runvec.push_back(subrun_vec);
      if(((count+1) % exe_on == 0) && count>0) {
        std::vector<std::vector<std::string>> runvec_temp = runvec;
        runvecs.push_back(runvec_temp);
        runvec.clear();
      }
      count++;
  }
  std::vector<std::vector<std::string> > runvec_temp = runvec;
  runvecs.push_back(runvec_temp);
  runvec.clear();
  
  for(unsigned j=0; j<runvecs.size(); j++){
    std::stringstream ss;
    ss << j;
    fs::create_directory(directory.c_str());
    std::string dirname = directory+ss.str()+"/";
    fs::create_directory(dirname.c_str());
    CreateRunPlots(runvecs[j],ntuple,dirname);
  }


  return 0;
}
